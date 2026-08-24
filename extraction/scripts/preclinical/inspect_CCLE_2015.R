suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(PharmacoGx)
})

# -------------------------------------------------------------------------
# Memory guard for macOS R vector heap
# -------------------------------------------------------------------------

if (exists("mem.maxVSize", mode = "function")) {
  try(mem.maxVSize(64000), silent = TRUE)
}

# -------------------------------------------------------------------------
# Config
# -------------------------------------------------------------------------

DATASET_ID <- 1

CCLE_2019_RDS_PATH <- "extraction/data/raw/preclinical/CCLE_2019.rds"
CCLE_2015_RDS_PATH <- "extraction/data/raw/preclinical/CCLE_2015.rds"

SHEET_CELL_LINE_QC_PATH <- "extraction/data/raw/preclinical/All_PSets_sarcoma_cell_line_QC.csv"

OUT_DIR <- "extraction/data/proc/preclinical/CCLE"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_CELL_LINES_RAW <- paste0(
  "A-204|D-247MG|NCI-H2373|NCI-H2596|RKN|S-117|SW684|VA-ES-BJ|",
  "A-204|CAL-78|CS-1 [Human chondrosarcoma]|Hs 633.T|OUMS-27|Rh41|",
  "S-117|SNU-1077|SYO-1|TE 125.T|VA-ES-BJ|CHSA0011|D-247MG|GI-1|",
  "KYM-1|NCI-H2373|RS-5|Rh41|SK-LMS-1|TE 441.T|TE 617.T|",
  "VA-ES-BJ|Aska-SS"
)

TARGET_CELL_LINES <- unique(trimws(
  unlist(strsplit(TARGET_CELL_LINES_RAW, "\\|"))
))

# -------------------------------------------------------------------------
# General helpers
# -------------------------------------------------------------------------

read_updated_rds <- function(path) {
  obj <- readRDS(path)

  obj <- tryCatch(
    updateObject(obj, verbose = TRUE),
    error = function(e) {
      warning("updateObject failed for ", path, ": ", conditionMessage(e))
      obj
    }
  )

  obj
}

clean_na <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | x == "NS" | is.na(x)] <- NA_character_
  x
}

parse_age_int <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | x == "NS" | x == "Adult" | is.na(x)] <- NA_character_
  suppressWarnings(as.integer(gsub("[^0-9]", "", x)))
}

strip_ensembl_version <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | is.na(x)] <- NA_character_
  sub("\\.[0-9]+$", "", x)
}

normalize_cell_line_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  # Handles values like "CS-1 [Human chondrosarcoma]" -> "CS-1"
  x <- gsub("\\s*\\[.*\\]\\s*$", "", x)

  # Makes A-204, A204, A 204 comparable.
  x <- toupper(gsub("[^A-Za-z0-9]", "", x))
  x[x == "" | is.na(x)] <- NA_character_

  x
}

normalize_category_value <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x <- gsub("[ -]+", "_", x)
  x[x == "" | x == "na" | is.na(x)] <- NA_character_
  x
}

normalize_array_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- basename(x)
  x <- sub("\\.CEL(\\.GZ)?$", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "" | x == "NA" | is.na(x)] <- NA_character_
  x
}

has_col <- function(dt, col) {
  col %in% colnames(dt)
}

first_existing_col <- function(dt, candidates) {
  for (candidate in candidates) {
    if (candidate %in% colnames(dt)) {
      return(candidate)
    }
  }

  NA_character_
}

safe_col <- function(dt, candidates, default = NA_character_) {
  col <- first_existing_col(dt, candidates)

  if (is.na(col)) {
    return(rep(default, nrow(dt)))
  }

  dt[[col]]
}

coalesce_dt_cols <- function(dt, candidates) {
  out <- rep(NA_character_, nrow(dt))

  for (col in candidates) {
    if (has_col(dt, col)) {
      vals <- clean_na(dt[[col]])
      fill_idx <- is.na(out) & !is.na(vals)
      out[fill_idx] <- vals[fill_idx]
    }
  }

  out
}

get_rowdata_col <- function(se, candidates) {
  rd <- as.data.frame(rowData(se))

  for (candidate in candidates) {
    if (candidate %in% colnames(rd)) {
      return(as.character(rd[[candidate]]))
    }
  }

  rep(NA_character_, nrow(se))
}

make_gene_mapping <- function(se, ensembl_candidates, name_candidates) {
  feature_id <- rownames(se)

  gene_id <- get_rowdata_col(se, ensembl_candidates)
  gene_name <- get_rowdata_col(se, name_candidates)

  gene_id <- strip_ensembl_version(gene_id)
  gene_name <- clean_na(gene_name)

  feature_as_ensembl <- strip_ensembl_version(feature_id)

  use_feature_id <- (
    is.na(gene_id) |
      gene_id == ""
  ) & grepl("^ENSG[0-9]+$", feature_as_ensembl)

  gene_id[use_feature_id] <- feature_as_ensembl[use_feature_id]

  data.table(
    feature_id = feature_id,
    gene_id = gene_id,
    gene_name = gene_name
  )
}

# -------------------------------------------------------------------------
# Uploaded QC-sheet helpers
# -------------------------------------------------------------------------

build_sheet_alias_table <- function(sheet_dt, row_subset = NULL) {
  if (is.null(row_subset)) {
    sheet_sub <- copy(sheet_dt)
  } else {
    sheet_sub <- copy(sheet_dt[row_subset])
  }

  if (!has_col(sheet_sub, "sheet_row_id")) {
    sheet_sub[, sheet_row_id := .I]
  }

  alias_cols <- c(
    "cell_line",
    "sampleid",
    "CCLE_rnaseq.sampleid",
    "Cell.line.primary.name",
    "CCLE.name"
  )

  alias_cols <- alias_cols[alias_cols %in% colnames(sheet_sub)]

  if (length(alias_cols) == 0) {
    stop(
      "No usable cell-line alias columns found in sheet. Expected one of: ",
      paste(
        c(
          "cell_line",
          "sampleid",
          "CCLE_rnaseq.sampleid",
          "Cell.line.primary.name",
          "CCLE.name"
        ),
        collapse = ", "
      )
    )
  }

  alias_pieces <- list()

  for (col in alias_cols) {
    alias_pieces[[col]] <- data.table(
      sheet_row_id = sheet_sub$sheet_row_id,
      dataset = if (has_col(sheet_sub, "dataset")) clean_na(sheet_sub[["dataset"]]) else NA_character_,
      object_type = if (has_col(sheet_sub, "object_type")) clean_na(sheet_sub[["object_type"]]) else NA_character_,

      alias_column = col,
      alias_value = clean_na(sheet_sub[[col]]),
      normalized_alias_value = normalize_cell_line_name(sheet_sub[[col]]),

      sheet_cell_line_name = if (has_col(sheet_sub, "cell_line")) clean_na(sheet_sub[["cell_line"]]) else NA_character_,
      sheet_sampleid = if (has_col(sheet_sub, "sampleid")) clean_na(sheet_sub[["sampleid"]]) else NA_character_,
      mod_tissueid = if (has_col(sheet_sub, "mod_tissueid")) clean_na(sheet_sub[["mod_tissueid"]]) else NA_character_,

      accession = clean_na(coalesce_dt_cols(
        sheet_sub,
        c(
          "Cellosaurus.Accession.id",
          "cellosaurus.cvcl_id",
          "cellosaurus.cellosaurus.cvcl_id"
        )
      )),

      category = clean_na(coalesce_dt_cols(
        sheet_sub,
        c(
          "cellosaurus.category",
          "cellosaurus.cellosaurus.category",
          "CellLine.Type"
        )
      )),

      sex = clean_na(coalesce_dt_cols(
        sheet_sub,
        c(
          "cellosaurus.sex",
          "cellosaurus.cellosaurus.sex",
          "Gender",
          "sex"
        )
      )),

      age = clean_na(coalesce_dt_cols(
        sheet_sub,
        c(
          "cellosaurus.age",
          "cellosaurus.cellosaurus.age",
          "Age",
          "age"
        )
      ))
    )
  }

  alias_dt <- rbindlist(alias_pieces, fill = TRUE)

  alias_dt <- alias_dt[
    !is.na(alias_value) &
      !is.na(normalized_alias_value)
  ]

  alias_dt[
    ,
    metadata_priority := fifelse(
      dataset == "CCLE_2019",
      1L,
      fifelse(dataset == "CCLE_2015", 2L, 99L)
    )
  ]

  alias_dt[
    ,
    metadata_score :=
      as.integer(!is.na(sheet_cell_line_name)) +
      as.integer(!is.na(accession)) +
      as.integer(!is.na(category)) +
      as.integer(!is.na(sex)) +
      as.integer(!is.na(age))
  ]

  unique(
    alias_dt,
    by = c("sheet_row_id", "alias_column", "normalized_alias_value")
  )
}

read_sheet_soft_tissue_targets <- function(path, out_dir) {
  if (!file.exists(path)) {
    warning(
      "Sheet file not found at: ",
      path,
      ". Sheet-derived soft tissue targets and cell-line metadata will be skipped."
    )

    return(list(
      sheet_dt = data.table(),
      sheet_soft_tissue_dt = data.table(),
      sheet_soft_tissue_alias_dt = data.table(),
      sheet_metadata_alias_dt = data.table()
    ))
  }

  sheet_dt <- fread(path)
  sheet_dt[, sheet_row_id := .I]

  fwrite(
    data.table(column_name = colnames(sheet_dt)),
    file.path(out_dir, "ccle_sheet_column_names.csv")
  )

  if (!has_col(sheet_dt, "mod_tissueid")) {
    stop(
      "Uploaded QC sheet must contain column 'mod_tissueid'. ",
      "See: ",
      file.path(out_dir, "ccle_sheet_column_names.csv")
    )
  }

  sheet_soft_tissue_dt <- sheet_dt[
    normalize_category_value(mod_tissueid) == "soft_tissue"
  ]

  sheet_soft_tissue_alias_dt <- build_sheet_alias_table(
    sheet_dt = sheet_soft_tissue_dt
  )

  sheet_metadata_alias_dt <- build_sheet_alias_table(
    sheet_dt = sheet_dt
  )

  fwrite(
    sheet_soft_tissue_dt,
    file.path(out_dir, "ccle_sheet_soft_tissue_rows.csv")
  )

  fwrite(
    sheet_soft_tissue_alias_dt,
    file.path(out_dir, "ccle_sheet_soft_tissue_aliases.csv")
  )

  fwrite(
    sheet_metadata_alias_dt,
    file.path(out_dir, "ccle_sheet_cell_line_metadata_aliases.csv")
  )

  cat("QC-sheet soft tissue rows:", nrow(sheet_soft_tissue_dt), "\n")
  cat("QC-sheet soft tissue aliases:", nrow(sheet_soft_tissue_alias_dt), "\n")
  cat("QC-sheet metadata aliases:", nrow(sheet_metadata_alias_dt), "\n")

  list(
    sheet_dt = sheet_dt,
    sheet_soft_tissue_dt = sheet_soft_tissue_dt,
    sheet_soft_tissue_alias_dt = sheet_soft_tissue_alias_dt,
    sheet_metadata_alias_dt = sheet_metadata_alias_dt
  )
}

attach_sheet_cell_line_metadata <- function(selected_sample_dt, sheet_metadata_alias_dt) {
  selected <- copy(selected_sample_dt)
  selected[, selected_row_id := .I]

  if (nrow(sheet_metadata_alias_dt) == 0) {
    selected[, final_cell_line_name := clean_na(`cellosaurus.cellLineName`)]
    selected[is.na(final_cell_line_name), final_cell_line_name := clean_na(sampleid)]

    selected[, sheet_metadata_match := FALSE]
    selected[, sheet_metadata_alias_column := NA_character_]
    selected[, sheet_metadata_alias_value := NA_character_]
    selected[, sheet_cell_line_name := NA_character_]
    selected[, sheet_accession := NA_character_]
    selected[, sheet_category := NA_character_]
    selected[, sheet_sex := NA_character_]
    selected[, sheet_age := NA_character_]
    selected[, sheet_metadata_dataset := NA_character_]
    selected[, sheet_metadata_object_type := NA_character_]

    return(selected)
  }

  selected_alias_cols <- c(
    "cellosaurus.cellLineName",
    "CCLE.name",
    "CCLE.sampleid",
    "sampleid",
    "unique.sampleid",
    "sample_rowname"
  )

  selected_alias_cols <- selected_alias_cols[
    selected_alias_cols %in% colnames(selected)
  ]

  selected_alias_pieces <- list()

  for (col in selected_alias_cols) {
    selected_alias_pieces[[col]] <- data.table(
      selected_row_id = selected$selected_row_id,
      selected_alias_column = col,
      selected_alias_value = clean_na(selected[[col]]),
      normalized_alias_value = normalize_cell_line_name(selected[[col]])
    )
  }

  selected_alias_dt <- rbindlist(selected_alias_pieces, fill = TRUE)

  selected_alias_dt <- selected_alias_dt[
    !is.na(selected_alias_value) &
      !is.na(normalized_alias_value)
  ]

  metadata_joined <- merge(
    selected_alias_dt,
    sheet_metadata_alias_dt,
    by = "normalized_alias_value",
    all = FALSE,
    allow.cartesian = TRUE
  )

  if (nrow(metadata_joined) > 0) {
    setorder(
      metadata_joined,
      selected_row_id,
      metadata_priority,
      -metadata_score
    )

    best_metadata <- metadata_joined[
      ,
      .SD[1],
      by = selected_row_id
    ]

    best_metadata <- best_metadata[
      ,
      .(
        selected_row_id,
        sheet_metadata_match = TRUE,
        sheet_metadata_alias_column = alias_column,
        sheet_metadata_alias_value = alias_value,
        sheet_cell_line_name,
        sheet_accession = accession,
        sheet_category = category,
        sheet_sex = sex,
        sheet_age = age,
        sheet_metadata_dataset = dataset,
        sheet_metadata_object_type = object_type
      )
    ]

    selected <- merge(
      selected,
      best_metadata,
      by = "selected_row_id",
      all.x = TRUE
    )
  } else {
    selected[, sheet_metadata_match := FALSE]
    selected[, sheet_metadata_alias_column := NA_character_]
    selected[, sheet_metadata_alias_value := NA_character_]
    selected[, sheet_cell_line_name := NA_character_]
    selected[, sheet_accession := NA_character_]
    selected[, sheet_category := NA_character_]
    selected[, sheet_sex := NA_character_]
    selected[, sheet_age := NA_character_]
    selected[, sheet_metadata_dataset := NA_character_]
    selected[, sheet_metadata_object_type := NA_character_]
  }

  selected[is.na(sheet_metadata_match), sheet_metadata_match := FALSE]

  selected[, final_cell_line_name := clean_na(`cellosaurus.cellLineName`)]
  selected[!is.na(sheet_cell_line_name), final_cell_line_name := sheet_cell_line_name]
  selected[is.na(final_cell_line_name), final_cell_line_name := clean_na(sampleid)]

  selected
}

write_cell_line_selection_audits <- function(
  selected_sample_dt,
  sheet_soft_tissue_alias_dt,
  out_dir
) {
  final_cell_line_list <- unique(
    data.table(
      sampleid = clean_na(selected_sample_dt[["sampleid"]]),
      cell_line_name = clean_na(selected_sample_dt[["final_cell_line_name"]]),
      ccle_cell_line_name = clean_na(selected_sample_dt[["cellosaurus.cellLineName"]]),
      ccle_name = clean_na(selected_sample_dt[["CCLE.name"]]),
      site_primary = clean_na(selected_sample_dt[["CCLE.site_Primary"]]),
      disease_type = clean_na(selected_sample_dt[["CCLE.type"]]),

      original_target_match = selected_sample_dt$target_match,
      original_soft_tissue_match = selected_sample_dt$soft_tissue_match,
      sheet_soft_tissue_match = selected_sample_dt$sheet_soft_tissue_match,

      sheet_metadata_match = selected_sample_dt$sheet_metadata_match,
      sheet_metadata_dataset = selected_sample_dt$sheet_metadata_dataset,
      sheet_metadata_alias_column = selected_sample_dt$sheet_metadata_alias_column,
      sheet_metadata_alias_value = selected_sample_dt$sheet_metadata_alias_value
    ),
    by = "sampleid"
  )

  final_cell_line_list[, normalized_cell_line_name := normalize_cell_line_name(cell_line_name)]
  final_cell_line_list[, normalized_sampleid := normalize_cell_line_name(sampleid)]

  final_cell_line_list[
    ,
    selection_criteria := paste(
      c(
        if (isTRUE(original_target_match)) "original_target_list" else NA_character_,
        if (isTRUE(original_soft_tissue_match)) "ccle_2019_soft_tissue_filter" else NA_character_,
        if (isTRUE(sheet_soft_tissue_match)) "qc_sheet_mod_tissueid_soft_tissue" else NA_character_
      )[!is.na(c(
        if (isTRUE(original_target_match)) "original_target_list" else NA_character_,
        if (isTRUE(original_soft_tissue_match)) "ccle_2019_soft_tissue_filter" else NA_character_,
        if (isTRUE(sheet_soft_tissue_match)) "qc_sheet_mod_tissueid_soft_tissue" else NA_character_
      ))],
      collapse = "|"
    ),
    by = sampleid
  ]

  fwrite(
    final_cell_line_list,
    file.path(out_dir, "ccle_final_extracted_cell_line_list.csv")
  )

  original_criteria_list <- final_cell_line_list[
    original_target_match == TRUE |
      original_soft_tissue_match == TRUE
  ]

  fwrite(
    original_criteria_list,
    file.path(out_dir, "ccle_original_criteria_extracted_cell_line_list.csv")
  )

  sheet_added_list <- final_cell_line_list[
    sheet_soft_tissue_match == TRUE &
      original_target_match == FALSE &
      original_soft_tissue_match == FALSE
  ]

  fwrite(
    sheet_added_list,
    file.path(out_dir, "ccle_sheet_soft_tissue_added_to_original_criteria.csv")
  )

  missing_sheet_metadata <- final_cell_line_list[
    sheet_metadata_match == FALSE |
      is.na(sheet_metadata_match)
  ]

  fwrite(
    missing_sheet_metadata,
    file.path(out_dir, "ccle_cell_line_metadata_missing_from_sheet.csv")
  )

  if (nrow(sheet_soft_tissue_alias_dt) > 0) {
    final_norm_values <- unique(c(
      final_cell_line_list$normalized_cell_line_name,
      final_cell_line_list$normalized_sampleid
    ))

    sheet_not_extracted <- sheet_soft_tissue_alias_dt[
      !(normalized_alias_value %in% final_norm_values)
    ]

    sheet_not_extracted <- unique(
      sheet_not_extracted,
      by = c(
        "dataset",
        "sheet_cell_line_name",
        "sheet_sampleid",
        "normalized_alias_value"
      )
    )

    fwrite(
      sheet_not_extracted,
      file.path(out_dir, "ccle_sheet_soft_tissue_not_extracted.csv")
    )
  }

  cat("Wrote cell-line selection audit CSVs\n")
}

# -------------------------------------------------------------------------
# CCLE 2019 selection helpers
# -------------------------------------------------------------------------

build_selected_sample_table_2019 <- function(
  pset_2019,
  target_cell_lines,
  sheet_soft_tissue_alias_dt = data.table()
) {
  sample_dt <- as.data.table(pset_2019@sample, keep.rownames = "sample_rowname")

  target_norm <- normalize_cell_line_name(target_cell_lines)

  candidate_cols <- c(
    "cellosaurus.cellLineName",
    "CCLE.name",
    "CCLE.sampleid",
    "sampleid",
    "unique.sampleid",
    "sample_rowname"
  )

  sample_dt[, target_match := FALSE]

  for (col in candidate_cols) {
    if (has_col(sample_dt, col)) {
      sample_dt[
        normalize_cell_line_name(get(col)) %in% target_norm,
        target_match := TRUE
      ]
    }
  }

  sample_dt[, soft_tissue_match := FALSE]

  if (has_col(sample_dt, "CCLE.site_Primary")) {
    sample_dt[
      normalize_category_value(`CCLE.site_Primary`) == "soft_tissue",
      soft_tissue_match := TRUE
    ]
  }

  if (has_col(sample_dt, "CCLE.type")) {
    sample_dt[
      normalize_category_value(`CCLE.type`) == "soft_tissue",
      soft_tissue_match := TRUE
    ]
  }

  sample_dt[, sheet_soft_tissue_match := FALSE]

  if (nrow(sheet_soft_tissue_alias_dt) > 0) {
    sheet_target_norm <- unique(sheet_soft_tissue_alias_dt$normalized_alias_value)

    for (col in candidate_cols) {
      if (has_col(sample_dt, col)) {
        sample_dt[
          normalize_cell_line_name(get(col)) %in% sheet_target_norm,
          sheet_soft_tissue_match := TRUE
        ]
      }
    }
  }

  selected <- sample_dt[
    target_match == TRUE |
      soft_tissue_match == TRUE |
      sheet_soft_tissue_match == TRUE
  ]

  if (nrow(selected) == 0) {
    stop(
      "No matching CCLE 2019 samples found for original target list, ",
      "CCLE soft_tissue filter, or QC-sheet mod_tissueid == 'Soft Tissue'."
    )
  }

  selected
}

build_canonical_lookup <- function(selected_sample_dt) {
  canonical_sample_id <- clean_na(selected_sample_dt[["sampleid"]])
  canonical_cell_line <- clean_na(selected_sample_dt[["final_cell_line_name"]])

  alt_cols <- c(
    "cellosaurus.cellLineName",
    "final_cell_line_name",
    "CCLE.name",
    "CCLE.sampleid",
    "sampleid",
    "unique.sampleid",
    "sample_rowname"
  )

  pieces <- list()

  for (col in alt_cols) {
    if (has_col(selected_sample_dt, col)) {
      pieces[[col]] <- data.table(
        lookup_name = clean_na(selected_sample_dt[[col]]),
        normalized_lookup_name = normalize_cell_line_name(selected_sample_dt[[col]]),
        sample_id = canonical_sample_id,
        cell_line_name = canonical_cell_line,
        source_column = col
      )
    }
  }

  lookup_dt <- rbindlist(pieces, fill = TRUE)

  lookup_dt <- lookup_dt[
    !is.na(normalized_lookup_name) &
      !is.na(sample_id) &
      !is.na(cell_line_name)
  ]

  unique(
    lookup_dt,
    by = c("normalized_lookup_name", "sample_id", "source_column")
  )
}

# -------------------------------------------------------------------------
# Generic molecular profile column mapping
# -------------------------------------------------------------------------

build_profile_column_map <- function(se, canonical_lookup) {
  cd <- as.data.table(
    as.data.frame(colData(se)),
    keep.rownames = "coldata_rownames"
  )

  cd[, colname := colnames(se)]

  candidate_cols <- c(
    "cellosaurus.cellLineName",
    "CCLE.name",
    "CCLE.sampleid",
    "sampleid",
    "unique.sampleid",
    "dataset_sample_id",
    "ccle_sample_id",
    "Cell.line.primary.name",
    "Cell_Line",
    "Name",
    "samplename",
    "filename",
    "Expression.arrays",
    "SNP.arrays",
    "rownames",
    "coldata_rownames",
    "colname"
  )

  pieces <- list()

  for (col in candidate_cols) {
    if (has_col(cd, col)) {
      pieces[[col]] <- data.table(
        colname = cd$colname,
        matched_column = col,
        normalized_lookup_name = normalize_cell_line_name(cd[[col]])
      )
    }
  }

  if (length(pieces) == 0) {
    return(data.table(
      colname = character(),
      sample_id = character(),
      cell_line_name = character(),
      matched_column = character()
    ))
  }

  candidate_map <- rbindlist(pieces, fill = TRUE)
  candidate_map <- candidate_map[!is.na(normalized_lookup_name)]

  canonical_lookup_unique <- unique(
    canonical_lookup[, .(normalized_lookup_name, sample_id, cell_line_name)],
    by = "normalized_lookup_name"
  )

  mapped <- merge(
    candidate_map,
    canonical_lookup_unique,
    by = "normalized_lookup_name",
    all = FALSE
  )

  mapped <- mapped[colname %in% colnames(se)]

  mapped <- unique(mapped, by = "colname")

  mapped[, .(colname, sample_id, cell_line_name, matched_column)]
}

# -------------------------------------------------------------------------
# CCLE 2015 RNA/microarray-specific column mapping
# -------------------------------------------------------------------------

build_2015_rna_column_map <- function(rna_se, pset_2015, canonical_lookup) {
  rna_cd <- as.data.table(
    as.data.frame(colData(rna_se)),
    keep.rownames = "coldata_rownames"
  )

  rna_cd[, colname := colnames(rna_se)]

  rna_direct_cols <- c(
    "Cell.line.primary.name",
    "sampleid",
    "CCLE.name",
    "samplename",
    "colname",
    "coldata_rownames",
    "rownames"
  )

  rna_direct_pieces <- list()

  for (col in rna_direct_cols) {
    if (has_col(rna_cd, col)) {
      rna_direct_pieces[[col]] <- data.table(
        colname = rna_cd$colname,
        rna_matched_column = col,
        normalized_lookup_name = normalize_cell_line_name(rna_cd[[col]])
      )
    }
  }

  canonical_lookup_unique <- unique(
    canonical_lookup[, .(normalized_lookup_name, sample_id, cell_line_name)],
    by = "normalized_lookup_name"
  )

  if (length(rna_direct_pieces) > 0) {
    rna_direct_candidates <- rbindlist(rna_direct_pieces, fill = TRUE)
    rna_direct_candidates <- rna_direct_candidates[!is.na(normalized_lookup_name)]

    direct_mapped <- merge(
      rna_direct_candidates,
      canonical_lookup_unique,
      by = "normalized_lookup_name",
      all = FALSE
    )

    if (nrow(direct_mapped) > 0) {
      direct_mapped[, matched_column := paste0("rna_colData.", rna_matched_column)]
    }
  } else {
    direct_mapped <- data.table()
  }

  if (!exists("direct_mapped") || nrow(direct_mapped) == 0) {
    direct_mapped <- data.table(
      colname = character(),
      sample_id = character(),
      cell_line_name = character(),
      matched_column = character()
    )
  } else {
    direct_mapped <- direct_mapped[
      ,
      .(colname, sample_id, cell_line_name, matched_column)
    ]
  }

  rna_array_cols <- c(
    "colname",
    "coldata_rownames",
    "rownames",
    "filename",
    "samplename",
    "Expression.arrays"
  )

  rna_array_pieces <- list()

  for (col in rna_array_cols) {
    if (has_col(rna_cd, col)) {
      rna_array_pieces[[col]] <- data.table(
        colname = rna_cd$colname,
        rna_matched_column = col,
        normalized_array_id = normalize_array_id(rna_cd[[col]])
      )
    }
  }

  if (length(rna_array_pieces) == 0) {
    via_sample_mapped <- data.table(
      colname = character(),
      sample_id = character(),
      cell_line_name = character(),
      matched_column = character()
    )
  } else {
    rna_array_candidates <- rbindlist(rna_array_pieces, fill = TRUE)
    rna_array_candidates <- rna_array_candidates[!is.na(normalized_array_id)]

    sample_dt <- as.data.table(
      pset_2015@sample,
      keep.rownames = "sample_rowname"
    )

    sample_cell_line_col <- first_existing_col(
      sample_dt,
      c(
        "Cell.line.primary.name",
        "cellosaurus.cellLineName",
        "CCLE.name",
        "sampleid"
      )
    )

    if (is.na(sample_cell_line_col)) {
      stop(
        "Could not find a cell-line-name column in CCLE_2015@sample. ",
        "Available columns are: ",
        paste(colnames(sample_dt), collapse = ", ")
      )
    }

    sample_array_cols <- c(
      "Expression.arrays",
      "sample_rowname",
      "rownames",
      "filename",
      "samplename"
    )

    sample_array_pieces <- list()

    for (col in sample_array_cols) {
      if (has_col(sample_dt, col)) {
        sample_array_pieces[[col]] <- data.table(
          sample_matched_column = col,
          normalized_array_id = normalize_array_id(sample_dt[[col]]),
          ccle_2015_cell_line_name = clean_na(sample_dt[[sample_cell_line_col]])
        )
      }
    }

    if (length(sample_array_pieces) == 0) {
      via_sample_mapped <- data.table(
        colname = character(),
        sample_id = character(),
        cell_line_name = character(),
        matched_column = character()
      )
    } else {
      sample_array_lookup <- rbindlist(sample_array_pieces, fill = TRUE)

      sample_array_lookup <- sample_array_lookup[
        !is.na(normalized_array_id) &
          !is.na(ccle_2015_cell_line_name)
      ]

      via_sample <- merge(
        rna_array_candidates,
        sample_array_lookup,
        by = "normalized_array_id",
        all = FALSE,
        allow.cartesian = TRUE
      )

      via_sample[
        ,
        normalized_lookup_name := normalize_cell_line_name(ccle_2015_cell_line_name)
      ]

      via_sample_mapped <- merge(
        via_sample,
        canonical_lookup_unique,
        by = "normalized_lookup_name",
        all = FALSE
      )

      if (nrow(via_sample_mapped) > 0) {
        via_sample_mapped[
          ,
          matched_column := paste0(
            "rna_colData.",
            rna_matched_column,
            " -> CCLE_2015@sample.",
            sample_matched_column,
            " -> ",
            sample_cell_line_col
          )
        ]

        via_sample_mapped <- via_sample_mapped[
          ,
          .(colname, sample_id, cell_line_name, matched_column)
        ]
      } else {
        via_sample_mapped <- data.table(
          colname = character(),
          sample_id = character(),
          cell_line_name = character(),
          matched_column = character()
        )
      }
    }
  }

  mapped <- rbindlist(
    list(
      direct_mapped,
      via_sample_mapped
    ),
    fill = TRUE
  )

  mapped <- mapped[
    !is.na(colname) &
      !is.na(sample_id) &
      !is.na(cell_line_name)
  ]

  mapped <- mapped[colname %in% colnames(rna_se)]

  duplicate_sample_map <- mapped[, .N, by = sample_id][N > 1]

  if (nrow(duplicate_sample_map) > 0) {
    cat(
      "CCLE 2015 RNA column map has duplicate source columns for ",
      nrow(duplicate_sample_map),
      " canonical samples. Keeping first per sample_id.\n",
      sep = ""
    )
  }

  mapped <- unique(mapped, by = "sample_id")

  mapped[, .(colname, sample_id, cell_line_name, matched_column)]
}

# -------------------------------------------------------------------------
# Chunked long-assay writer
# -------------------------------------------------------------------------

write_long_assay_with_column_map <- function(
  se,
  gene_map,
  column_map,
  out_path,
  value_col,
  assay_name = "exprs",
  value_as_character = FALSE,
  column_chunk_size = 5
) {
  if (nrow(column_map) == 0) {
    empty_dt <- data.table(
      sample_id = character(),
      gene_id = character()
    )
    empty_dt[[value_col]] <- if (value_as_character) character() else numeric()
    fwrite(empty_dt, out_path)
    return(invisible(NULL))
  }

  if (!(assay_name %in% assayNames(se))) {
    assay_name <- assayNames(se)[1]
  }

  gene_map <- copy(gene_map)
  gene_map[, gene_id := strip_ensembl_version(gene_id)]

  gene_map <- gene_map[
    !is.na(feature_id) &
      !is.na(gene_id) &
      gene_id != ""
  ]

  duplicate_samples <- column_map[, .N, by = sample_id][N > 1]

  if (nrow(duplicate_samples) > 0) {
    duplicate_path <- sub("\\.csv$", "_duplicate_sample_columns.csv", out_path)

    fwrite(
      duplicate_samples,
      duplicate_path
    )

    warning(
      "Some source columns mapped to the same sample_id for ",
      out_path,
      ". Keeping the first source column per sample_id. ",
      "See: ",
      duplicate_path
    )
  }

  column_map <- unique(column_map, by = "sample_id")
  column_map <- column_map[colname %in% colnames(se)]

  selected_cols <- column_map$colname

  if (length(selected_cols) == 0) {
    empty_dt <- data.table(
      sample_id = character(),
      gene_id = character()
    )
    empty_dt[[value_col]] <- if (value_as_character) character() else numeric()
    fwrite(empty_dt, out_path)
    return(invisible(NULL))
  }

  if (file.exists(out_path)) {
    file.remove(out_path)
  }

  first_write <- TRUE

  for (start_idx in seq(1, length(selected_cols), by = column_chunk_size)) {
    end_idx <- min(start_idx + column_chunk_size - 1, length(selected_cols))
    chunk_cols <- selected_cols[start_idx:end_idx]

    cat(
      "Writing ",
      basename(out_path),
      " columns ",
      start_idx,
      "-",
      end_idx,
      " of ",
      length(selected_cols),
      "\n",
      sep = ""
    )

    chunk_map <- column_map[colname %in% chunk_cols]

    mat_chunk <- assay(se, assay_name)[, chunk_cols, drop = FALSE]

    dt <- as.data.table(as.table(mat_chunk))
    setnames(dt, c("feature_id", "source_colname", value_col))

    dt[, feature_id := as.character(feature_id)]
    dt[, source_colname := as.character(source_colname)]

    dt <- merge(
      dt,
      chunk_map[, .(source_colname = colname, sample_id)],
      by = "source_colname",
      all.x = TRUE
    )

    dt <- merge(
      dt,
      gene_map[, .(feature_id, gene_id)],
      by = "feature_id",
      all.x = TRUE
    )

    dt[, source_colname := NULL]
    dt[, feature_id := NULL]
    dt[, gene_id := strip_ensembl_version(gene_id)]

    dt <- dt[
      !is.na(sample_id) &
        !is.na(gene_id) &
        gene_id != "" &
        !is.na(get(value_col))
    ]

    if (value_as_character) {
      dt[, (value_col) := as.character(get(value_col))]

      dt <- dt[
        ,
        .SD[1],
        by = .(sample_id, gene_id)
      ]
    } else {
      dt[, (value_col) := suppressWarnings(as.numeric(get(value_col)))]

      dt <- dt[
        ,
        .(
          value_tmp = mean(get(value_col), na.rm = TRUE)
        ),
        by = .(sample_id, gene_id)
      ]

      setnames(dt, "value_tmp", value_col)
    }

    setcolorder(dt, c("sample_id", "gene_id", value_col))

    fwrite(
      dt,
      out_path,
      append = !first_write,
      col.names = first_write
    )

    first_write <- FALSE

    rm(mat_chunk, dt, chunk_map)
    gc(verbose = FALSE)
  }

  invisible(NULL)
}

# -------------------------------------------------------------------------
# Gene table helpers
# -------------------------------------------------------------------------

write_gene_part <- function(gene_maps, out_path) {
  gene_dt <- rbindlist(
    lapply(
      gene_maps,
      function(gene_map) {
        gene_map[, .(id = gene_id, name = gene_name)]
      }
    ),
    fill = TRUE
  )

  gene_dt[, id := strip_ensembl_version(id)]
  gene_dt[, name := clean_na(name)]

  gene_dt <- gene_dt[!is.na(id) & id != ""]

  fwrite(gene_dt, out_path)

  cat("Wrote gene part:", out_path, "rows:", nrow(gene_dt), "\n")
}

finalize_gene_table <- function(gene_part_paths, out_dir) {
  existing_paths <- gene_part_paths[file.exists(gene_part_paths)]

  if (length(existing_paths) == 0) {
    stop("No gene part files found to finalize pre_clinical_gene.csv.")
  }

  gene_dt <- rbindlist(
    lapply(existing_paths, fread),
    fill = TRUE
  )

  gene_dt[, id := strip_ensembl_version(id)]
  gene_dt[, name := clean_na(name)]

  gene_dt <- gene_dt[!is.na(id) & id != ""]

  gene_name_conflicts <- gene_dt[
    !is.na(id) & !is.na(name),
    .(
      names = paste(sort(unique(name)), collapse = "|"),
      n_names = uniqueN(name)
    ),
    by = id
  ][n_names > 1]

  fwrite(
    gene_name_conflicts,
    file.path(out_dir, "pre_clinical_gene_name_conflicts.csv")
  )

  cat("Gene IDs with multiple names:", nrow(gene_name_conflicts), "\n")

  gene_dt <- unique(gene_dt, by = "id")

  fwrite(
    gene_dt,
    file.path(out_dir, "pre_clinical_gene.csv")
  )

  cat("Wrote final genes:", nrow(gene_dt), "\n")
}

# -------------------------------------------------------------------------
# CCLE 2015 treatment response using summarizeSensitivityProfiles
# -------------------------------------------------------------------------

summarized_sensitivity_to_long <- function(
  pset,
  sensitivity_measure,
  cell_lines,
  out_dir
) {
  cat(
    "Summarizing sensitivity measure ",
    sensitivity_measure,
    " with summary.stat='mean'\n",
    sep = ""
  )

  mat <- summarizeSensitivityProfiles(
    object = pset,
    sensitivity.measure = sensitivity_measure,
    cell.lines = cell_lines,
    summary.stat = "mean",
    fill.missing = TRUE,
    verbose = TRUE
  )

  mat <- as.matrix(mat)

  fwrite(
    data.table(
      treatment_id = rownames(mat)
    ),
    file.path(
      out_dir,
      paste0("ccle_2015_treatment_ids_from_", sensitivity_measure, ".csv")
    )
  )

  fwrite(
    data.table(
      sampleid = colnames(mat)
    ),
    file.path(
      out_dir,
      paste0("ccle_2015_cell_lines_from_", sensitivity_measure, ".csv")
    )
  )

  long_dt <- as.data.table(as.table(mat))

  setnames(
    long_dt,
    c("treatment_id", "sampleid", sensitivity_measure)
  )

  long_dt[, treatment_id := clean_na(treatment_id)]
  long_dt[, sampleid := clean_na(sampleid)]
  long_dt[, (sensitivity_measure) := suppressWarnings(as.numeric(get(sensitivity_measure)))]

  long_dt[
    !is.na(treatment_id) &
      !is.na(sampleid)
  ]
}

extract_treatment_response_2015 <- function(pset_2015, canonical_lookup, out_path) {
  canonical_sample_lookup <- unique(
    canonical_lookup[
      source_column %in% c("sampleid", "CCLE.sampleid", "unique.sampleid"),
      .(
        sample_id,
        cell_line_name,
        normalized_sample_id = normalize_cell_line_name(sample_id)
      )
    ],
    by = "sample_id"
  )

  canonical_sample_lookup <- canonical_sample_lookup[
    !is.na(sample_id) &
      !is.na(cell_line_name)
  ]

  if (nrow(canonical_sample_lookup) == 0) {
    stop("No selected CCLE sample IDs found in canonical lookup.")
  }

  ccle_2015_sample_dt <- as.data.table(
    pset_2015@sample,
    keep.rownames = "sample_rowname"
  )

  if (!has_col(ccle_2015_sample_dt, "sampleid")) {
    stop(
      "CCLE_2015@sample must contain sampleid to call ",
      "summarizeSensitivityProfiles by cell.lines."
    )
  }

  available_cell_lines <- unique(
    data.table(
      ccle_2015_sampleid = clean_na(ccle_2015_sample_dt[["sampleid"]]),
      normalized_sample_id = normalize_cell_line_name(ccle_2015_sample_dt[["sampleid"]])
    )
  )

  cell_line_match <- merge(
    canonical_sample_lookup,
    available_cell_lines,
    by = "normalized_sample_id",
    all.x = TRUE
  )

  fwrite(
    cell_line_match,
    sub("\\.csv$", "_selected_cell_lines_available_in_ccle_2015.csv", out_path)
  )

  missing_cell_lines <- cell_line_match[is.na(ccle_2015_sampleid)]

  fwrite(
    missing_cell_lines,
    sub("\\.csv$", "_selected_cell_lines_missing_from_ccle_2015.csv", out_path)
  )

  cell_lines_to_use <- unique(clean_na(cell_line_match$ccle_2015_sampleid))
  cell_lines_to_use <- cell_lines_to_use[!is.na(cell_lines_to_use)]

  if (length(cell_lines_to_use) == 0) {
    stop(
      "None of the selected CCLE 2019 sample IDs matched CCLE 2015 sampleid values. ",
      "See: ",
      sub("\\.csv$", "_selected_cell_lines_missing_from_ccle_2015.csv", out_path)
    )
  }

  aac_long <- summarized_sensitivity_to_long(
    pset = pset_2015,
    sensitivity_measure = "aac_recomputed",
    cell_lines = cell_lines_to_use,
    out_dir = OUT_DIR
  )

  ic50_long <- summarized_sensitivity_to_long(
    pset = pset_2015,
    sensitivity_measure = "ic50_recomputed",
    cell_lines = cell_lines_to_use,
    out_dir = OUT_DIR
  )

  merged_response <- merge(
    aac_long,
    ic50_long,
    by = c("treatment_id", "sampleid"),
    all = TRUE
  )

  merged_response[, normalized_sample_id := normalize_cell_line_name(sampleid)]

  response_sample_lookup <- unique(
    cell_line_match[
      !is.na(ccle_2015_sampleid),
      .(
        normalized_sample_id = normalize_cell_line_name(ccle_2015_sampleid),
        ccle_2015_sampleid,
        canonical_sample_id = sample_id,
        cell_line_name
      )
    ],
    by = "normalized_sample_id"
  )

  treatment_response_mapped <- merge(
    merged_response,
    response_sample_lookup,
    by = "normalized_sample_id",
    all = FALSE
  )

  treatment_response_mapped <- treatment_response_mapped[
    !is.na(aac_recomputed) |
      !is.na(ic50_recomputed)
  ]

  fwrite(
    treatment_response_mapped,
    sub("\\.csv$", "_summarized_mapped_raw_rows.csv", out_path)
  )

  treatment_response_dt <- data.table(
    cell_line_name = clean_na(treatment_response_mapped[["cell_line_name"]]),
    treatment_id = clean_na(treatment_response_mapped[["treatment_id"]]),

    ic50_recomputed = suppressWarnings(
      as.numeric(treatment_response_mapped[["ic50_recomputed"]])
    ),

    # Keeping the existing DB/schema field name.
    # This value comes from summarizeSensitivityProfiles(..., "aac_recomputed").
    acc_recomputed = suppressWarnings(
      as.numeric(treatment_response_mapped[["aac_recomputed"]])
    ),

    mechanism_of_action = NA_character_
  )

  treatment_response_dt <- treatment_response_dt[
    !is.na(cell_line_name) &
      !is.na(treatment_id)
  ]

  setcolorder(
    treatment_response_dt,
    c(
      "cell_line_name",
      "treatment_id",
      "ic50_recomputed",
      "acc_recomputed",
      "mechanism_of_action"
    )
  )

  fwrite(treatment_response_dt, out_path)

  cat("Wrote CCLE 2015 summarized treatment response rows:", nrow(treatment_response_dt), "\n")
}

# -------------------------------------------------------------------------
# Phase 1: Load CCLE 2019 only
# -------------------------------------------------------------------------

cat("\n==============================\n")
cat("Phase 1: Loading CCLE 2019\n")
cat("==============================\n")

ccle_2019 <- read_updated_rds(CCLE_2019_RDS_PATH)

sheet_target_data <- read_sheet_soft_tissue_targets(
  path = SHEET_CELL_LINE_QC_PATH,
  out_dir = OUT_DIR
)

selected_sample_dt <- build_selected_sample_table_2019(
  pset_2019 = ccle_2019,
  target_cell_lines = TARGET_CELL_LINES,
  sheet_soft_tissue_alias_dt = sheet_target_data$sheet_soft_tissue_alias_dt
)

selected_sample_dt <- attach_sheet_cell_line_metadata(
  selected_sample_dt = selected_sample_dt,
  sheet_metadata_alias_dt = sheet_target_data$sheet_metadata_alias_dt
)

canonical_lookup <- build_canonical_lookup(selected_sample_dt)

cat("Requested unique target cell lines:", length(unique(TARGET_CELL_LINES)), "\n")
cat("Matched selected CCLE 2019 sample rows:", nrow(selected_sample_dt), "\n")
cat("Explicit target matches:", selected_sample_dt[target_match == TRUE, .N], "\n")
cat("Soft-tissue matches:", selected_sample_dt[soft_tissue_match == TRUE, .N], "\n")
cat("QC-sheet soft-tissue matches:", selected_sample_dt[sheet_soft_tissue_match == TRUE, .N], "\n")
cat("QC-sheet cell-line metadata matches:", selected_sample_dt[sheet_metadata_match == TRUE, .N], "\n")
cat("Canonical lookup rows:", nrow(canonical_lookup), "\n")

canonical_lookup_path <- file.path(OUT_DIR, "ccle_2019_canonical_lookup.csv")

fwrite(
  canonical_lookup,
  canonical_lookup_path
)

cat("Wrote canonical lookup:", canonical_lookup_path, "\n")

write_cell_line_selection_audits(
  selected_sample_dt = selected_sample_dt,
  sheet_soft_tissue_alias_dt = sheet_target_data$sheet_soft_tissue_alias_dt,
  out_dir = OUT_DIR
)

# -------------------------------------------------------------------------
# pre_clinical_cell_line.csv from QC CSV metadata
# -------------------------------------------------------------------------

cell_line_dt <- data.table(
  cell_line_name = clean_na(selected_sample_dt[["final_cell_line_name"]]),
  accession = clean_na(selected_sample_dt[["sheet_accession"]]),
  category = clean_na(selected_sample_dt[["sheet_category"]]),
  sex = clean_na(selected_sample_dt[["sheet_sex"]]),
  age = parse_age_int(selected_sample_dt[["sheet_age"]])
)

cell_line_dt <- cell_line_dt[!is.na(cell_line_name)]

cell_line_dt[
  ,
  metadata_score :=
    as.integer(!is.na(accession)) +
    as.integer(!is.na(category)) +
    as.integer(!is.na(sex)) +
    as.integer(!is.na(age))
]

setorder(
  cell_line_dt,
  cell_line_name,
  -metadata_score
)

cell_line_dt <- cell_line_dt[
  ,
  .SD[1],
  by = cell_line_name
]

cell_line_dt[, metadata_score := NULL]

setcolorder(
  cell_line_dt,
  c(
    "cell_line_name",
    "accession",
    "category",
    "sex",
    "age"
  )
)

fwrite(
  cell_line_dt,
  file.path(OUT_DIR, "pre_clinical_cell_line.csv")
)

cat("Wrote cell lines:", nrow(cell_line_dt), "\n")

rm(cell_line_dt)
gc(verbose = FALSE)

# -------------------------------------------------------------------------
# pre_clinical_sample.csv from CCLE 2019
# -------------------------------------------------------------------------

sample_out_dt <- data.table(
  id = clean_na(selected_sample_dt[["sampleid"]]),
  cell_line_name = clean_na(selected_sample_dt[["final_cell_line_name"]]),

  site_primary = clean_na(selected_sample_dt[["CCLE.site_Primary"]]),
  site_subtype1 = clean_na(selected_sample_dt[["CCLE.site_Subtype1"]]),
  site_subtype2 = clean_na(selected_sample_dt[["CCLE.site_Subtype2"]]),
  site_subtype3 = clean_na(selected_sample_dt[["CCLE.site_Subtype3"]]),

  histology = clean_na(selected_sample_dt[["CCLE.histology"]]),
  histology_subtype1 = clean_na(selected_sample_dt[["CCLE.histology_Subtype1"]]),
  histology_subtype2 = clean_na(selected_sample_dt[["CCLE.histology_Subtype2"]]),
  histology_subtype3 = clean_na(selected_sample_dt[["CCLE.histology_Subtype3"]]),

  gender = clean_na(selected_sample_dt[["CCLE.gender"]]),
  age = parse_age_int(selected_sample_dt[["CCLE.age"]]),
  race = clean_na(selected_sample_dt[["CCLE.race"]]),
  diseases = clean_na(selected_sample_dt[["cellosaurus.diseases"]]),
  disease_type = clean_na(selected_sample_dt[["CCLE.type"]])
)

sample_out_dt <- sample_out_dt[!is.na(id)]
sample_out_dt <- unique(sample_out_dt, by = "id")

fwrite(
  sample_out_dt,
  file.path(OUT_DIR, "pre_clinical_sample.csv")
)

cat("Wrote samples:", nrow(sample_out_dt), "\n")

rm(sample_out_dt)
gc(verbose = FALSE)

# -------------------------------------------------------------------------
# CCLE 2019 molecular profiles
# -------------------------------------------------------------------------

rnaseq_se <- ccle_2019@molecularProfiles[["rnaseq.gene_tpm"]]
cnv_se <- ccle_2019@molecularProfiles[["cnv.gene_log2"]]
mutation_se <- ccle_2019@molecularProfiles[["mutation.gene_binary"]]

rnaseq_gene_map <- make_gene_mapping(
  rnaseq_se,
  ensembl_candidates = c("gene_id", "EnsemblGeneID", "EnsemblGeneId"),
  name_candidates = c("gene_name", "gene_symbol", "Symbol")
)

cnv_gene_map <- make_gene_mapping(
  cnv_se,
  ensembl_candidates = c("gene_id", "EnsemblGeneID", "EnsemblGeneId"),
  name_candidates = c("gene_name", "gene_symbol", "Symbol")
)

mutation_gene_map <- make_gene_mapping(
  mutation_se,
  ensembl_candidates = c("EnsemblGeneID", "EnsemblGeneId", "gene_id"),
  name_candidates = c("gene_symbol", "Symbol", "gene_name")
)

write_gene_part(
  gene_maps = list(
    rnaseq_gene_map,
    cnv_gene_map,
    mutation_gene_map
  ),
  out_path = file.path(OUT_DIR, "pre_clinical_gene_2019_part.csv")
)

rnaseq_column_map <- build_profile_column_map(rnaseq_se, canonical_lookup)
mutation_column_map <- build_profile_column_map(mutation_se, canonical_lookup)
cnv_column_map <- build_profile_column_map(cnv_se, canonical_lookup)

cat("Matched CCLE 2019 RNA-seq columns:", nrow(rnaseq_column_map), "\n")
cat("Matched CCLE 2019 mutation columns:", nrow(mutation_column_map), "\n")
cat("Matched CCLE 2019 CNV columns:", nrow(cnv_column_map), "\n")

fwrite(
  rnaseq_column_map,
  file.path(OUT_DIR, "ccle_2019_rnaseq_column_map.csv")
)

fwrite(
  mutation_column_map,
  file.path(OUT_DIR, "ccle_2019_mutation_column_map.csv")
)

fwrite(
  cnv_column_map,
  file.path(OUT_DIR, "ccle_2019_cnv_column_map.csv")
)

write_long_assay_with_column_map(
  se = rnaseq_se,
  gene_map = rnaseq_gene_map,
  column_map = rnaseq_column_map,
  out_path = file.path(OUT_DIR, "pre_clinical_rna_seq.csv"),
  value_col = "expression_value",
  assay_name = "exprs",
  value_as_character = FALSE,
  column_chunk_size = 2
)

cat("Wrote CCLE 2019 RNA-seq assay CSV\n")

write_long_assay_with_column_map(
  se = mutation_se,
  gene_map = mutation_gene_map,
  column_map = mutation_column_map,
  out_path = file.path(OUT_DIR, "pre_clinical_mutation.csv"),
  value_col = "value",
  assay_name = "exprs",
  value_as_character = TRUE,
  column_chunk_size = 10
)

cat("Wrote CCLE 2019 mutation assay CSV\n")

write_long_assay_with_column_map(
  se = cnv_se,
  gene_map = cnv_gene_map,
  column_map = cnv_column_map,
  out_path = file.path(OUT_DIR, "pre_clinical_copy_number_variation.csv"),
  value_col = "value",
  assay_name = "exprs",
  value_as_character = FALSE,
  column_chunk_size = 5
)

cat("Wrote CCLE 2019 CNV assay CSV\n")

# -------------------------------------------------------------------------
# Free CCLE 2019 memory before loading CCLE 2015
# -------------------------------------------------------------------------

cat("\nFreeing CCLE 2019 objects from memory\n")

rm(
  ccle_2019,
  selected_sample_dt,
  sheet_target_data,
  canonical_lookup,
  rnaseq_se,
  cnv_se,
  mutation_se,
  rnaseq_gene_map,
  cnv_gene_map,
  mutation_gene_map,
  rnaseq_column_map,
  mutation_column_map,
  cnv_column_map
)

gc(verbose = TRUE)

# -------------------------------------------------------------------------
# Phase 2: Load CCLE 2015 only
# -------------------------------------------------------------------------

cat("\n==============================\n")
cat("Phase 2: Loading CCLE 2015\n")
cat("==============================\n")

ccle_2015 <- read_updated_rds(CCLE_2015_RDS_PATH)

canonical_lookup <- fread(canonical_lookup_path)

# -------------------------------------------------------------------------
# CCLE 2015 treatment response using summarizeSensitivityProfiles
# -------------------------------------------------------------------------

extract_treatment_response_2015(
  pset_2015 = ccle_2015,
  canonical_lookup = canonical_lookup,
  out_path = file.path(OUT_DIR, "pre_clinical_treatment_response.csv")
)

# -------------------------------------------------------------------------
# CCLE 2015 rna profile only: older microarray expression
# -------------------------------------------------------------------------

microarray_se <- ccle_2015@molecularProfiles[["rna"]]

if (is.null(microarray_se)) {
  stop(
    "Could not find CCLE 2015 molecular profile named 'rna'. ",
    "Available profiles are: ",
    paste(names(ccle_2015@molecularProfiles), collapse = ", ")
  )
}

cat("Selected CCLE 2015 microarray RNA profile: rna\n")

microarray_gene_map <- make_gene_mapping(
  microarray_se,
  ensembl_candidates = c("EnsemblGeneId", "EnsemblGeneID", "gene_id"),
  name_candidates = c("Symbol", "gene_symbol", "gene_name")
)

microarray_unmapped <- microarray_gene_map[
  is.na(gene_id) | gene_id == "",
  .(feature_id, gene_name)
]

fwrite(
  microarray_unmapped,
  file.path(OUT_DIR, "ccle_2015_microarray_unmapped_features.csv")
)

cat("CCLE 2015 microarray features without Ensembl ID:", nrow(microarray_unmapped), "\n")

write_gene_part(
  gene_maps = list(
    microarray_gene_map
  ),
  out_path = file.path(OUT_DIR, "pre_clinical_gene_2015_microarray_part.csv")
)

microarray_column_map <- build_2015_rna_column_map(
  rna_se = microarray_se,
  pset_2015 = ccle_2015,
  canonical_lookup = canonical_lookup
)

cat("Matched CCLE 2015 microarray RNA columns:", nrow(microarray_column_map), "\n")

fwrite(
  microarray_column_map,
  file.path(OUT_DIR, "ccle_2015_microarray_column_map.csv")
)

write_long_assay_with_column_map(
  se = microarray_se,
  gene_map = microarray_gene_map,
  column_map = microarray_column_map,
  out_path = file.path(OUT_DIR, "pre_clinical_microarray.csv"),
  value_col = "expression_value",
  assay_name = "exprs",
  value_as_character = FALSE,
  column_chunk_size = 5
)

cat("Wrote CCLE 2015 microarray RNA assay CSV\n")

# -------------------------------------------------------------------------
# Free CCLE 2015 memory
# -------------------------------------------------------------------------

cat("\nFreeing CCLE 2015 objects from memory\n")

rm(
  ccle_2015,
  canonical_lookup,
  microarray_se,
  microarray_gene_map,
  microarray_unmapped,
  microarray_column_map
)

gc(verbose = TRUE)

# -------------------------------------------------------------------------
# Finalize pre_clinical_gene.csv from small gene parts only
# -------------------------------------------------------------------------

finalize_gene_table(
  gene_part_paths = c(
    file.path(OUT_DIR, "pre_clinical_gene_2019_part.csv"),
    file.path(OUT_DIR, "pre_clinical_gene_2015_microarray_part.csv")
  ),
  out_dir = OUT_DIR
)

cat("Finished extracting selected CCLE preclinical CSVs into:", OUT_DIR, "\n")