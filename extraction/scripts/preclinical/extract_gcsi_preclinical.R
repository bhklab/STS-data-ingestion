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

GCSI_2019_RDS_PATH <- "extraction/data/raw/preclinical/GCSI_2019.rds"
SHEET_CELL_LINE_QC_PATH <- "extraction/data/raw/preclinical/All_PSets_sarcoma_cell_line_QC.csv"

OUT_DIR <- "extraction/data/proc/preclinical/GCSI"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# New extraction set:
#   1. Read external QC sheet.
#   2. Keep rows where mod_tissueid == "Soft Tissue".
#   3. Match sheet$cell_line against ONE resolved GCSI sample-slot cell-line column.
SHEET_TARGET_CELL_LINE_COL <- "cell_line"

# Prefer an explicit cell_line_name column if present. Fall back through known
# GCSI/PharmacoSet sample-slot aliases. Only one resolved column is used for
# the new extraction match.
GCSI_CELL_LINE_NAME_CANDIDATES <- c(
  "cell_line_name",
  "Cell.line.primary.name",
  "cellosaurus.cellLineName",
  "sampleid",
  "sample_id",
  "cellLineName",
  "cellline",
  "cellid",
  "id",
  "Name",
  "CCLE.name",
  "sample_rowname"
)

# Original criteria are now audit-only.
# They DO NOT control the final extraction set anymore.
TARGET_CELL_LINES_RAW <- paste0(
  "105KC|SNU-1077|TE 159.T|Aska-SS|GI-1|MES-SA|NCI-H2373|",
  "NCI-H2596|Rh30|SW982|CAL-78|JJ012|TE 125.T"
)

TARGET_CELL_LINES <- unique(trimws(
  unlist(strsplit(TARGET_CELL_LINES_RAW, "\\|"))
))

SAMPLE_ID_PREFIX <- "gcsi_"

RNA_PROFILE_NAME <- "Kallisto_0.46.1.rnaseq"
CNV_PROFILE_NAME <- "cnv"
MUTATION_PROFILE_NAME <- "mutation"

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

make_prefixed_sample_id <- function(cell_line_name) {
  cell_line_name <- clean_na(cell_line_name)
  out <- rep(NA_character_, length(cell_line_name))
  keep_idx <- !is.na(cell_line_name)

  out[keep_idx] <- ifelse(
    startsWith(cell_line_name[keep_idx], SAMPLE_ID_PREFIX),
    cell_line_name[keep_idx],
    paste0(SAMPLE_ID_PREFIX, cell_line_name[keep_idx])
  )

  out
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

safe_col <- function(dt, candidates, n = nrow(dt), default = NA_character_) {
  col <- first_existing_col(dt, candidates)

  if (is.na(col)) {
    return(rep(default, n))
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

  # Fallback: if rownames are Ensembl IDs, use rownames as gene_id.
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

get_assay_matrix <- function(se, assay_name = "exprs") {
  if (!(assay_name %in% assayNames(se))) {
    assay_name <- assayNames(se)[1]
  }

  assay(se, assay_name)
}

require_profile <- function(pset, profile_name) {
  profile <- pset@molecularProfiles[[profile_name]]

  if (is.null(profile)) {
    stop(
      "Could not find molecular profile named '",
      profile_name,
      "'. Available profiles are: ",
      paste(names(pset@molecularProfiles), collapse = ", ")
    )
  }

  profile
}

resolve_gcsi_cell_line_col <- function(sample_dt) {
  col <- first_existing_col(
    sample_dt,
    GCSI_CELL_LINE_NAME_CANDIDATES
  )

  if (is.na(col)) {
    stop(
      "Could not find a GCSI sample-slot cell-line-name column. ",
      "Expected one of: ",
      paste(GCSI_CELL_LINE_NAME_CANDIDATES, collapse = ", "),
      ". Available columns are: ",
      paste(colnames(sample_dt), collapse = ", ")
    )
  }

  col
}

# -------------------------------------------------------------------------
# QC-sheet helpers
# -------------------------------------------------------------------------

read_qc_sheet_soft_tissue_targets <- function(path, out_dir) {
  if (!file.exists(path)) {
    stop(
      "QC sheet file not found at: ",
      path,
      ". Copy All_PSets_sarcoma_cell_line_QC.csv to that path."
    )
  }

  sheet_dt <- fread(path)
  sheet_dt[, sheet_row_id := .I]

  fwrite(
    data.table(column_name = colnames(sheet_dt)),
    file.path(out_dir, "gcsi_sheet_column_names.csv")
  )

  required_cols <- c(
    SHEET_TARGET_CELL_LINE_COL,
    "tissueid",
    "mod_tissueid"
  )

  missing_required <- setdiff(required_cols, colnames(sheet_dt))

  if (length(missing_required) > 0) {
    stop(
      "QC sheet is missing required columns: ",
      paste(missing_required, collapse = ", "),
      ". See: ",
      file.path(out_dir, "gcsi_sheet_column_names.csv")
    )
  }

  sheet_soft_tissue_dt <- sheet_dt[
    normalize_category_value(mod_tissueid) == "soft_tissue"
  ]

  if (nrow(sheet_soft_tissue_dt) == 0) {
    stop(
      "No rows found in QC sheet with mod_tissueid == 'Soft Tissue'."
    )
  }

  sheet_soft_tissue_dt[
    ,
    sheet_cell_line_name := clean_na(get(SHEET_TARGET_CELL_LINE_COL))
  ]

  sheet_soft_tissue_dt[
    ,
    normalized_sheet_cell_line_name := normalize_cell_line_name(sheet_cell_line_name)
  ]

  sheet_soft_tissue_dt <- sheet_soft_tissue_dt[
    !is.na(sheet_cell_line_name) &
      !is.na(normalized_sheet_cell_line_name)
  ]

  if (nrow(sheet_soft_tissue_dt) == 0) {
    stop(
      "QC sheet has mod_tissueid == 'Soft Tissue' rows, but none have a usable ",
      SHEET_TARGET_CELL_LINE_COL,
      " value."
    )
  }

  sheet_soft_tissue_dt[, tissueid := clean_na(tissueid)]
  sheet_soft_tissue_dt[, mod_tissueid := clean_na(mod_tissueid)]

  sheet_soft_tissue_dt[
    ,
    accession := clean_na(coalesce_dt_cols(
      .SD,
      c(
        "Cellosaurus.Accession.id",
        "cellosaurus.cvcl_id",
        "cellosaurus.cellosaurus.cvcl_id",
        "cellosaurus.accession",
        "accession"
      )
    ))
  ]

  sheet_soft_tissue_dt[
    ,
    category := clean_na(coalesce_dt_cols(
      .SD,
      c(
        "cellosaurus.category",
        "cellosaurus.cellosaurus.category",
        "CellLine.Type",
        "category"
      )
    ))
  ]

  sheet_soft_tissue_dt[
    ,
    sex := clean_na(coalesce_dt_cols(
      .SD,
      c(
        "cellosaurus.sex",
        "cellosaurus.cellosaurus.sex",
        "cellosaurus.sexOfCell",
        "Gender",
        "gender",
        "Sex",
        "sex"
      )
    ))
  ]

  sheet_soft_tissue_dt[
    ,
    age_raw := clean_na(coalesce_dt_cols(
      .SD,
      c(
        "cellosaurus.age",
        "cellosaurus.cellosaurus.age",
        "cellosaurus.ageAtSampling",
        "Age",
        "age"
      )
    ))
  ]

  sheet_soft_tissue_dt[
    ,
    sheet_dataset_value := clean_na(safe_col(.SD, c("dataset"), n = .N))
  ]

  sheet_soft_tissue_dt[
    ,
    sheet_object_type_value := clean_na(safe_col(.SD, c("object_type"), n = .N))
  ]

  sheet_soft_tissue_dt[
    ,
    sheet_sampleid_value := clean_na(safe_col(.SD, c("sampleid"), n = .N))
  ]

  sheet_soft_tissue_dt[
    ,
    metadata_priority := fifelse(
      grepl("GCSI", toupper(sheet_dataset_value)),
      1L,
      fifelse(grepl("GCSI", toupper(sheet_object_type_value)), 2L, 99L)
    )
  ]

  sheet_soft_tissue_dt[
    ,
    metadata_score :=
      as.integer(!is.na(sheet_cell_line_name)) +
      as.integer(!is.na(tissueid)) +
      as.integer(!is.na(mod_tissueid)) +
      as.integer(!is.na(accession)) +
      as.integer(!is.na(category)) +
      as.integer(!is.na(sex)) +
      as.integer(!is.na(age_raw))
  ]

  setorder(
    sheet_soft_tissue_dt,
    normalized_sheet_cell_line_name,
    metadata_priority,
    -metadata_score
  )

  sheet_soft_tissue_metadata_dt <- sheet_soft_tissue_dt[
    ,
    .SD[1],
    by = normalized_sheet_cell_line_name
  ]

  fwrite(
    sheet_soft_tissue_dt,
    file.path(out_dir, "gcsi_sheet_soft_tissue_rows.csv")
  )

  fwrite(
    sheet_soft_tissue_metadata_dt,
    file.path(out_dir, "gcsi_sheet_soft_tissue_cell_line_metadata.csv")
  )

  fwrite(
    unique(
      sheet_soft_tissue_metadata_dt[
        ,
        .(
          sheet_cell_line_name,
          normalized_sheet_cell_line_name,
          dataset = sheet_dataset_value,
          object_type = sheet_object_type_value,
          sampleid = sheet_sampleid_value,
          tissueid,
          mod_tissueid,
          accession,
          category,
          sex,
          age_raw
        )
      ]
    ),
    file.path(out_dir, "gcsi_sheet_soft_tissue_cell_line_targets.csv")
  )

  cat("QC-sheet soft tissue rows:", nrow(sheet_soft_tissue_dt), "\n")
  cat("QC-sheet unique soft tissue cell lines:", nrow(sheet_soft_tissue_metadata_dt), "\n")

  list(
    sheet_soft_tissue_dt = sheet_soft_tissue_dt,
    sheet_soft_tissue_metadata_dt = sheet_soft_tissue_metadata_dt
  )
}

# -------------------------------------------------------------------------
# GCSI 2019 sample selection and original-criteria audit
# -------------------------------------------------------------------------

build_gcsi_selection <- function(
  pset_gcsi,
  target_cell_lines,
  sheet_soft_tissue_metadata_dt,
  out_dir
) {
  sample_dt <- as.data.table(pset_gcsi@sample, keep.rownames = "sample_rowname")

  fwrite(
    data.table(column_name = colnames(sample_dt)),
    file.path(out_dir, "gcsi_sample_slot_columns.csv")
  )

  gcsi_cell_line_col <- resolve_gcsi_cell_line_col(sample_dt)

  fwrite(
    data.table(
      gcsi_cell_line_name_column_used = gcsi_cell_line_col
    ),
    file.path(out_dir, "gcsi_cell_line_name_column_used.csv")
  )

  sample_dt[
    ,
    gcsi_cell_line_name := clean_na(get(gcsi_cell_line_col))
  ]

  sample_dt[
    ,
    normalized_gcsi_cell_line_name := normalize_cell_line_name(gcsi_cell_line_name)
  ]

  # New extraction criterion only:
  # QC sheet mod_tissueid == "Soft Tissue", matched by sheet$cell_line to the
  # one resolved GCSI sample-slot cell-line-name column.
  sample_dt[
    ,
    new_sheet_soft_tissue_match :=
      normalized_gcsi_cell_line_name %in%
        sheet_soft_tissue_metadata_dt$normalized_sheet_cell_line_name
  ]

  selected_sample_dt <- sample_dt[
    new_sheet_soft_tissue_match == TRUE
  ]

  if (nrow(selected_sample_dt) == 0) {
    stop(
      "No GCSI samples matched QC-sheet mod_tissueid == 'Soft Tissue' ",
      "cell lines using GCSI sample-slot column: ",
      gcsi_cell_line_col
    )
  }

  selected_sample_dt <- merge(
    selected_sample_dt,
    sheet_soft_tissue_metadata_dt[
      ,
      .(
        normalized_sheet_cell_line_name,
        sheet_cell_line_name,
        sheet_dataset = sheet_dataset_value,
        sheet_object_type = sheet_object_type_value,
        sheet_sampleid = sheet_sampleid_value,
        sheet_tissueid = tissueid,
        sheet_mod_tissueid = mod_tissueid,
        sheet_accession = accession,
        sheet_category = category,
        sheet_sex = sex,
        sheet_age = age_raw
      )
    ],
    by.x = "normalized_gcsi_cell_line_name",
    by.y = "normalized_sheet_cell_line_name",
    all.x = TRUE
  )

  selected_sample_dt[
    ,
    canonical_cell_line_name := clean_na(sheet_cell_line_name)
  ]

  selected_sample_dt[
    is.na(canonical_cell_line_name),
    canonical_cell_line_name := gcsi_cell_line_name
  ]

  selected_sample_dt[
    ,
    final_cell_line_name := canonical_cell_line_name
  ]

  selected_sample_dt[
    ,
    source_sampleid := canonical_cell_line_name
  ]

  selected_sample_dt[
    ,
    canonical_sample_id := make_prefixed_sample_id(source_sampleid)
  ]

  # -----------------------------------------------------------------------
  # Original criteria audit only.
  # These rows do not change selected_sample_dt.
  # -----------------------------------------------------------------------

  original_target_norm <- normalize_cell_line_name(target_cell_lines)

  original_candidate_cols <- c(
    "cellosaurus.cellLineName",
    "Cell.line.primary.name",
    "cell_line_name",
    "cellLineName",
    "cellline",
    "cellid",
    "sampleid",
    "sample_id",
    "id",
    "Name",
    "CCLE.name",
    "rownames",
    "sample_rowname"
  )

  sample_dt[, original_target_match := FALSE]

  for (col in original_candidate_cols) {
    if (has_col(sample_dt, col)) {
      sample_dt[
        normalize_cell_line_name(get(col)) %in% original_target_norm,
        original_target_match := TRUE
      ]
    }
  }

  sample_dt[, original_soft_tissue_match := FALSE]

  if (has_col(sample_dt, "tissueid")) {
    sample_dt[
      normalize_category_value(tissueid) == "soft_tissue",
      original_soft_tissue_match := TRUE
    ]
  }

  original_criteria_dt <- sample_dt[
    original_target_match == TRUE |
      original_soft_tissue_match == TRUE
  ]

  original_criteria_dt[
    ,
    original_cell_line_name := gcsi_cell_line_name
  ]

  original_criteria_dt[
    ,
    normalized_original_cell_line_name := normalize_cell_line_name(original_cell_line_name)
  ]

  final_norm <- unique(
    normalize_cell_line_name(selected_sample_dt$final_cell_line_name)
  )

  original_not_in_final_dt <- original_criteria_dt[
    !(normalized_original_cell_line_name %in% final_norm)
  ]

  original_audit_dt <- data.table(
    source_sampleid = clean_na(coalesce_dt_cols(
      original_criteria_dt,
      c("sampleid", "sample_id", "cellid", "id", "sample_rowname")
    )),
    sample_id = make_prefixed_sample_id(clean_na(original_criteria_dt$original_cell_line_name)),
    cell_line_name = clean_na(original_criteria_dt$original_cell_line_name),
    gcsi_cell_line_name_column_used = gcsi_cell_line_col,
    tissueid = clean_na(safe_col(original_criteria_dt, c("tissueid"))),
    original_target_match = original_criteria_dt$original_target_match,
    original_soft_tissue_match = original_criteria_dt$original_soft_tissue_match
  )

  original_not_in_final_audit_dt <- data.table(
    source_sampleid = clean_na(coalesce_dt_cols(
      original_not_in_final_dt,
      c("sampleid", "sample_id", "cellid", "id", "sample_rowname")
    )),
    sample_id = make_prefixed_sample_id(clean_na(original_not_in_final_dt$original_cell_line_name)),
    cell_line_name = clean_na(original_not_in_final_dt$original_cell_line_name),
    gcsi_cell_line_name_column_used = gcsi_cell_line_col,
    tissueid = clean_na(safe_col(original_not_in_final_dt, c("tissueid"))),
    original_target_match = original_not_in_final_dt$original_target_match,
    original_soft_tissue_match = original_not_in_final_dt$original_soft_tissue_match
  )

  fwrite(
    original_audit_dt,
    file.path(out_dir, "gcsi_original_criteria_cell_line_list_audit_only.csv")
  )

  fwrite(
    original_not_in_final_audit_dt,
    file.path(out_dir, "gcsi_original_criteria_not_in_final_new_criteria.csv")
  )

  fwrite(
    selected_sample_dt[
      ,
      .(
        sample_id = canonical_sample_id,
        source_sampleid,
        cell_line_name = final_cell_line_name,
        gcsi_cell_line_name,
        gcsi_cell_line_name_column_used = gcsi_cell_line_col,
        sheet_cell_line_name,
        sheet_dataset,
        sheet_object_type,
        sheet_sampleid,
        sheet_tissueid,
        sheet_mod_tissueid,
        sheet_accession,
        sheet_category,
        sheet_sex,
        sheet_age
      )
    ],
    file.path(out_dir, "gcsi_final_extracted_cell_line_list_new_criteria.csv")
  )

  fwrite(
    selected_sample_dt,
    file.path(out_dir, "gcsi_2019_selected_sample_slot_rows.csv")
  )

  cat("GCSI cell-line-name column used:", gcsi_cell_line_col, "\n")
  cat("Final new-criteria selected GCSI sample rows:", nrow(selected_sample_dt), "\n")
  cat("Original-criteria audit rows:", nrow(original_criteria_dt), "\n")
  cat("Original-criteria rows not in final new criteria:", nrow(original_not_in_final_dt), "\n")

  list(
    selected_sample_dt = selected_sample_dt,
    original_criteria_dt = original_criteria_dt,
    original_not_in_final_dt = original_not_in_final_dt,
    gcsi_cell_line_col = gcsi_cell_line_col
  )
}

build_canonical_lookup_gcsi <- function(selected_sample_dt) {
  alt_cols <- c(
    "canonical_cell_line_name",
    "final_cell_line_name",
    "canonical_sample_id",
    "source_sampleid",
    "gcsi_cell_line_name",
    "cellosaurus.cellLineName",
    "Cell.line.primary.name",
    "cell_line_name",
    "cellLineName",
    "cellline",
    "cellid",
    "sampleid",
    "sample_id",
    "id",
    "Name",
    "CCLE.name",
    "rownames",
    "sample_rowname"
  )

  pieces <- list()

  for (col in alt_cols) {
    if (has_col(selected_sample_dt, col)) {
      pieces[[col]] <- data.table(
        lookup_name = clean_na(selected_sample_dt[[col]]),
        normalized_lookup_name = normalize_cell_line_name(selected_sample_dt[[col]]),
        sample_id = selected_sample_dt$canonical_sample_id,
        source_sampleid = selected_sample_dt$source_sampleid,
        cell_line_name = selected_sample_dt$canonical_cell_line_name,
        source_column = col
      )
    }
  }

  lookup_dt <- rbindlist(pieces, fill = TRUE)

  lookup_dt <- lookup_dt[
    !is.na(normalized_lookup_name) &
      !is.na(sample_id) &
      !is.na(source_sampleid) &
      !is.na(cell_line_name)
  ]

  unique(
    lookup_dt,
    by = c("normalized_lookup_name", "sample_id", "source_column")
  )
}

# -------------------------------------------------------------------------
# Molecular profile column mapping
# -------------------------------------------------------------------------

build_profile_column_map <- function(se, canonical_lookup, profile_label) {
  cd <- as.data.table(
    as.data.frame(colData(se)),
    keep.rownames = "coldata_rownames"
  )

  cd[, colname := colnames(se)]

  fwrite(
    data.table(column_name = colnames(cd)),
    file.path(OUT_DIR, paste0("gcsi_", profile_label, "_coldata_columns.csv"))
  )

  candidate_cols <- c(
    "cellosaurus.cellLineName",
    "Cell.line.primary.name",
    "cell_line_name",
    "cellLineName",
    "cellline",
    "Cell_Line",
    "sampleid",
    "sample_id",
    "ccle_sample_id",
    "dataset_sample_id",
    "cellid",
    "id",
    "Name",
    "CCLE.name",
    "samplename",
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
    warning(
      "No candidate mapping columns found in colData for profile ",
      profile_label,
      ". See gcsi_",
      profile_label,
      "_coldata_columns.csv."
    )

    return(data.table(
      colname = character(),
      sample_id = character(),
      source_sampleid = character(),
      cell_line_name = character(),
      matched_column = character()
    ))
  }

  candidate_map <- rbindlist(pieces, fill = TRUE)
  candidate_map <- candidate_map[!is.na(normalized_lookup_name)]

  canonical_lookup_unique <- unique(
    canonical_lookup[, .(normalized_lookup_name, sample_id, source_sampleid, cell_line_name)],
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

  mapped[, .(colname, sample_id, source_sampleid, cell_line_name, matched_column)]
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
  duplicate_pair_path <- sub("\\.csv$", "_duplicate_sample_gene_pairs.csv", out_path)

  if (file.exists(duplicate_pair_path)) {
    file.remove(duplicate_pair_path)
  }

  duplicate_pair_first_write <- TRUE

  if (nrow(column_map) == 0) {
    empty_dt <- data.table(
      sample_id = character(),
      gene_id = character()
    )

    empty_dt[[value_col]] <- if (value_as_character) character() else numeric()

    fwrite(empty_dt, out_path)

    fwrite(
      data.table(
        sample_id = character(),
        gene_id = character(),
        N = integer()
      ),
      duplicate_pair_path
    )

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

    fwrite(
      data.table(
        sample_id = character(),
        gene_id = character(),
        N = integer()
      ),
      duplicate_pair_path
    )

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

    duplicate_pairs_chunk <- dt[
      ,
      .N,
      by = .(sample_id, gene_id)
    ][N > 1]

    fwrite(
      duplicate_pairs_chunk,
      duplicate_pair_path,
      append = !duplicate_pair_first_write,
      col.names = duplicate_pair_first_write
    )

    duplicate_pair_first_write <- FALSE

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

    rm(mat_chunk, dt, chunk_map, duplicate_pairs_chunk)
    gc(verbose = FALSE)
  }

  if (!file.exists(duplicate_pair_path)) {
    fwrite(
      data.table(
        sample_id = character(),
        gene_id = character(),
        N = integer()
      ),
      duplicate_pair_path
    )
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
# GCSI treatment response using summarizeSensitivityProfiles
# -------------------------------------------------------------------------

get_available_cell_lines_for_sensitivity <- function(pset_gcsi) {
  out <- tryCatch(
    cellNames(pset_gcsi),
    error = function(e) {
      character()
    }
  )

  out <- clean_na(out)
  out <- out[!is.na(out)]

  if (length(out) > 0) {
    return(unique(out))
  }

  sample_dt <- as.data.table(
    pset_gcsi@sample,
    keep.rownames = "sample_rowname"
  )

  candidate_col <- first_existing_col(
    sample_dt,
    GCSI_CELL_LINE_NAME_CANDIDATES
  )

  if (is.na(candidate_col)) {
    stop(
      "Could not determine available GCSI cell lines. ",
      "cellNames(pset_gcsi) failed and no usable sample column was found."
    )
  }

  unique(clean_na(sample_dt[[candidate_col]]))
}

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

  # Expected PharmacoGx shape is treatment/drug IDs as rows and cell lines as
  # columns. If an object returns the reverse orientation, transpose it.
  input_norm <- normalize_cell_line_name(cell_lines)
  row_cell_count <- sum(normalize_cell_line_name(rownames(mat)) %in% input_norm, na.rm = TRUE)
  col_cell_count <- sum(normalize_cell_line_name(colnames(mat)) %in% input_norm, na.rm = TRUE)

  if (row_cell_count > col_cell_count) {
    mat <- t(mat)
  }

  fwrite(
    data.table(
      treatment_id = rownames(mat)
    ),
    file.path(
      out_dir,
      paste0("gcsi_treatment_ids_from_", sensitivity_measure, ".csv")
    )
  )

  fwrite(
    data.table(
      cell_line_name_for_summary = colnames(mat)
    ),
    file.path(
      out_dir,
      paste0("gcsi_cell_lines_from_", sensitivity_measure, ".csv")
    )
  )

  long_dt <- as.data.table(as.table(mat))

  setnames(
    long_dt,
    c("treatment_id", "cell_line_name_for_summary", sensitivity_measure)
  )

  long_dt[, treatment_id := clean_na(treatment_id)]
  long_dt[, cell_line_name_for_summary := clean_na(cell_line_name_for_summary)]
  long_dt[, (sensitivity_measure) := suppressWarnings(as.numeric(get(sensitivity_measure)))]

  long_dt[
    !is.na(treatment_id) &
      !is.na(cell_line_name_for_summary)
  ]
}

extract_treatment_response_gcsi <- function(pset_gcsi, canonical_lookup, out_path) {
  selected_cell_line_lookup <- unique(
    canonical_lookup[
      source_column %in% c("canonical_cell_line_name", "final_cell_line_name"),
      .(
        sample_id,
        source_sampleid,
        cell_line_name,
        normalized_cell_line_name = normalize_cell_line_name(cell_line_name)
      )
    ],
    by = c("sample_id", "cell_line_name")
  )

  selected_cell_line_lookup <- selected_cell_line_lookup[
    !is.na(sample_id) &
      !is.na(cell_line_name) &
      !is.na(normalized_cell_line_name)
  ]

  if (nrow(selected_cell_line_lookup) == 0) {
    stop(
      "No selected cell_line_name values found in canonical lookup. ",
      "Treatment response requires selected final/canonical cell-line names."
    )
  }

  available_cell_lines <- get_available_cell_lines_for_sensitivity(pset_gcsi)

  available_lookup <- data.table(
    gcsi_cell_line_name_for_summary = available_cell_lines,
    normalized_cell_line_name = normalize_cell_line_name(available_cell_lines)
  )

  available_lookup <- unique(
    available_lookup[
      !is.na(gcsi_cell_line_name_for_summary) &
        !is.na(normalized_cell_line_name)
    ],
    by = "normalized_cell_line_name"
  )

  cell_line_match <- merge(
    selected_cell_line_lookup,
    available_lookup,
    by = "normalized_cell_line_name",
    all.x = TRUE
  )

  fwrite(
    cell_line_match,
    sub("\\.csv$", "_selected_cell_lines_available_in_gcsi.csv", out_path)
  )

  missing_cell_lines <- cell_line_match[is.na(gcsi_cell_line_name_for_summary)]

  fwrite(
    missing_cell_lines,
    sub("\\.csv$", "_selected_cell_lines_missing_from_gcsi.csv", out_path)
  )

  cell_lines_to_use <- unique(clean_na(cell_line_match$gcsi_cell_line_name_for_summary))
  cell_lines_to_use <- cell_lines_to_use[!is.na(cell_lines_to_use)]

  if (length(cell_lines_to_use) == 0) {
    stop(
      "None of the selected GCSI cell_line_name values matched sensitivity cell lines. ",
      "See: ",
      sub("\\.csv$", "_selected_cell_lines_missing_from_gcsi.csv", out_path)
    )
  }

  aac_long <- summarized_sensitivity_to_long(
    pset = pset_gcsi,
    sensitivity_measure = "aac_recomputed",
    cell_lines = cell_lines_to_use,
    out_dir = OUT_DIR
  )

  ic50_long <- summarized_sensitivity_to_long(
    pset = pset_gcsi,
    sensitivity_measure = "ic50_recomputed",
    cell_lines = cell_lines_to_use,
    out_dir = OUT_DIR
  )

  merged_response <- merge(
    aac_long,
    ic50_long,
    by = c("treatment_id", "cell_line_name_for_summary"),
    all = TRUE
  )

  merged_response[
    ,
    normalized_cell_line_name := normalize_cell_line_name(cell_line_name_for_summary)
  ]

  response_cell_line_lookup <- unique(
    cell_line_match[
      !is.na(gcsi_cell_line_name_for_summary),
      .(
        normalized_cell_line_name,
        gcsi_cell_line_name_for_summary,
        sample_id,
        source_sampleid,
        cell_line_name
      )
    ],
    by = "normalized_cell_line_name"
  )

  treatment_response_mapped <- merge(
    merged_response,
    response_cell_line_lookup,
    by = "normalized_cell_line_name",
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

    # Keeping the shared DB/schema field name.
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

  # summarizeSensitivityProfiles(summary.stat = "mean") should already collapse
  # raw replicate experiments. This check only catches accidental duplicate
  # mappings caused by duplicated selected cell-line names.
  duplicate_after_summary <- treatment_response_dt[
    ,
    .N,
    by = .(cell_line_name, treatment_id)
  ][N > 1]

  fwrite(
    duplicate_after_summary,
    sub("\\.csv$", "_duplicates_after_summarizeSensitivityProfiles.csv", out_path)
  )

  if (nrow(duplicate_after_summary) > 0) {
    warning(
      "Duplicates remained after summarizeSensitivityProfiles because multiple selected samples map to the same cell_line_name. ",
      "Keeping first row per cell_line_name/treatment_id. See: ",
      sub("\\.csv$", "_duplicates_after_summarizeSensitivityProfiles.csv", out_path)
    )

    treatment_response_dt <- unique(
      treatment_response_dt,
      by = c("cell_line_name", "treatment_id")
    )
  }

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

  cat("Wrote GCSI summarized treatment response rows:", nrow(treatment_response_dt), "\n")
}

# -------------------------------------------------------------------------
# Load GCSI 2019
# -------------------------------------------------------------------------

cat("\n==============================\n")
cat("Loading GCSI 2019\n")
cat("==============================\n")

gcsi_2019 <- read_updated_rds(GCSI_2019_RDS_PATH)

# -------------------------------------------------------------------------
# Select samples/cell lines using external QC sheet new criteria only
# -------------------------------------------------------------------------

sheet_target_data <- read_qc_sheet_soft_tissue_targets(
  path = SHEET_CELL_LINE_QC_PATH,
  out_dir = OUT_DIR
)

selection_data <- build_gcsi_selection(
  pset_gcsi = gcsi_2019,
  target_cell_lines = TARGET_CELL_LINES,
  sheet_soft_tissue_metadata_dt = sheet_target_data$sheet_soft_tissue_metadata_dt,
  out_dir = OUT_DIR
)

selected_sample_dt <- selection_data$selected_sample_dt

canonical_lookup <- build_canonical_lookup_gcsi(selected_sample_dt)

cat("Final new-criteria selected sample rows:", nrow(selected_sample_dt), "\n")
cat("Canonical lookup rows:", nrow(canonical_lookup), "\n")

fwrite(
  canonical_lookup,
  file.path(OUT_DIR, "gcsi_2019_canonical_lookup.csv")
)

# -------------------------------------------------------------------------
# pre_clinical_cell_line.csv from QC CSV metadata
# -------------------------------------------------------------------------

cell_line_dt <- data.table(
  cell_line_name = clean_na(selected_sample_dt[["canonical_cell_line_name"]]),
  tissueid = clean_na(selected_sample_dt[["sheet_tissueid"]]),
  mod_tissueid = clean_na(selected_sample_dt[["sheet_mod_tissueid"]]),
  accession = clean_na(selected_sample_dt[["sheet_accession"]]),
  category = clean_na(selected_sample_dt[["sheet_category"]]),
  sex = clean_na(selected_sample_dt[["sheet_sex"]]),
  age = parse_age_int(selected_sample_dt[["sheet_age"]])
)

cell_line_dt <- cell_line_dt[!is.na(cell_line_name)]

cell_line_dt[
  ,
  metadata_score :=
    as.integer(!is.na(tissueid)) +
    as.integer(!is.na(mod_tissueid)) +
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
    "tissueid",
    "mod_tissueid",
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

# -------------------------------------------------------------------------
# pre_clinical_sample.csv
# -------------------------------------------------------------------------

sample_out_dt <- data.table(
  id = clean_na(selected_sample_dt[["canonical_sample_id"]]),
  cell_line_name = clean_na(selected_sample_dt[["canonical_cell_line_name"]]),

  site_primary = clean_na(safe_col(
    selected_sample_dt,
    c("site_primary", "Site.Primary", "primary_site", "tissueid")
  )),

  site_subtype1 = clean_na(safe_col(
    selected_sample_dt,
    c("site_subtype1", "Site.Subtype1", "Site_Subtype1")
  )),

  site_subtype2 = clean_na(safe_col(
    selected_sample_dt,
    c("site_subtype2", "Site.Subtype2", "Site_Subtype2")
  )),

  site_subtype3 = clean_na(safe_col(
    selected_sample_dt,
    c("site_subtype3", "Site.Subtype3", "Site_Subtype3")
  )),

  histology = clean_na(safe_col(
    selected_sample_dt,
    c("histology", "Histology")
  )),

  histology_subtype1 = clean_na(safe_col(
    selected_sample_dt,
    c("histology_subtype1", "Hist.Subtype1", "Histology_Subtype1")
  )),

  histology_subtype2 = clean_na(safe_col(
    selected_sample_dt,
    c("histology_subtype2", "Hist.Subtype2", "Histology_Subtype2")
  )),

  histology_subtype3 = clean_na(safe_col(
    selected_sample_dt,
    c("histology_subtype3", "Hist.Subtype3", "Histology_Subtype3")
  )),

  gender = clean_na(safe_col(
    selected_sample_dt,
    c("gender", "Gender", "sex", "Sex")
  )),

  age = parse_age_int(safe_col(
    selected_sample_dt,
    c("age", "Age")
  )),

  race = clean_na(safe_col(
    selected_sample_dt,
    c("race", "Race")
  )),

  diseases = clean_na(safe_col(
    selected_sample_dt,
    c("diseases", "cellosaurus.diseases", "Disease", "disease")
  )),

  disease_type = clean_na(safe_col(
    selected_sample_dt,
    c("disease_type", "type", "tissueid")
  ))
)

sample_out_dt <- sample_out_dt[!is.na(id)]
sample_out_dt <- unique(sample_out_dt, by = "id")

fwrite(
  sample_out_dt,
  file.path(OUT_DIR, "pre_clinical_sample.csv")
)

cat("Wrote samples:", nrow(sample_out_dt), "\n")

# -------------------------------------------------------------------------
# Treatment response using summarizeSensitivityProfiles
# -------------------------------------------------------------------------

extract_treatment_response_gcsi(
  pset_gcsi = gcsi_2019,
  canonical_lookup = canonical_lookup,
  out_path = file.path(OUT_DIR, "pre_clinical_treatment_response.csv")
)

# -------------------------------------------------------------------------
# Existing GCSI molecular profiles only. No new molecular profiles are added.
# -------------------------------------------------------------------------

rnaseq_se <- require_profile(gcsi_2019, RNA_PROFILE_NAME)
cnv_se <- require_profile(gcsi_2019, CNV_PROFILE_NAME)
mutation_se <- require_profile(gcsi_2019, MUTATION_PROFILE_NAME)

# -------------------------------------------------------------------------
# Gene mappings
# -------------------------------------------------------------------------

rnaseq_gene_map <- make_gene_mapping(
  rnaseq_se,
  ensembl_candidates = c(
    "gene_id",
    "EnsemblGeneID",
    "EnsemblGeneId",
    "ensembl_gene_id"
  ),
  name_candidates = c(
    "gene_name",
    "gene_symbol",
    "Symbol",
    "symbol",
    "GeneSymbol"
  )
)

cnv_gene_map <- make_gene_mapping(
  cnv_se,
  ensembl_candidates = c(
    "gene_id",
    "EnsemblGeneID",
    "EnsemblGeneId",
    "ensembl_gene_id"
  ),
  name_candidates = c(
    "gene_name",
    "gene_symbol",
    "Symbol",
    "symbol",
    "GeneSymbol"
  )
)

mutation_gene_map <- make_gene_mapping(
  mutation_se,
  ensembl_candidates = c(
    "EnsemblGeneID",
    "EnsemblGeneId",
    "gene_id",
    "ensembl_gene_id"
  ),
  name_candidates = c(
    "gene_symbol",
    "Symbol",
    "symbol",
    "gene_name",
    "GeneSymbol"
  )
)

write_gene_part(
  gene_maps = list(
    rnaseq_gene_map,
    cnv_gene_map,
    mutation_gene_map
  ),
  out_path = file.path(OUT_DIR, "pre_clinical_gene_gcsi_part.csv")
)

finalize_gene_table(
  gene_part_paths = c(
    file.path(OUT_DIR, "pre_clinical_gene_gcsi_part.csv")
  ),
  out_dir = OUT_DIR
)

# -------------------------------------------------------------------------
# Molecular column maps
# -------------------------------------------------------------------------

rnaseq_column_map <- build_profile_column_map(
  se = rnaseq_se,
  canonical_lookup = canonical_lookup,
  profile_label = "rnaseq"
)

cnv_column_map <- build_profile_column_map(
  se = cnv_se,
  canonical_lookup = canonical_lookup,
  profile_label = "cnv"
)

mutation_column_map <- build_profile_column_map(
  se = mutation_se,
  canonical_lookup = canonical_lookup,
  profile_label = "mutation"
)

cat("Matched GCSI RNA-seq columns:", nrow(rnaseq_column_map), "\n")
cat("Matched GCSI CNV columns:", nrow(cnv_column_map), "\n")
cat("Matched GCSI mutation columns:", nrow(mutation_column_map), "\n")

fwrite(
  rnaseq_column_map,
  file.path(OUT_DIR, "gcsi_2019_rnaseq_column_map.csv")
)

fwrite(
  cnv_column_map,
  file.path(OUT_DIR, "gcsi_2019_cnv_column_map.csv")
)

fwrite(
  mutation_column_map,
  file.path(OUT_DIR, "gcsi_2019_mutation_column_map.csv")
)

# -------------------------------------------------------------------------
# pre_clinical_rna_seq.csv
# -------------------------------------------------------------------------

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

cat("Wrote GCSI RNA-seq assay CSV\n")

# -------------------------------------------------------------------------
# pre_clinical_copy_number_variation.csv
# -------------------------------------------------------------------------

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

cat("Wrote GCSI CNV assay CSV\n")

# -------------------------------------------------------------------------
# pre_clinical_mutation.csv
# -------------------------------------------------------------------------

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

cat("Wrote GCSI mutation assay CSV\n")

# -------------------------------------------------------------------------
# Cleanup
# -------------------------------------------------------------------------

cat("\nFreeing GCSI objects from memory\n")

rm(
  gcsi_2019,
  selected_sample_dt,
  sheet_target_data,
  selection_data,
  canonical_lookup,
  rnaseq_se,
  cnv_se,
  mutation_se,
  rnaseq_gene_map,
  cnv_gene_map,
  mutation_gene_map,
  rnaseq_column_map,
  cnv_column_map,
  mutation_column_map
)

gc(verbose = TRUE)

cat("Finished extracting selected GCSI preclinical CSVs into:", OUT_DIR, "\n")