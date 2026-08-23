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

# Prefix used only for database/output sample IDs.
# Raw CCLE sample IDs are preserved separately as source_sampleid for matching
# CCLE_2019 sample IDs to CCLE_2015 microarray colData(rna)$CCLE.name.
SAMPLE_ID_PREFIX <- "CCLE_"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# The new extraction set comes ONLY from the QC sheet where:
#   mod_tissueid == "Soft Tissue"
# and then matching sheet$cell_line against the CCLE_2019 sample slot
# resolved cell-line-name column.
SHEET_TARGET_CELL_LINE_COL <- "cell_line"

# Prefer an explicit cell_line_name column if it exists in CCLE_2019@sample.
# Fall back to cellosaurus.cellLineName because the previous CCLE script used it.
# Only ONE resolved column is used for the actual new extraction match.
CCLE_2019_CELL_LINE_NAME_CANDIDATES <- c(
  "cell_line_name",
  "cellosaurus.cellLineName"
)

# Original criteria are now audit-only.
# They do NOT control the final extraction set anymore.
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

prefix_sample_id <- function(x, prefix = SAMPLE_ID_PREFIX) {
  x <- clean_na(x)

  out <- rep(NA_character_, length(x))
  keep_idx <- !is.na(x)

  out[keep_idx] <- ifelse(
    startsWith(x[keep_idx], prefix),
    x[keep_idx],
    paste0(prefix, x[keep_idx])
  )

  out
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

resolve_ccle_2019_cell_line_col <- function(sample_dt) {
  col <- first_existing_col(
    sample_dt,
    CCLE_2019_CELL_LINE_NAME_CANDIDATES
  )

  if (is.na(col)) {
    stop(
      "Could not find a CCLE_2019 sample-slot cell-line-name column. ",
      "Expected one of: ",
      paste(CCLE_2019_CELL_LINE_NAME_CANDIDATES, collapse = ", "),
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
    file.path(out_dir, "ccle_sheet_column_names.csv")
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
      file.path(out_dir, "ccle_sheet_column_names.csv")
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

  sheet_soft_tissue_dt[
    ,
    tissueid := clean_na(tissueid)
  ]

  sheet_soft_tissue_dt[
    ,
    mod_tissueid := clean_na(mod_tissueid)
  ]

  sheet_soft_tissue_dt[
    ,
    accession := clean_na(coalesce_dt_cols(
      .SD,
      c(
        "Cellosaurus.Accession.id",
        "cellosaurus.cvcl_id",
        "cellosaurus.cellosaurus.cvcl_id"
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
        "CellLine.Type"
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
        "Gender",
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
        "Age",
        "age"
      )
    ))
  ]

  sheet_soft_tissue_dt[
    ,
    metadata_priority := fifelse(
      dataset == "CCLE_2019",
      1L,
      fifelse(dataset == "CCLE_2015", 2L, 99L)
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
    file.path(out_dir, "ccle_sheet_soft_tissue_rows.csv")
  )

  fwrite(
    sheet_soft_tissue_metadata_dt,
    file.path(out_dir, "ccle_sheet_soft_tissue_cell_line_metadata.csv")
  )

  fwrite(
    unique(
      sheet_soft_tissue_metadata_dt[
        ,
        .(
          sheet_cell_line_name,
          normalized_sheet_cell_line_name,
          dataset,
          object_type,
          sampleid,
          tissueid,
          mod_tissueid,
          accession,
          category,
          sex,
          age_raw
        )
      ]
    ),
    file.path(out_dir, "ccle_sheet_soft_tissue_cell_line_targets.csv")
  )

  cat("QC-sheet soft tissue rows:", nrow(sheet_soft_tissue_dt), "\n")
  cat("QC-sheet unique soft tissue cell lines:", nrow(sheet_soft_tissue_metadata_dt), "\n")

  list(
    sheet_soft_tissue_dt = sheet_soft_tissue_dt,
    sheet_soft_tissue_metadata_dt = sheet_soft_tissue_metadata_dt
  )
}

# -------------------------------------------------------------------------
# CCLE 2019 selection and audit helpers
# -------------------------------------------------------------------------

build_ccle_2019_selection <- function(
  pset_2019,
  target_cell_lines,
  sheet_soft_tissue_metadata_dt,
  out_dir
) {
  sample_dt <- as.data.table(
    pset_2019@sample,
    keep.rownames = "sample_rowname"
  )

  ccle_cell_line_col <- resolve_ccle_2019_cell_line_col(sample_dt)

  fwrite(
    data.table(
      ccle_2019_cell_line_name_column_used = ccle_cell_line_col
    ),
    file.path(out_dir, "ccle_2019_cell_line_name_column_used.csv")
  )

  sample_dt[
    ,
    ccle_2019_cell_line_name := clean_na(get(ccle_cell_line_col))
  ]

  sample_dt[
    ,
    normalized_ccle_2019_cell_line_name :=
      normalize_cell_line_name(ccle_2019_cell_line_name)
  ]

  # -----------------------------------------------------------------------
  # New extraction criterion:
  # QC sheet mod_tissueid == "Soft Tissue" cell lines only.
  # Match sheet$cell_line to ONE resolved CCLE_2019 sample-slot cell-line-name column.
  # -----------------------------------------------------------------------

  sample_dt[
    ,
    new_sheet_soft_tissue_match :=
      normalized_ccle_2019_cell_line_name %in%
        sheet_soft_tissue_metadata_dt$normalized_sheet_cell_line_name
  ]

  selected_sample_dt <- sample_dt[
    new_sheet_soft_tissue_match == TRUE
  ]

  if (nrow(selected_sample_dt) == 0) {
    stop(
      "No CCLE_2019 samples matched QC-sheet mod_tissueid == 'Soft Tissue' ",
      "cell lines using CCLE_2019 sample-slot column: ",
      ccle_cell_line_col
    )
  }

  selected_sample_dt <- merge(
    selected_sample_dt,
    sheet_soft_tissue_metadata_dt[
      ,
      .(
        normalized_sheet_cell_line_name,
        sheet_cell_line_name,
        sheet_dataset = dataset,
        sheet_object_type = object_type,
        sheet_sampleid = sampleid,
        sheet_tissueid = tissueid,
        sheet_mod_tissueid = mod_tissueid,
        sheet_accession = accession,
        sheet_category = category,
        sheet_sex = sex,
        sheet_age = age_raw
      )
    ],
    by.x = "normalized_ccle_2019_cell_line_name",
    by.y = "normalized_sheet_cell_line_name",
    all.x = TRUE
  )

  selected_sample_dt[
    ,
    final_cell_line_name := clean_na(sheet_cell_line_name)
  ]

  selected_sample_dt[
    is.na(final_cell_line_name),
    final_cell_line_name := ccle_2019_cell_line_name
  ]

  if (!has_col(selected_sample_dt, "sampleid")) {
    stop(
      "CCLE_2019 sample slot must contain sampleid because this is the ",
      "final pre_clinical_sample.id and the bridge to CCLE_2015 microarray colData$CCLE.name."
    )
  }

  selected_sample_dt[
    ,
    source_sampleid := clean_na(sampleid)
  ]

  selected_sample_dt[
    ,
    prefixed_sampleid := prefix_sample_id(source_sampleid)
  ]

  # -----------------------------------------------------------------------
  # Original criteria are audit-only now.
  # They do not affect selected_sample_dt.
  # -----------------------------------------------------------------------

  original_target_norm <- normalize_cell_line_name(target_cell_lines)

  original_candidate_cols <- c(
    "cellosaurus.cellLineName",
    "CCLE.name",
    "CCLE.sampleid",
    "sampleid",
    "unique.sampleid",
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

  if (has_col(sample_dt, "CCLE.site_Primary")) {
    sample_dt[
      normalize_category_value(`CCLE.site_Primary`) == "soft_tissue",
      original_soft_tissue_match := TRUE
    ]
  }

  if (has_col(sample_dt, "CCLE.type")) {
    sample_dt[
      normalize_category_value(`CCLE.type`) == "soft_tissue",
      original_soft_tissue_match := TRUE
    ]
  }

  original_criteria_dt <- sample_dt[
    original_target_match == TRUE |
      original_soft_tissue_match == TRUE
  ]

  original_criteria_dt[
    ,
    original_cell_line_name := ccle_2019_cell_line_name
  ]

  original_criteria_dt[
    ,
    normalized_original_cell_line_name :=
      normalize_cell_line_name(original_cell_line_name)
  ]

  final_norm <- unique(
    normalize_cell_line_name(selected_sample_dt$final_cell_line_name)
  )

  original_not_in_final_dt <- original_criteria_dt[
    !(normalized_original_cell_line_name %in% final_norm)
  ]

  fwrite(
    original_criteria_dt[
      ,
      .(
        sampleid = prefix_sample_id(sampleid),
        source_sampleid = clean_na(sampleid),
        cell_line_name = original_cell_line_name,
        ccle_2019_cell_line_name_column_used = ccle_cell_line_col,
        CCLE.name = if (has_col(.SD, "CCLE.name")) clean_na(`CCLE.name`) else NA_character_,
        CCLE.sampleid = if (has_col(.SD, "CCLE.sampleid")) clean_na(`CCLE.sampleid`) else NA_character_,
        CCLE.site_Primary = if (has_col(.SD, "CCLE.site_Primary")) clean_na(`CCLE.site_Primary`) else NA_character_,
        CCLE.type = if (has_col(.SD, "CCLE.type")) clean_na(`CCLE.type`) else NA_character_,
        original_target_match,
        original_soft_tissue_match
      )
    ],
    file.path(out_dir, "ccle_original_criteria_cell_line_list_audit_only.csv")
  )

  fwrite(
    original_not_in_final_dt[
      ,
      .(
        sampleid = prefix_sample_id(sampleid),
        source_sampleid = clean_na(sampleid),
        cell_line_name = original_cell_line_name,
        ccle_2019_cell_line_name_column_used = ccle_cell_line_col,
        CCLE.name = if (has_col(.SD, "CCLE.name")) clean_na(`CCLE.name`) else NA_character_,
        CCLE.sampleid = if (has_col(.SD, "CCLE.sampleid")) clean_na(`CCLE.sampleid`) else NA_character_,
        CCLE.site_Primary = if (has_col(.SD, "CCLE.site_Primary")) clean_na(`CCLE.site_Primary`) else NA_character_,
        CCLE.type = if (has_col(.SD, "CCLE.type")) clean_na(`CCLE.type`) else NA_character_,
        original_target_match,
        original_soft_tissue_match
      )
    ],
    file.path(out_dir, "ccle_original_criteria_not_in_final_new_criteria.csv")
  )

  fwrite(
    selected_sample_dt[
      ,
      .(
        sampleid = prefix_sample_id(source_sampleid),
        source_sampleid = clean_na(source_sampleid),
        cell_line_name = final_cell_line_name,
        ccle_2019_cell_line_name = ccle_2019_cell_line_name,
        ccle_2019_cell_line_name_column_used = ccle_cell_line_col,
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
    file.path(out_dir, "ccle_final_extracted_cell_line_list_new_criteria.csv")
  )

  cat("CCLE_2019 cell-line-name column used:", ccle_cell_line_col, "\n")
  cat("Final new-criteria selected CCLE_2019 sample rows:", nrow(selected_sample_dt), "\n")
  cat("Original-criteria audit rows:", nrow(original_criteria_dt), "\n")
  cat("Original-criteria rows not in final new criteria:", nrow(original_not_in_final_dt), "\n")

  list(
    selected_sample_dt = selected_sample_dt,
    original_criteria_dt = original_criteria_dt,
    original_not_in_final_dt = original_not_in_final_dt,
    ccle_cell_line_col = ccle_cell_line_col
  )
}

build_canonical_lookup <- function(selected_sample_dt) {
  canonical_source_sample_id <- clean_na(selected_sample_dt[["source_sampleid"]])

  if (all(is.na(canonical_source_sample_id))) {
    canonical_source_sample_id <- clean_na(selected_sample_dt[["sampleid"]])
  }

  canonical_sample_id <- prefix_sample_id(canonical_source_sample_id)
  canonical_cell_line <- clean_na(selected_sample_dt[["final_cell_line_name"]])

  alt_cols <- c(
    "ccle_2019_cell_line_name",
    "final_cell_line_name",
    "CCLE.name",
    "CCLE.sampleid",
    "source_sampleid",
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
        source_sampleid = canonical_source_sample_id,
        cell_line_name = canonical_cell_line,
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

  fwrite(
    rna_cd,
    file.path(OUT_DIR, "ccle_2015_microarray_coldata.csv")
  )

  fwrite(
    data.table(column_name = colnames(rna_cd)),
    file.path(OUT_DIR, "ccle_2015_microarray_coldata_columns.csv")
  )

  if (!has_col(rna_cd, "CCLE.name")) {
    stop(
      "CCLE_2015 RNA colData must contain 'CCLE.name'. ",
      "The requested mapping is CCLE_2019@sample$sampleid -> colData(rna)$CCLE.name. ",
      "See: ",
      file.path(OUT_DIR, "ccle_2015_microarray_coldata_columns.csv")
    )
  }

  canonical_sample_lookup <- unique(
    canonical_lookup[
      source_column %in% c("source_sampleid", "sampleid"),
      .(
        sample_id,
        source_sampleid,
        normalized_source_sampleid = normalize_cell_line_name(source_sampleid),
        cell_line_name
      )
    ],
    by = "source_sampleid"
  )

  canonical_sample_lookup <- canonical_sample_lookup[
    !is.na(sample_id) &
      !is.na(source_sampleid) &
      !is.na(normalized_source_sampleid) &
      !is.na(cell_line_name)
  ]

  rna_candidates <- data.table(
    colname = rna_cd$colname,
    coldata_rownames = rna_cd$coldata_rownames,
    rna_coldata_ccle_name = clean_na(rna_cd[["CCLE.name"]]),
    normalized_source_sampleid = normalize_cell_line_name(rna_cd[["CCLE.name"]])
  )

  direct_mapped <- merge(
    rna_candidates,
    canonical_sample_lookup,
    by = "normalized_source_sampleid",
    all = FALSE
  )

  if (nrow(direct_mapped) == 0) {
    stop(
      "No CCLE 2015 RNA microarray columns matched using raw ",
      "CCLE_2019@sample$sampleid -> colData(rna)$CCLE.name. ",
      "The database sample_id is prefixed, but matching still uses source_sampleid. ",
      "Check ccle_2015_microarray_coldata.csv and ccle_2019_canonical_lookup.csv."
    )
  }

  direct_mapped[
    ,
    matched_column := "raw CCLE_2019@sample$sampleid -> colData(rna)$CCLE.name; output sample_id is prefixed"
  ]

  duplicate_sample_map <- direct_mapped[, .N, by = sample_id][N > 1]

  fwrite(
    duplicate_sample_map,
    file.path(OUT_DIR, "ccle_2015_microarray_duplicate_sample_columns.csv")
  )

  if (nrow(duplicate_sample_map) > 0) {
    warning(
      "CCLE 2015 RNA column map has duplicate source columns for ",
      nrow(duplicate_sample_map),
      " canonical sample_id values. Keeping first per sample_id. See: ",
      file.path(OUT_DIR, "ccle_2015_microarray_duplicate_sample_columns.csv")
    )
  }

  direct_mapped <- direct_mapped[
    colname %in% colnames(rna_se)
  ]

  direct_mapped <- unique(direct_mapped, by = "sample_id")

  direct_mapped[
    ,
    .(
      colname,
      sample_id,
      source_sampleid,
      cell_line_name,
      matched_column,
      rna_coldata_ccle_name,
      coldata_rownames
    )
  ]
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
# CCLE 2015 treatment response using summarizeSensitivityProfiles

# -------------------------------------------------------------------------
# Treatment CID helpers
# -------------------------------------------------------------------------

build_treatment_cid_lookup <- function(pset, out_dir = OUT_DIR, label = "dataset") {
  empty_lookup <- data.table(
    treatment_id = character(),
    cid = character(),
    treatment_id_source_column = character()
  )

  treatment_obj <- tryCatch(
    pset@treatment,
    error = function(e) NULL
  )

  if (is.null(treatment_obj)) {
    warning("No treatment slot found for ", label, "; cid will be NA.")
    return(empty_lookup)
  }

  treatment_dt <- as.data.table(
    treatment_obj,
    keep.rownames = "treatment_rowname"
  )

  fwrite(
    data.table(column_name = colnames(treatment_dt)),
    file.path(out_dir, paste0(tolower(label), "_treatment_slot_columns.csv"))
  )

  if (nrow(treatment_dt) == 0) {
    warning("Treatment slot is empty for ", label, "; cid will be NA.")
    return(empty_lookup)
  }

  cid_col <- first_existing_col(
    treatment_dt,
    c(
      "cid",
      "CID",
      "PubChem.CID",
      "PubChem_CID",
      "pubchem_cid",
      "pubchem.cid",
      "compound_cid",
      "compound.cid"
    )
  )

  if (is.na(cid_col)) {
    warning(
      "No cid column found in the treatment slot for ", label,
      ". Wrote treatment slot columns for inspection; cid will be NA."
    )
    return(empty_lookup)
  }

  treatment_id_candidate_cols <- c(
    "treatment_rowname",
    "treatment_id",
    "treatmentid",
    "treatment.id",
    "TreatmentID",
    "Treatment.ID",
    "drug_id",
    "drugid",
    "drug.id",
    "compound_id",
    "compoundid",
    "compound.id",
    "master_cpd_id",
    "cpd_name",
    "drug_name",
    "treatment_name",
    "name"
  )

  pieces <- list()

  for (col in treatment_id_candidate_cols) {
    if (has_col(treatment_dt, col)) {
      pieces[[col]] <- data.table(
        treatment_id = clean_na(treatment_dt[[col]]),
        cid = clean_na(treatment_dt[[cid_col]]),
        treatment_id_source_column = col
      )
    }
  }

  if (length(pieces) == 0) {
    warning("No treatment ID candidate columns found for ", label, "; cid will be NA.")
    return(empty_lookup)
  }

  lookup <- rbindlist(pieces, fill = TRUE)
  lookup <- lookup[!is.na(treatment_id) & treatment_id != ""]
  lookup <- unique(lookup, by = "treatment_id")

  fwrite(
    lookup,
    file.path(out_dir, paste0(tolower(label), "_treatment_cid_lookup.csv"))
  )

  lookup
}

add_treatment_cid <- function(treatment_response_dt, pset, out_path, label = "dataset") {
  cid_lookup <- build_treatment_cid_lookup(pset, out_dir = OUT_DIR, label = label)

  treatment_response_dt <- copy(treatment_response_dt)
  treatment_response_dt[, treatment_id := clean_na(treatment_id)]

  if (nrow(cid_lookup) == 0) {
    treatment_response_dt[, cid := NA_character_]
    return(treatment_response_dt)
  }

  treatment_response_dt <- merge(
    treatment_response_dt,
    cid_lookup[, .(treatment_id, cid)],
    by = "treatment_id",
    all.x = TRUE
  )

  setcolorder(
    treatment_response_dt,
    c(
      "cell_line_name",
      "treatment_id",
      "cid",
      setdiff(colnames(treatment_response_dt), c("cell_line_name", "treatment_id", "cid"))
    )
  )

  fwrite(
    treatment_response_dt,
    sub("\\.csv$", "_with_cid_mapped_raw_rows.csv", out_path)
  )

  treatment_response_dt
}

# -------------------------------------------------------------------------

get_available_cell_lines_for_sensitivity <- function(pset_2015) {
  out <- tryCatch(
    cellNames(pset_2015),
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
    pset_2015@sample,
    keep.rownames = "sample_rowname"
  )

  candidate_col <- first_existing_col(
    sample_dt,
    c(
      "Cell.line.primary.name",
      "cellosaurus.cellLineName",
      "CCLE.name",
      "sampleid",
      "sample_rowname"
    )
  )

  if (is.na(candidate_col)) {
    stop(
      "Could not determine available CCLE 2015 cell lines. ",
      "cellNames(pset_2015) failed and no usable sample column was found."
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
      cell_line_name_for_summary = colnames(mat)
    ),
    file.path(
      out_dir,
      paste0("ccle_2015_cell_lines_from_", sensitivity_measure, ".csv")
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

extract_treatment_response_2015 <- function(pset_2015, canonical_lookup, out_path) {
  # Important:
  # Treatment response uses selected CCLE_2019 final cell_line_name values,
  # not CCLE_2019 sampleid values.
  selected_cell_line_lookup <- unique(
    canonical_lookup[
      source_column == "final_cell_line_name",
      .(
        canonical_sample_id = sample_id,
        cell_line_name,
        normalized_cell_line_name = normalize_cell_line_name(cell_line_name)
      )
    ],
    by = c("canonical_sample_id", "cell_line_name")
  )

  selected_cell_line_lookup <- selected_cell_line_lookup[
    !is.na(canonical_sample_id) &
      !is.na(cell_line_name) &
      !is.na(normalized_cell_line_name)
  ]

  if (nrow(selected_cell_line_lookup) == 0) {
    stop(
      "No selected cell_line_name values found in canonical lookup. ",
      "Treatment response requires selected CCLE_2019 final_cell_line_name values."
    )
  }

  available_cell_lines <- get_available_cell_lines_for_sensitivity(pset_2015)

  available_lookup <- data.table(
    ccle_2015_cell_line_name = available_cell_lines,
    normalized_cell_line_name = normalize_cell_line_name(available_cell_lines)
  )

  available_lookup <- unique(
    available_lookup[
      !is.na(ccle_2015_cell_line_name) &
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
    sub("\\.csv$", "_selected_cell_lines_available_in_ccle_2015.csv", out_path)
  )

  missing_cell_lines <- cell_line_match[is.na(ccle_2015_cell_line_name)]

  fwrite(
    missing_cell_lines,
    sub("\\.csv$", "_selected_cell_lines_missing_from_ccle_2015.csv", out_path)
  )

  cell_lines_to_use <- unique(clean_na(cell_line_match$ccle_2015_cell_line_name))
  cell_lines_to_use <- cell_lines_to_use[!is.na(cell_lines_to_use)]

  if (length(cell_lines_to_use) == 0) {
    stop(
      "None of the selected CCLE_2019 cell_line_name values matched CCLE_2015 sensitivity cell lines. ",
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
    by = c("treatment_id", "cell_line_name_for_summary"),
    all = TRUE
  )

  merged_response[
    ,
    normalized_cell_line_name := normalize_cell_line_name(cell_line_name_for_summary)
  ]

  response_cell_line_lookup <- unique(
    cell_line_match[
      !is.na(ccle_2015_cell_line_name),
      .(
        normalized_cell_line_name,
        ccle_2015_cell_line_name,
        canonical_sample_id,
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

  treatment_response_dt <- add_treatment_cid(
    treatment_response_dt = treatment_response_dt,
    pset = pset_2015,
    out_path = out_path,
    label = "CCLE"
  )

  # summarizeSensitivityProfiles(summary.stat = "mean") should already collapse
  # replicate/raw experiment duplicates. This check only catches accidental
  # duplicate mappings caused by duplicated selected cell-line names.
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
      "cid",
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

sheet_target_data <- read_qc_sheet_soft_tissue_targets(
  path = SHEET_CELL_LINE_QC_PATH,
  out_dir = OUT_DIR
)

selection_data <- build_ccle_2019_selection(
  pset_2019 = ccle_2019,
  target_cell_lines = TARGET_CELL_LINES,
  sheet_soft_tissue_metadata_dt = sheet_target_data$sheet_soft_tissue_metadata_dt,
  out_dir = OUT_DIR
)

selected_sample_dt <- selection_data$selected_sample_dt

canonical_lookup <- build_canonical_lookup(selected_sample_dt)

cat("Final new-criteria selected sample rows:", nrow(selected_sample_dt), "\n")
cat("Canonical lookup rows:", nrow(canonical_lookup), "\n")

canonical_lookup_path <- file.path(OUT_DIR, "ccle_2019_canonical_lookup.csv")

fwrite(
  canonical_lookup,
  canonical_lookup_path
)

cat("Wrote canonical lookup:", canonical_lookup_path, "\n")

# -------------------------------------------------------------------------
# pre_clinical_cell_line.csv from QC CSV metadata
# -------------------------------------------------------------------------

cell_line_dt <- data.table(
  cell_line_name = clean_na(selected_sample_dt[["final_cell_line_name"]]),
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

rm(cell_line_dt)
gc(verbose = FALSE)

# -------------------------------------------------------------------------
# pre_clinical_sample.csv from CCLE 2019
# -------------------------------------------------------------------------

sample_out_dt <- data.table(
  id = prefix_sample_id(selected_sample_dt[["source_sampleid"]]),
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
  value_col = "value",
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
  selection_data,
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
  value_col = "value",
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
