suppressPackageStartupMessages({
  library(data.table)
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

PRISM_RDS_PATH <- "extraction/data/raw/preclinical/PSet_PRISM.rds"
SHEET_CELL_LINE_QC_PATH <- "extraction/data/raw/preclinical/All_PSets_sarcoma_cell_line_QC.csv"

OUT_DIR <- "extraction/data/proc/preclinical/PRISM"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Final sample IDs will look like:
#   PRISM_<PRISM.sampleid>
SAMPLE_ID_PREFIX <- "PRISM_"

# New extraction set:
#   1. Read external QC sheet.
#   2. Keep rows where mod_tissueid == "Soft Tissue".
#   3. Match sheet$cell_line against ONE resolved PRISM sample-slot cell-line column.
SHEET_TARGET_CELL_LINE_COL <- "cell_line"

# PRISM @sample$sampleid has been used as the cell-line / treatment-response key.
# Keep it first. Only one resolved column is used for the new extraction match.
PRISM_CELL_LINE_NAME_CANDIDATES <- c(
  "sampleid",
  "cell_line_name",
  "Cell.line.primary.name",
  "cellosaurus.cellLineName",
  "cellLineName",
  "cellline",
  "Name",
  "sample_rowname"
)

# Original criteria are now audit-only.
# They DO NOT control the final extraction set anymore.
TARGET_CELL_LINES_RAW <- paste0(
  "CHSA8926|GI-1|RD|Rh41|SKN|TE 125.T|105KC|GCT|KYM-1|Rh18|",
  "SKN|SW872|Aska-SS|GCT|Hs 819.T|OUMS-27|Rh30|SF539|SKN"
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

make_prefixed_sample_id <- function(sampleid) {
  sampleid <- clean_na(sampleid)
  out <- rep(NA_character_, length(sampleid))
  keep_idx <- !is.na(sampleid)

  out[keep_idx] <- ifelse(
    startsWith(sampleid[keep_idx], SAMPLE_ID_PREFIX),
    sampleid[keep_idx],
    paste0(SAMPLE_ID_PREFIX, sampleid[keep_idx])
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

resolve_prism_cell_line_col <- function(sample_dt) {
  col <- first_existing_col(
    sample_dt,
    PRISM_CELL_LINE_NAME_CANDIDATES
  )

  if (is.na(col)) {
    stop(
      "Could not find a PRISM sample-slot cell-line-name column. ",
      "Expected one of: ",
      paste(PRISM_CELL_LINE_NAME_CANDIDATES, collapse = ", "),
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
    file.path(out_dir, "prism_sheet_column_names.csv")
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
      file.path(out_dir, "prism_sheet_column_names.csv")
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
        "Cellosaurus.Accession.ID",
        "cellosaurus.cvcl_id",
        "cellosaurus.cellosaurus.cvcl_id",
        "cellosaurus.accession",
        "cellosaurus.accession.id",
        "accession",
        "cell_line_accession"
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
        "Cellosaurus.Category",
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
        "Cellosaurus.Sex",
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
        "Cellosaurus.Age",
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
      !is.na(sheet_dataset_value) & grepl("PRISM", toupper(sheet_dataset_value)),
      1L,
      fifelse(
        !is.na(sheet_object_type_value) & grepl("PRISM", toupper(sheet_object_type_value)),
        2L,
        99L
      )
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
    file.path(out_dir, "prism_sheet_soft_tissue_rows.csv")
  )

  fwrite(
    sheet_soft_tissue_metadata_dt,
    file.path(out_dir, "prism_sheet_soft_tissue_cell_line_metadata.csv")
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
    file.path(out_dir, "prism_sheet_soft_tissue_cell_line_targets.csv")
  )

  cat("QC-sheet soft tissue rows:", nrow(sheet_soft_tissue_dt), "\n")
  cat("QC-sheet unique soft tissue cell lines:", nrow(sheet_soft_tissue_metadata_dt), "\n")

  list(
    sheet_soft_tissue_dt = sheet_soft_tissue_dt,
    sheet_soft_tissue_metadata_dt = sheet_soft_tissue_metadata_dt
  )
}

# -------------------------------------------------------------------------
# PRISM sample selection and original-criteria audit
# -------------------------------------------------------------------------

build_prism_selection <- function(
  pset_prism,
  target_cell_lines,
  sheet_soft_tissue_metadata_dt,
  out_dir
) {
  sample_dt <- as.data.table(
    pset_prism@sample,
    keep.rownames = "sample_rowname"
  )

  fwrite(
    data.table(column_name = colnames(sample_dt)),
    file.path(out_dir, "prism_sample_slot_columns.csv")
  )

  prism_cell_line_col <- resolve_prism_cell_line_col(sample_dt)

  fwrite(
    data.table(
      prism_cell_line_name_column_used = prism_cell_line_col
    ),
    file.path(out_dir, "prism_cell_line_name_column_used.csv")
  )

  sample_dt[
    ,
    prism_cell_line_name := clean_na(get(prism_cell_line_col))
  ]

  sample_dt[
    ,
    normalized_prism_cell_line_name := normalize_cell_line_name(prism_cell_line_name)
  ]

  # New extraction criterion only:
  # QC sheet mod_tissueid == "Soft Tissue", matched by sheet$cell_line to the
  # one resolved PRISM sample-slot cell-line-name column.
  sample_dt[
    ,
    new_sheet_soft_tissue_match :=
      normalized_prism_cell_line_name %in%
        sheet_soft_tissue_metadata_dt$normalized_sheet_cell_line_name
  ]

  selected_sample_dt <- sample_dt[
    new_sheet_soft_tissue_match == TRUE
  ]

  if (nrow(selected_sample_dt) == 0) {
    stop(
      "No PRISM samples matched QC-sheet mod_tissueid == 'Soft Tissue' ",
      "cell lines using PRISM sample-slot column: ",
      prism_cell_line_col
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
    by.x = "normalized_prism_cell_line_name",
    by.y = "normalized_sheet_cell_line_name",
    all.x = TRUE
  )

  selected_sample_dt[
    ,
    canonical_cell_line_name := clean_na(sheet_cell_line_name)
  ]

  selected_sample_dt[
    is.na(canonical_cell_line_name),
    canonical_cell_line_name := prism_cell_line_name
  ]

  selected_sample_dt[
    ,
    final_cell_line_name := canonical_cell_line_name
  ]

  # For PRISM, the sample slot's sampleid is the cell-line / treatment-response key.
  selected_sample_dt[
    ,
    source_sampleid := clean_na(safe_col(.SD, c("sampleid", prism_cell_line_col), n = .N))
  ]

  selected_sample_dt[
    is.na(source_sampleid),
    source_sampleid := canonical_cell_line_name
  ]

  # For PRISM database sample ID, use PRISM.sampleid when present, otherwise the
  # resolved source_sampleid. Keep this separate from the treatment key.
  selected_sample_dt[
    ,
    canonical_sample_source_id := clean_na(safe_col(.SD, c("PRISM.sampleid"), n = .N))
  ]

  selected_sample_dt[
    is.na(canonical_sample_source_id),
    canonical_sample_source_id := source_sampleid
  ]

  selected_sample_dt[
    ,
    canonical_sample_id := make_prefixed_sample_id(canonical_sample_source_id)
  ]

  selected_sample_dt[
    ,
    treatment_response_sample_key := source_sampleid
  ]

  selected_sample_dt[
    ,
    normalized_treatment_response_sample_key := normalize_cell_line_name(treatment_response_sample_key)
  ]

  selected_sample_dt <- selected_sample_dt[
    !is.na(canonical_sample_source_id) &
      !is.na(canonical_sample_id) &
      !is.na(canonical_cell_line_name) &
      !is.na(treatment_response_sample_key)
  ]

  # -----------------------------------------------------------------------
  # Original criteria audit only.
  # These rows do not change selected_sample_dt.
  # -----------------------------------------------------------------------

  original_target_norm <- normalize_cell_line_name(target_cell_lines)

  sample_dt[, original_target_match := FALSE]

  if (has_col(sample_dt, "sampleid")) {
    sample_dt[
      normalize_cell_line_name(sampleid) %in% original_target_norm,
      original_target_match := TRUE
    ]
  }

  sample_dt[, original_soft_tissue_match := FALSE]

  if (has_col(sample_dt, "tissueid")) {
    sample_dt[
      normalize_category_value(tissueid) == "soft_tissue",
      original_soft_tissue_match := TRUE
    ]
  }

  if (has_col(sample_dt, "PRISM.tissueid")) {
    sample_dt[
      normalize_category_value(`PRISM.tissueid`) == "soft_tissue",
      original_soft_tissue_match := TRUE
    ]
  }

  if (has_col(sample_dt, "primary_tissue")) {
    sample_dt[
      normalize_category_value(primary_tissue) == "soft_tissue",
      original_soft_tissue_match := TRUE
    ]
  }

  original_criteria_dt <- sample_dt[
    original_target_match == TRUE |
      original_soft_tissue_match == TRUE
  ]

  original_criteria_dt[
    ,
    original_cell_line_name := prism_cell_line_name
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
      c("sampleid", prism_cell_line_col, "sample_rowname")
    )),
    source_prism_sampleid = clean_na(safe_col(original_criteria_dt, c("PRISM.sampleid"))),
    sample_id = make_prefixed_sample_id(clean_na(safe_col(original_criteria_dt, c("PRISM.sampleid", "sampleid", prism_cell_line_col)))),
    cell_line_name = clean_na(original_criteria_dt$original_cell_line_name),
    prism_cell_line_name_column_used = prism_cell_line_col,
    tissueid = clean_na(safe_col(original_criteria_dt, c("tissueid"))),
    PRISM.tissueid = clean_na(safe_col(original_criteria_dt, c("PRISM.tissueid"))),
    primary_tissue = clean_na(safe_col(original_criteria_dt, c("primary_tissue"))),
    original_target_match = original_criteria_dt$original_target_match,
    original_soft_tissue_match = original_criteria_dt$original_soft_tissue_match
  )

  original_not_in_final_audit_dt <- data.table(
    source_sampleid = clean_na(coalesce_dt_cols(
      original_not_in_final_dt,
      c("sampleid", prism_cell_line_col, "sample_rowname")
    )),
    source_prism_sampleid = clean_na(safe_col(original_not_in_final_dt, c("PRISM.sampleid"))),
    sample_id = make_prefixed_sample_id(clean_na(safe_col(original_not_in_final_dt, c("PRISM.sampleid", "sampleid", prism_cell_line_col)))),
    cell_line_name = clean_na(original_not_in_final_dt$original_cell_line_name),
    prism_cell_line_name_column_used = prism_cell_line_col,
    tissueid = clean_na(safe_col(original_not_in_final_dt, c("tissueid"))),
    PRISM.tissueid = clean_na(safe_col(original_not_in_final_dt, c("PRISM.tissueid"))),
    primary_tissue = clean_na(safe_col(original_not_in_final_dt, c("primary_tissue"))),
    original_target_match = original_not_in_final_dt$original_target_match,
    original_soft_tissue_match = original_not_in_final_dt$original_soft_tissue_match
  )

  fwrite(
    original_audit_dt,
    file.path(out_dir, "prism_original_criteria_cell_line_list_audit_only.csv")
  )

  fwrite(
    original_not_in_final_audit_dt,
    file.path(out_dir, "prism_original_criteria_not_in_final_new_criteria.csv")
  )

  fwrite(
    selected_sample_dt[
      ,
      .(
        sample_id = canonical_sample_id,
        source_prism_sampleid = canonical_sample_source_id,
        source_sampleid,
        cell_line_name = final_cell_line_name,
        prism_cell_line_name,
        prism_cell_line_name_column_used = prism_cell_line_col,
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
    file.path(out_dir, "prism_final_extracted_cell_line_list_new_criteria.csv")
  )

  fwrite(
    selected_sample_dt,
    file.path(out_dir, "prism_selected_sample_slot_rows.csv")
  )

  cat("PRISM cell-line-name column used:", prism_cell_line_col, "\n")
  cat("Final new-criteria selected PRISM sample rows:", nrow(selected_sample_dt), "\n")
  cat("Original-criteria audit rows:", nrow(original_criteria_dt), "\n")
  cat("Original-criteria rows not in final new criteria:", nrow(original_not_in_final_dt), "\n")

  list(
    selected_sample_dt = selected_sample_dt,
    original_criteria_dt = original_criteria_dt,
    original_not_in_final_dt = original_not_in_final_dt,
    prism_cell_line_col = prism_cell_line_col
  )
}

build_canonical_lookup_prism <- function(selected_sample_dt) {
  alt_cols <- c(
    # Most important for treatment-response sample/cell-line matching.
    "sampleid",
    "treatment_response_sample_key",
    "source_sampleid",

    # Final PRISM output IDs.
    "PRISM.sampleid",
    "canonical_sample_source_id",
    "canonical_sample_id",

    # Final/canonical cell-line names.
    "canonical_cell_line_name",
    "final_cell_line_name",
    "prism_cell_line_name",
    "sheet_cell_line_name",

    # Other possible aliases.
    "sample_rowname",
    "cellosaurus.cellLineName",
    "Cell.line.primary.name",
    "cell_line_name",
    "cellLineName",
    "cellline",
    "Name"
  )

  pieces <- list()

  for (col in alt_cols) {
    if (has_col(selected_sample_dt, col)) {
      pieces[[col]] <- data.table(
        lookup_name = clean_na(selected_sample_dt[[col]]),
        normalized_lookup_name = normalize_cell_line_name(selected_sample_dt[[col]]),

        # Final prefixed sample ID used in pre_clinical_sample.csv.
        sample_id = selected_sample_dt$canonical_sample_id,

        # Original PRISM.sampleid used to construct final sample ID.
        source_prism_sampleid = selected_sample_dt$canonical_sample_source_id,

        # sampleid / resolved sample-slot cell-line key.
        source_sampleid = selected_sample_dt$treatment_response_sample_key,

        cell_line_name = selected_sample_dt$canonical_cell_line_name,

        source_column = col
      )
    }
  }

  lookup_dt <- rbindlist(pieces, fill = TRUE)

  lookup_dt <- lookup_dt[
    !is.na(normalized_lookup_name) &
      !is.na(sample_id) &
      !is.na(source_prism_sampleid) &
      !is.na(source_sampleid) &
      !is.na(cell_line_name)
  ]

  unique(
    lookup_dt,
    by = c("normalized_lookup_name", "sample_id", "source_column")
  )
}

# -------------------------------------------------------------------------
# PRISM treatment response using summarizeSensitivityProfiles
# -------------------------------------------------------------------------

get_available_cell_lines_for_sensitivity <- function(pset_prism) {
  out <- tryCatch(
    cellNames(pset_prism),
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
    pset_prism@sample,
    keep.rownames = "sample_rowname"
  )

  candidate_col <- first_existing_col(
    sample_dt,
    PRISM_CELL_LINE_NAME_CANDIDATES
  )

  if (is.na(candidate_col)) {
    stop(
      "Could not determine available PRISM cell lines. ",
      "cellNames(pset_prism) failed and no usable sample column was found."
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
      paste0("prism_treatment_ids_from_", sensitivity_measure, ".csv")
    )
  )

  fwrite(
    data.table(
      cell_line_name_for_summary = colnames(mat)
    ),
    file.path(
      out_dir,
      paste0("prism_cell_lines_from_", sensitivity_measure, ".csv")
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

extract_treatment_response_prism <- function(pset_prism, canonical_lookup, out_path) {
  selected_cell_line_lookup <- unique(
    canonical_lookup[
      source_column %in% c(
        "canonical_cell_line_name",
        "final_cell_line_name",
        "treatment_response_sample_key",
        "source_sampleid",
        "sampleid",
        "prism_cell_line_name"
      ),
      .(
        sample_id,
        source_prism_sampleid,
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

  available_cell_lines <- get_available_cell_lines_for_sensitivity(pset_prism)

  available_lookup <- data.table(
    prism_cell_line_name_for_summary = available_cell_lines,
    normalized_cell_line_name = normalize_cell_line_name(available_cell_lines)
  )

  available_lookup <- unique(
    available_lookup[
      !is.na(prism_cell_line_name_for_summary) &
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
    sub("\\.csv$", "_selected_cell_lines_available_in_prism.csv", out_path)
  )

  missing_cell_lines <- cell_line_match[is.na(prism_cell_line_name_for_summary)]

  fwrite(
    missing_cell_lines,
    sub("\\.csv$", "_selected_cell_lines_missing_from_prism.csv", out_path)
  )

  cell_lines_to_use <- unique(clean_na(cell_line_match$prism_cell_line_name_for_summary))
  cell_lines_to_use <- cell_lines_to_use[!is.na(cell_lines_to_use)]

  if (length(cell_lines_to_use) == 0) {
    stop(
      "None of the selected PRISM cell_line_name values matched sensitivity cell lines. ",
      "See: ",
      sub("\\.csv$", "_selected_cell_lines_missing_from_prism.csv", out_path)
    )
  }

  aac_long <- summarized_sensitivity_to_long(
    pset = pset_prism,
    sensitivity_measure = "aac_recomputed",
    cell_lines = cell_lines_to_use,
    out_dir = OUT_DIR
  )

  ic50_long <- summarized_sensitivity_to_long(
    pset = pset_prism,
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
      !is.na(prism_cell_line_name_for_summary),
      .(
        normalized_cell_line_name,
        prism_cell_line_name_for_summary,
        sample_id,
        source_prism_sampleid,
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

  cat("Wrote PRISM summarized treatment response rows:", nrow(treatment_response_dt), "\n")
}

# -------------------------------------------------------------------------
# Load PRISM
# -------------------------------------------------------------------------

cat("\n==============================\n")
cat("Loading PRISM\n")
cat("==============================\n")

prism <- read_updated_rds(PRISM_RDS_PATH)

# -------------------------------------------------------------------------
# Select samples/cell lines using external QC sheet new criteria only
# -------------------------------------------------------------------------

sheet_target_data <- read_qc_sheet_soft_tissue_targets(
  path = SHEET_CELL_LINE_QC_PATH,
  out_dir = OUT_DIR
)

selection_data <- build_prism_selection(
  pset_prism = prism,
  target_cell_lines = TARGET_CELL_LINES,
  sheet_soft_tissue_metadata_dt = sheet_target_data$sheet_soft_tissue_metadata_dt,
  out_dir = OUT_DIR
)

selected_sample_dt <- selection_data$selected_sample_dt

canonical_lookup <- build_canonical_lookup_prism(selected_sample_dt)

cat("Final new-criteria selected sample rows:", nrow(selected_sample_dt), "\n")
cat("Canonical lookup rows:", nrow(canonical_lookup), "\n")

fwrite(
  canonical_lookup,
  file.path(OUT_DIR, "prism_canonical_lookup.csv")
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
    c(
      "primary_tissue",
      "PRISM.tissueid",
      "tissueid",
      "ccle_primary_site",
      "site_primary",
      "Site.Primary",
      "primary_site"
    )
  )),

  site_subtype1 = clean_na(safe_col(
    selected_sample_dt,
    c(
      "site_subtype1",
      "Site.Subtype1",
      "Site_Subtype1"
    )
  )),

  site_subtype2 = clean_na(safe_col(
    selected_sample_dt,
    c(
      "site_subtype2",
      "Site.Subtype2",
      "Site_Subtype2"
    )
  )),

  site_subtype3 = clean_na(safe_col(
    selected_sample_dt,
    c(
      "site_subtype3",
      "Site.Subtype3",
      "Site_Subtype3"
    )
  )),

  histology = clean_na(safe_col(
    selected_sample_dt,
    c(
      "ccle_primary_hist",
      "histology",
      "Histology"
    )
  )),

  histology_subtype1 = clean_na(safe_col(
    selected_sample_dt,
    c(
      "ccle_hist_subtype_1",
      "histology_subtype1",
      "Hist.Subtype1",
      "Histology_Subtype1"
    )
  )),

  histology_subtype2 = clean_na(safe_col(
    selected_sample_dt,
    c(
      "ccle_hist_subtype_2",
      "histology_subtype2",
      "Hist.Subtype2",
      "Histology_Subtype2"
    )
  )),

  histology_subtype3 = clean_na(safe_col(
    selected_sample_dt,
    c(
      "ccle_hist_subtype_3",
      "histology_subtype3",
      "Hist.Subtype3",
      "Histology_Subtype3"
    )
  )),

  gender = clean_na(safe_col(
    selected_sample_dt,
    c(
      "gender",
      "Gender",
      "sex",
      "Sex"
    )
  )),

  age = parse_age_int(safe_col(
    selected_sample_dt,
    c(
      "age",
      "Age"
    )
  )),

  race = clean_na(safe_col(
    selected_sample_dt,
    c(
      "race",
      "Race"
    )
  )),

  diseases = clean_na(safe_col(
    selected_sample_dt,
    c(
      "Cellosaurus.Disease.Type",
      "diseases",
      "cellosaurus.diseases",
      "Disease",
      "disease"
    )
  )),

  disease_type = clean_na(safe_col(
    selected_sample_dt,
    c(
      "Cellosaurus.Disease.Type",
      "disease_type",
      "type",
      "PRISM.tissueid",
      "tissueid",
      "primary_tissue"
    )
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
# pre_clinical_treatment_response.csv using summarizeSensitivityProfiles
# -------------------------------------------------------------------------

extract_treatment_response_prism(
  pset_prism = prism,
  canonical_lookup = canonical_lookup,
  out_path = file.path(OUT_DIR, "pre_clinical_treatment_response.csv")
)

# -------------------------------------------------------------------------
# Cleanup
# -------------------------------------------------------------------------

cat("\nFreeing PRISM objects from memory\n")

rm(
  prism,
  selected_sample_dt,
  sheet_target_data,
  selection_data,
  canonical_lookup,
  cell_line_dt,
  sample_out_dt
)

gc(verbose = TRUE)

cat("Finished extracting selected PRISM preclinical CSVs into:", OUT_DIR, "\n")