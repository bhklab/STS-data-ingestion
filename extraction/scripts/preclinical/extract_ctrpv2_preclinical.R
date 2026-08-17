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

CTRPV2_RDS_PATH <- "extraction/data/raw/preclinical/Pset_CTRPv2.rds"

OUT_DIR <- "extraction/data/proc/preclinical/CTRPv2"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Final output sample IDs will look like:
#   CTRP_22Rv1
#   CTRP_CAS1
SAMPLE_ID_PREFIX <- "CTRP_"

TARGET_CELL_LINES_RAW <- paste0(
  "Aska-SS|CS-1 [Human chondrosarcoma]|DM-3|NCI-H2731|RS-5|SF539|",
  "SW872|Yamato-SS|CHSA0011|D-247MG|Hs 729.T|RD|SF539|SNU-685|",
  "TE 159.T|Yamato-SS|CHSA0108|H-EMC-SS|MES-SA|NCI-H2596|Rh18|",
  "S-117|SK-UT-1|SW872|TE 159.T|Yamato-SS|CAL-78|CHSA8926"
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
  out <- paste0(SAMPLE_ID_PREFIX, sampleid)
  out[is.na(sampleid)] <- NA_character_
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

find_matching_key_col <- function(dt, profile_ids, candidates) {
  profile_ids <- as.character(profile_ids)

  for (candidate in candidates) {
    if (!has_col(dt, candidate)) {
      next
    }

    vals <- as.character(dt[[candidate]])
    n_matches <- sum(vals %in% profile_ids, na.rm = TRUE)

    if (n_matches > 0) {
      return(candidate)
    }
  }

  NA_character_
}

# -------------------------------------------------------------------------
# CTRPv2 sample selection
# -------------------------------------------------------------------------

build_selected_sample_table_ctrpv2 <- function(pset_ctrpv2, target_cell_lines) {
  sample_dt <- as.data.table(
    pset_ctrpv2@sample,
    keep.rownames = "sample_rowname"
  )

  fwrite(
    data.table(column_name = colnames(sample_dt)),
    file.path(OUT_DIR, "ctrpv2_sample_slot_columns.csv")
  )

  target_norm <- normalize_cell_line_name(target_cell_lines)

  cell_line_candidate_cols <- c(
    "sampleid",
    "sample_rowname",
    "ccl_name",
    "cellosaurus.cellLineName",
    "Cell.line.primary.name",
    "cell_line_name",
    "cellLineName",
    "cellline",
    "Name"
  )

  sample_dt[, target_match := FALSE]

  for (col in cell_line_candidate_cols) {
    if (has_col(sample_dt, col)) {
      sample_dt[
        normalize_cell_line_name(get(col)) %in% target_norm,
        target_match := TRUE
      ]
    }
  }

  sample_dt[, soft_tissue_match := FALSE]

  if (has_col(sample_dt, "tissueid")) {
    sample_dt[
      normalize_category_value(tissueid) == "soft_tissue",
      soft_tissue_match := TRUE
    ]
  } else {
    warning(
      "Column 'tissueid' was not found in CTRPv2 sample slot. ",
      "Soft Tissue matching from tissueid was skipped. ",
      "See ctrpv2_sample_slot_columns.csv."
    )
  }

  if (has_col(sample_dt, "ccle_primary_site")) {
    sample_dt[
      normalize_category_value(ccle_primary_site) == "soft_tissue",
      soft_tissue_match := TRUE
    ]
  } else {
    warning(
      "Column 'ccle_primary_site' was not found in CTRPv2 sample slot. ",
      "Soft Tissue matching from ccle_primary_site was skipped. ",
      "See ctrpv2_sample_slot_columns.csv."
    )
  }

  selected <- sample_dt[target_match == TRUE | soft_tissue_match == TRUE]

  if (nrow(selected) == 0) {
    stop(
      "No matching CTRPv2 samples found for target list, ",
      "tissueid == 'Soft Tissue', or ccle_primary_site == 'soft_tissue'. ",
      "See: ",
      file.path(OUT_DIR, "ctrpv2_sample_slot_columns.csv")
    )
  }

  if (!has_col(selected, "sampleid")) {
    stop(
      "CTRPv2 sample slot must contain 'sampleid' because ",
      "treatmentResponse$info$sampleid maps back to CTR@sample$sampleid."
    )
  }

  selected[, canonical_sample_source_id := clean_na(sampleid)]

  selected[, canonical_cell_line_name := coalesce_dt_cols(.SD, c(
    "sampleid",
    "sample_rowname",
    "ccl_name"
  ))]

  selected <- selected[
    !is.na(canonical_sample_source_id) &
      !is.na(canonical_cell_line_name)
  ]

  selected[, canonical_sample_id := make_prefixed_sample_id(canonical_sample_source_id)]
  selected[, normalized_sampleid_key := normalize_cell_line_name(canonical_sample_source_id)]

  selected
}

build_canonical_lookup_ctrpv2 <- function(selected_sample_dt) {
  alt_cols <- c(
    # Most important for treatmentResponse$info$sampleid mapping
    "sampleid",

    # Output identifiers
    "canonical_sample_source_id",
    "canonical_cell_line_name",
    "canonical_sample_id",

    # Other aliases for auditing / future use
    "sample_rowname",
    "ccl_name",
    "master_ccl_id",
    "PharmacoDB.id",
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

        # Final prefixed sample ID used in pre_clinical_sample.csv
        sample_id = selected_sample_dt$canonical_sample_id,

        # Original CTR@sample$sampleid
        source_sampleid = selected_sample_dt$canonical_sample_source_id,

        # Final cell-line name used in pre_clinical_cell_line.csv
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
# CTRPv2 treatment response extraction using treatmentResponse$info
# -------------------------------------------------------------------------

extract_treatment_response_ctrpv2 <- function(pset_ctrpv2, canonical_lookup, out_path) {
  profiles_obj <- pset_ctrpv2@treatmentResponse$profiles
  info_obj <- pset_ctrpv2@treatmentResponse$info

  if (is.null(profiles_obj)) {
    stop("CTRPv2 treatmentResponse$profiles is NULL.")
  }

  if (is.null(info_obj)) {
    stop("CTRPv2 treatmentResponse$info is NULL.")
  }

  profiles_dt <- as.data.table(
    profiles_obj,
    keep.rownames = "profile_id"
  )

  info_dt <- as.data.table(
    info_obj,
    keep.rownames = "info_rowname"
  )

  fwrite(
    data.table(column_name = colnames(profiles_dt)),
    file.path(OUT_DIR, "ctrpv2_treatment_response_profile_columns.csv")
  )

  fwrite(
    data.table(column_name = colnames(info_dt)),
    file.path(OUT_DIR, "ctrpv2_treatment_response_info_columns.csv")
  )

  cat("CTRPv2 treatment profiles columns:\n")
  print(colnames(profiles_dt))

  cat("CTRPv2 treatment info columns:\n")
  print(colnames(info_dt))

  if (!("profile_id" %in% colnames(profiles_dt))) {
    stop("Could not preserve rownames from CTRPv2 treatmentResponse$profiles.")
  }

  profile_ids <- as.character(profiles_dt$profile_id)

  info_key_col <- find_matching_key_col(
    dt = info_dt,
    profile_ids = profile_ids,
    candidates = c(
      "profile_id",
      "profileid",
      "profileID",
      "response_id",
      "responseid",
      "experiment_id",
      "experimentid",
      "rownames",
      "rowname",
      "info_rowname"
    )
  )

  if (is.na(info_key_col)) {
    stop(
      "Could not find a column/rownames in treatmentResponse$info that maps ",
      "back to treatmentResponse$profiles rownames. See: ",
      file.path(OUT_DIR, "ctrpv2_treatment_response_info_columns.csv")
    )
  }

  cat("Using treatmentResponse$info key column: ", info_key_col, "\n", sep = "")

  info_dt[, profile_id := as.character(get(info_key_col))]

  info_dt <- unique(info_dt, by = "profile_id")

  # Prefix info columns so we know they came from treatmentResponse$info.
  info_non_key_cols <- setdiff(colnames(info_dt), "profile_id")

  setnames(
    info_dt,
    info_non_key_cols,
    paste0("info__", info_non_key_cols)
  )

  tr_dt <- merge(
    profiles_dt,
    info_dt,
    by = "profile_id",
    all.x = TRUE
  )

  fwrite(
    head(tr_dt, 100),
    sub("\\.csv$", "_profiles_info_joined_head.csv", out_path)
  )

  sample_info_col <- first_existing_col(
    tr_dt,
    paste0(
      "info__",
      c(
        "sampleid",
        "sample_id",
        "sample.id",
        "SampleID",
        "Sample.ID",
        "cellid",
        "cell_id",
        "ccl_name"
      )
    )
  )

  treatment_info_col <- first_existing_col(
    tr_dt,
    paste0(
      "info__",
      c(
        "treatmentid",
        "treatment_id",
        "treatment.id",
        "TreatmentID",
        "Treatment.ID",
        "drugid",
        "drug_id",
        "compoundid",
        "compound_id",
        "master_cpd_id",
        "cpd_name",
        "drug_name",
        "treatment_name"
      )
    )
  )

  if (is.na(sample_info_col)) {
    stop(
      "Could not find sample ID column in treatmentResponse$info after joining. ",
      "Expected something like info__sampleid. See: ",
      file.path(OUT_DIR, "ctrpv2_treatment_response_info_columns.csv")
    )
  }

  if (is.na(treatment_info_col)) {
    stop(
      "Could not find treatment ID column in treatmentResponse$info after joining. ",
      "Expected something like info__treatmentid. See: ",
      file.path(OUT_DIR, "ctrpv2_treatment_response_info_columns.csv")
    )
  }

  cat("Using treatmentResponse$info sample column: ", sample_info_col, "\n", sep = "")
  cat("Using treatmentResponse$info treatment column: ", treatment_info_col, "\n", sep = "")

  tr_dt[, info_sampleid := clean_na(get(sample_info_col))]
  tr_dt[, treatment_id := clean_na(get(treatment_info_col))]
  tr_dt[, normalized_info_sampleid := normalize_cell_line_name(info_sampleid)]

  # Important: info$sampleid maps to CTR@sample$sampleid.
  sampleid_lookup <- canonical_lookup[
    source_column %in% c("sampleid", "canonical_sample_source_id"),
    .(
      normalized_lookup_name,
      source_sampleid,
      sample_id,
      cell_line_name
    )
  ]

  sampleid_lookup <- unique(sampleid_lookup, by = "normalized_lookup_name")

  if (nrow(sampleid_lookup) == 0) {
    stop(
      "No sampleid lookup rows found in canonical_lookup. ",
      "Cannot map treatmentResponse$info$sampleid back to CTR@sample$sampleid."
    )
  }

  unmatched_info_rows <- tr_dt[
    is.na(info_sampleid) |
      is.na(treatment_id) |
      !(normalized_info_sampleid %in% sampleid_lookup$normalized_lookup_name),
    .(
      profile_id,
      info_sampleid,
      treatment_id,
      normalized_info_sampleid
    )
  ]

  fwrite(
    unmatched_info_rows,
    sub("\\.csv$", "_unmatched_info_rows.csv", out_path)
  )

  cat("CTRPv2 treatment profile rows:", nrow(tr_dt), "\n")
  cat("Unmatched treatment info rows:", nrow(unmatched_info_rows), "\n")

  tr_mapped <- merge(
    tr_dt,
    sampleid_lookup,
    by.x = "normalized_info_sampleid",
    by.y = "normalized_lookup_name",
    all = FALSE
  )

  ic50_col <- first_existing_col(
    tr_mapped,
    c(
      "ic50_recomputed",
      "IC50",
      "ic50",
      "published_IC50",
      "published_ic50",
      "ic50_published"
    )
  )

  aac_col <- first_existing_col(
    tr_mapped,
    c(
      "aac_recomputed",
      "acc_recomputed",
      "AUC",
      "auc",
      "auc_recomputed",
      "published_ActArea",
      "published_AUC",
      "aac_published"
    )
  )

  moa_col <- first_existing_col(
    tr_mapped,
    paste0(
      "info__",
      c(
        "mechanism_of_action",
        "moa",
        "MOA",
        "target",
        "drug_target"
      )
    )
  )

  treatment_response_raw_mapped <- data.table(
    profile_id = tr_mapped$profile_id,

    # From treatmentResponse$info
    info_sampleid = tr_mapped$info_sampleid,
    treatment_id = tr_mapped$treatment_id,

    # Mapped back to CTR@sample$sampleid
    source_sampleid = tr_mapped$source_sampleid,

    # Final prefixed sample ID used in pre_clinical_sample.csv
    prefixed_sample_id = tr_mapped$sample_id,

    cell_line_name = tr_mapped$cell_line_name,

    ic50_recomputed = if (!is.na(ic50_col)) {
      suppressWarnings(as.numeric(tr_mapped[[ic50_col]]))
    } else {
      NA_real_
    },

    acc_recomputed = if (!is.na(aac_col)) {
      suppressWarnings(as.numeric(tr_mapped[[aac_col]]))
    } else {
      NA_real_
    },

    mechanism_of_action = if (!is.na(moa_col)) {
      clean_na(tr_mapped[[moa_col]])
    } else {
      NA_character_
    }
  )

  fwrite(
    treatment_response_raw_mapped,
    sub("\\.csv$", "_mapped_raw_rows.csv", out_path)
  )

  treatment_response_dt <- data.table(
    cell_line_name = clean_na(tr_mapped[["cell_line_name"]]),
    treatment_id = clean_na(tr_mapped[["treatment_id"]]),

    ic50_recomputed = if (!is.na(ic50_col)) {
      suppressWarnings(as.numeric(tr_mapped[[ic50_col]]))
    } else {
      NA_real_
    },

    # Keeping same output field as previous extractors / DB schema.
    # Source column is usually aac_recomputed.
    acc_recomputed = if (!is.na(aac_col)) {
      suppressWarnings(as.numeric(tr_mapped[[aac_col]]))
    } else {
      NA_real_
    },

    mechanism_of_action = if (!is.na(moa_col)) {
      clean_na(tr_mapped[[moa_col]])
    } else {
      NA_character_
    }
  )

  treatment_response_dt <- treatment_response_dt[
    !is.na(cell_line_name) &
      !is.na(treatment_id)
  ]

  duplicate_pairs <- treatment_response_dt[
    ,
    .N,
    by = .(cell_line_name, treatment_id)
  ][N > 1]

  fwrite(
    duplicate_pairs,
    sub("\\.csv$", "_duplicate_cell_line_treatment_pairs.csv", out_path)
  )

  # Prefer rows with the most non-NA numeric response values if duplicates exist.
  treatment_response_dt[, completeness_score :=
                          as.integer(!is.na(ic50_recomputed)) +
                          as.integer(!is.na(acc_recomputed))]

  setorder(
    treatment_response_dt,
    cell_line_name,
    treatment_id,
    -completeness_score
  )

  treatment_response_dt <- treatment_response_dt[
    ,
    .SD[1],
    by = .(cell_line_name, treatment_id)
  ]

  treatment_response_dt[, completeness_score := NULL]

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

  cat("Wrote CTRPv2 treatment response rows:", nrow(treatment_response_dt), "\n")
}

# -------------------------------------------------------------------------
# Load CTRPv2
# -------------------------------------------------------------------------

cat("\n==============================\n")
cat("Loading CTRPv2\n")
cat("==============================\n")

ctrpv2 <- read_updated_rds(CTRPV2_RDS_PATH)

# -------------------------------------------------------------------------
# Select samples/cell lines
# -------------------------------------------------------------------------

selected_sample_dt <- build_selected_sample_table_ctrpv2(
  pset_ctrpv2 = ctrpv2,
  target_cell_lines = TARGET_CELL_LINES
)

canonical_lookup <- build_canonical_lookup_ctrpv2(selected_sample_dt)

cat("Requested unique target cell lines:", length(unique(TARGET_CELL_LINES)), "\n")
cat("Matched selected CTRPv2 sample rows:", nrow(selected_sample_dt), "\n")
cat("Explicit target matches:", selected_sample_dt[target_match == TRUE, .N], "\n")
cat("Soft Tissue matches:", selected_sample_dt[soft_tissue_match == TRUE, .N], "\n")
cat("Canonical lookup rows:", nrow(canonical_lookup), "\n")

fwrite(
  canonical_lookup,
  file.path(OUT_DIR, "ctrpv2_canonical_lookup.csv")
)

fwrite(
  selected_sample_dt,
  file.path(OUT_DIR, "ctrpv2_selected_sample_slot_rows.csv")
)

# -------------------------------------------------------------------------
# pre_clinical_cell_line.csv
# Same fields as previous extractors:
#   cell_line_name, accession, category, sex, age
# -------------------------------------------------------------------------

cell_line_dt <- unique(
  data.table(
    cell_line_name = clean_na(selected_sample_dt[["canonical_cell_line_name"]]),

    accession = clean_na(safe_col(
      selected_sample_dt,
      c(
        "Cellosaurus.Accession.id",
        "Cellosaurus.Accession.ID",
        "cellosaurus.accession",
        "cellosaurus.accession.id",
        "accession",
        "cell_line_accession"
      )
    )),

    category = clean_na(safe_col(
      selected_sample_dt,
      c(
        "CellLine.Type",
        "cellosaurus.category",
        "Cellosaurus.Category",
        "category"
      )
    )),

    sex = clean_na(safe_col(
      selected_sample_dt,
      c(
        "cellosaurus.sexOfCell",
        "Cellosaurus.Sex",
        "sex",
        "Sex",
        "Gender",
        "gender"
      )
    )),

    age = parse_age_int(safe_col(
      selected_sample_dt,
      c(
        "cellosaurus.ageAtSampling",
        "Cellosaurus.Age",
        "age",
        "Age"
      )
    ))
  ),
  by = "cell_line_name"
)

cell_line_dt <- cell_line_dt[!is.na(cell_line_name)]

fwrite(
  cell_line_dt,
  file.path(OUT_DIR, "pre_clinical_cell_line.csv")
)

cat("Wrote cell lines:", nrow(cell_line_dt), "\n")

# -------------------------------------------------------------------------
# pre_clinical_sample.csv
# Same fields as previous extractors:
#   id, cell_line_name, site_primary, site_subtype1, site_subtype2,
#   site_subtype3, histology, histology_subtype1, histology_subtype2,
#   histology_subtype3, gender, age, race, diseases, disease_type
# -------------------------------------------------------------------------

sample_out_dt <- data.table(
  id = clean_na(selected_sample_dt[["canonical_sample_id"]]),
  cell_line_name = clean_na(selected_sample_dt[["canonical_cell_line_name"]]),

  site_primary = clean_na(safe_col(
    selected_sample_dt,
    c(
      "ccle_primary_site",
      "site_primary",
      "Site.Primary",
      "primary_site",
      "tissueid"
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
      "tissueid"
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
# pre_clinical_treatment_response.csv
# Same fields as previous extractors:
#   cell_line_name, treatment_id, ic50_recomputed, acc_recomputed,
#   mechanism_of_action
# -------------------------------------------------------------------------

extract_treatment_response_ctrpv2(
  pset_ctrpv2 = ctrpv2,
  canonical_lookup = canonical_lookup,
  out_path = file.path(OUT_DIR, "pre_clinical_treatment_response.csv")
)

# -------------------------------------------------------------------------
# Cleanup
# -------------------------------------------------------------------------

cat("\nFreeing CTRPv2 objects from memory\n")

rm(
  ctrpv2,
  selected_sample_dt,
  canonical_lookup,
  cell_line_dt,
  sample_out_dt
)

gc(verbose = TRUE)

cat("Finished extracting selected CTRPv2 preclinical CSVs into:", OUT_DIR, "\n")