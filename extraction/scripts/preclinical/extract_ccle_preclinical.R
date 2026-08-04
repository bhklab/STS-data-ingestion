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

  # Keep basename in case filenames contain paths.
  x <- basename(x)

  # Remove common Affymetrix file extensions.
  x <- sub("\\.CEL(\\.GZ)?$", "", x, ignore.case = TRUE)

  # Normalize case because old sample metadata may use NIECE_p_...
  # while matrix colnames use NIECE_P_...
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

# -------------------------------------------------------------------------
# CCLE 2019 selection helpers
# -------------------------------------------------------------------------

build_selected_sample_table_2019 <- function(pset_2019, target_cell_lines) {
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

  selected <- sample_dt[target_match == TRUE | soft_tissue_match == TRUE]

  if (nrow(selected) == 0) {
    stop(
      "No matching CCLE 2019 samples found for target list or soft_tissue filter."
    )
  }

  selected
}

build_canonical_lookup <- function(selected_sample_dt) {
  canonical_sample_id <- clean_na(selected_sample_dt[["sampleid"]])
  canonical_cell_line <- clean_na(selected_sample_dt[["cellosaurus.cellLineName"]])

  alt_cols <- c(
    "cellosaurus.cellLineName",
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
        cell_line_name = canonical_cell_line
      )
    }
  }

  lookup_dt <- rbindlist(pieces, fill = TRUE)

  lookup_dt <- lookup_dt[
    !is.na(normalized_lookup_name) &
      !is.na(sample_id) &
      !is.na(cell_line_name)
  ]

  unique(lookup_dt, by = "normalized_lookup_name")
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
    # CCLE 2019-style columns
    "cellosaurus.cellLineName",
    "CCLE.name",
    "CCLE.sampleid",
    "sampleid",
    "unique.sampleid",
    "dataset_sample_id",
    "ccle_sample_id",

    # Some CCLE 2015-style columns
    "Cell.line.primary.name",
    "Cell_Line",
    "Name",
    "samplename",
    "filename",
    "Expression.arrays",
    "SNP.arrays",

    # fallbacks
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

  mapped <- merge(
    candidate_map,
    canonical_lookup[, .(normalized_lookup_name, sample_id, cell_line_name)],
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

  # Direct mapping is preferred because the CCLE_2015 rna colData contains
  # sampleid, CCLE.name, and Cell.line.primary.name for each expression array.
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

  if (length(rna_direct_pieces) > 0) {
    rna_direct_candidates <- rbindlist(rna_direct_pieces, fill = TRUE)
    rna_direct_candidates <- rna_direct_candidates[!is.na(normalized_lookup_name)]

    direct_mapped <- merge(
      rna_direct_candidates,
      canonical_lookup[, .(normalized_lookup_name, sample_id, cell_line_name)],
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

  # Array-ID fallback. This maps:
  # rna colname/filename/Expression.arrays -> CCLE_2015@sample Expression.arrays
  # -> CCLE_2015@sample Cell.line.primary.name -> CCLE_2019 canonical lookup.
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

      via_sample[, normalized_lookup_name := normalize_cell_line_name(ccle_2015_cell_line_name)]

      via_sample_mapped <- merge(
        via_sample,
        canonical_lookup[, .(normalized_lookup_name, sample_id, cell_line_name)],
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

  # direct_mapped is bound first, and direct cols are ordered with
  # Cell.line.primary.name first, so this keeps the best direct match first.
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
# CCLE 2015 treatment response
# -------------------------------------------------------------------------

extract_treatment_response_2015 <- function(pset_2015, canonical_lookup, out_path) {
  profiles_obj <- pset_2015@treatmentResponse$profiles

  # Older CCLE_2015 stores treatment/sample identifiers in rownames:
  #   drugid_AEW541_A-204
  #   drugid_Nilotinib_TE 617.T
  tr_dt <- as.data.table(
    profiles_obj,
    keep.rownames = "profile_id"
  )

  cat("CCLE 2015 treatment profile columns:\n")
  print(colnames(tr_dt))

  if (!("profile_id" %in% colnames(tr_dt))) {
    stop("Could not preserve rownames from CCLE 2015 treatmentResponse$profiles.")
  }

  fwrite(
    head(tr_dt, 100),
    sub("\\.csv$", "_raw_profile_head.csv", out_path)
  )

  parse_one_profile_id <- function(profile_id, valid_norm_names) {
    profile_id <- as.character(profile_id)

    if (is.na(profile_id) || profile_id == "") {
      return(list(treatment_id = NA_character_, sampleid = NA_character_))
    }

    without_prefix <- sub("^drugid_", "", profile_id)

    parts <- strsplit(without_prefix, "_", fixed = TRUE)[[1]]

    if (length(parts) < 2) {
      return(list(treatment_id = NA_character_, sampleid = NA_character_))
    }

    best_treatment <- NA_character_
    best_sample <- NA_character_

    for (split_idx in seq_len(length(parts) - 1)) {
      candidate_treatment <- paste(parts[seq_len(split_idx)], collapse = "_")
      candidate_sample <- paste(parts[(split_idx + 1):length(parts)], collapse = "_")

      if (normalize_cell_line_name(candidate_sample) %in% valid_norm_names) {
        best_treatment <- candidate_treatment
        best_sample <- candidate_sample
        break
      }
    }

    list(
      treatment_id = best_treatment,
      sampleid = best_sample
    )
  }

  valid_norm_names <- unique(canonical_lookup$normalized_lookup_name)

  parsed <- rbindlist(
    lapply(
      tr_dt$profile_id,
      function(profile_id) {
        parsed_one <- parse_one_profile_id(profile_id, valid_norm_names)

        data.table(
          treatment_id = parsed_one$treatment_id,
          sampleid = parsed_one$sampleid
        )
      }
    )
  )

  tr_dt[, treatment_id := parsed$treatment_id]
  tr_dt[, sampleid := parsed$sampleid]
  tr_dt[, normalized_sample := normalize_cell_line_name(sampleid)]

  unmatched_profile_ids <- tr_dt[
    is.na(treatment_id) |
      is.na(sampleid) |
      !(normalized_sample %in% valid_norm_names),
    .(profile_id, treatment_id, sampleid, normalized_sample)
  ]

  fwrite(
    unmatched_profile_ids,
    sub("\\.csv$", "_unmatched_profile_ids.csv", out_path)
  )

  cat("CCLE 2015 treatment profile rows:", nrow(tr_dt), "\n")
  cat("Unmatched treatment profile IDs:", nrow(unmatched_profile_ids), "\n")

  tr_mapped <- merge(
    tr_dt,
    canonical_lookup[, .(normalized_lookup_name, cell_line_name)],
    by.x = "normalized_sample",
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
      "AUC",
      "auc",
      "auc_recomputed",
      "published_ActArea",
      "published_AUC",
      "aac_published"
    )
  )

  treatment_response_dt <- data.table(
    cell_line_name = clean_na(tr_mapped[["cell_line_name"]]),
    treatment_id = clean_na(tr_mapped[["treatment_id"]]),
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
    mechanism_of_action = NA_character_
  )

  treatment_response_dt <- treatment_response_dt[
    !is.na(cell_line_name) &
      !is.na(treatment_id)
  ]

  treatment_response_dt <- unique(
    treatment_response_dt,
    by = c("cell_line_name", "treatment_id")
  )

  fwrite(treatment_response_dt, out_path)

  cat("Wrote CCLE 2015 treatment response rows:", nrow(treatment_response_dt), "\n")
}

# -------------------------------------------------------------------------
# Phase 1: Load CCLE 2019 only
# -------------------------------------------------------------------------

cat("\n==============================\n")
cat("Phase 1: Loading CCLE 2019\n")
cat("==============================\n")

ccle_2019 <- read_updated_rds(CCLE_2019_RDS_PATH)

selected_sample_dt <- build_selected_sample_table_2019(
  pset_2019 = ccle_2019,
  target_cell_lines = TARGET_CELL_LINES
)

canonical_lookup <- build_canonical_lookup(selected_sample_dt)

cat("Requested unique target cell lines:", length(unique(TARGET_CELL_LINES)), "\n")
cat("Matched selected CCLE 2019 sample rows:", nrow(selected_sample_dt), "\n")
cat("Explicit target matches:", selected_sample_dt[target_match == TRUE, .N], "\n")
cat("Soft-tissue matches:", selected_sample_dt[soft_tissue_match == TRUE, .N], "\n")
cat("Canonical lookup rows:", nrow(canonical_lookup), "\n")

canonical_lookup_path <- file.path(OUT_DIR, "ccle_2019_canonical_lookup.csv")

fwrite(
  canonical_lookup,
  canonical_lookup_path
)

cat("Wrote canonical lookup:", canonical_lookup_path, "\n")

# -------------------------------------------------------------------------
# pre_clinical_cell_line.csv from CCLE 2019
# -------------------------------------------------------------------------

cell_line_dt <- unique(
  data.table(
    cell_line_name = clean_na(selected_sample_dt[["cellosaurus.cellLineName"]]),
    accession = clean_na(selected_sample_dt[["cellosaurus.accession"]]),
    category = clean_na(selected_sample_dt[["cellosaurus.category"]]),
    sex = clean_na(selected_sample_dt[["cellosaurus.sexOfCell"]]),
    age = parse_age_int(selected_sample_dt[["cellosaurus.ageAtSampling"]])
  ),
  by = "cell_line_name"
)

cell_line_dt <- cell_line_dt[!is.na(cell_line_name)]

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
  cell_line_name = clean_na(selected_sample_dt[["cellosaurus.cellLineName"]]),

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
# CCLE 2015 treatment response only
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