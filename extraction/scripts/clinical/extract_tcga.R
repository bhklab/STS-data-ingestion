#!/usr/bin/env Rscript
#
# Extract clinical-table CSVs from a TCGA-SARC MultiAssayExperiment (.qs) for
# loading into the STS-data-ingestion clinical schema defined in
# seeding/models/tables.py.
#
# Style/conventions follow the existing preclinical extraction pipeline
# (extraction/scripts/preclinical/extract_ctrpv2_preclinical.R) and the
# preclinical seeder (seeding/scripts/preclinical_seed_ctrpv2.py): explicit
# config constants, small `clean_*` helpers, one function per output CSV,
# and audit-friendly `cat()` progress messages.
#
# Run from the repository root:
#   R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/clinical/extract_tcga.R
#
# Required raw input (see README section "1. Required raw input files" for
# the equivalent preclinical convention):
#   extraction/data/raw/clinical/TCGA_SARC_mae_final.qs
#
# -------------------------------------------------------------------------
# What this produces (all under OUT_DIR)
# -------------------------------------------------------------------------
#
# clinical_sample.csv
#   id,race,ethnicity,sex,age,histology,tissue,tissue_origin
#   -> seeds clinical_sample. dataset_id is added during seeding, as with
#      pre_clinical_cell_line.csv in the preclinical pipeline.
#   One row per distinct sample barcode (colname) across all six assays,
#   e.g. TCGA-DX-A6B7-01A. A handful of patients contribute more than one
#   sample within a given assay (extra portions/relapse/normal samples);
#   demographic fields are copied from that sample's patient-level colData
#   row via the MAE sampleMap.
#
# pre_clinical_gene.csv
#   id,name
#   -> seeds the shared pre_clinical_gene table (Ensembl gene ID with the
#      version suffix stripped, e.g. ENSG00000002586.20_PAR_Y ->
#      ENSG00000002586). Built from the union of RNA and CNV rowData, which
#      is also used to resolve Hugo symbols (Mutation rownames) to Ensembl
#      IDs. Existing gene rows from the preclinical load are left as-is;
#      the seeder inserts only genes that are not already present.
#
# clinical_rna.csv           sample_id,gene_id,value          (tpm_unstrand)
# clinical_cnv.csv           sample_id,gene_id,value          (copy_number)
# clinical_mutation.csv      sample_id,gene_id,mutation,oncoprint
# clinical_mirna.csv         id,sample_id,gene_id,value       (rpm; gene_id
#                                                               always NA --
#                                                               no reliable
#                                                               mature-miRNA
#                                                               -> Ensembl
#                                                               mapping)
# clinical_antigen.csv       id,catalogue_number,peptide_target,peptide_target_gene
# clinical_rppa.csv          sample_id,antigen_id,value       (expression)
# clinical_probe.csv         id,name
# clinical_methylation.csv   sample_id,probe_id,value         (beta values)
#
# -------------------------------------------------------------------------
# Notes on scale
# -------------------------------------------------------------------------
#
# RNA (60,660 genes x ~265 samples) and Methylation (485,577 probes x ~269
# samples) are large in long format (tens of millions of rows each). Both
# are written in column chunks with fwrite(..., append = TRUE) to keep
# memory bounded. Use --limit-genes / --limit-samples for a fast smoke test
# before a full run, and the EXTRACT_* toggles below to skip assays you
# don't need yet.

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

# -------------------------------------------------------------------------
# Config
# -------------------------------------------------------------------------

MAE_PATH <- "extraction/data/raw/clinical/TCGA_SARC_mae_final.qs"
OUT_DIR <- "extraction/data/proc/clinical/TCGA_SARC"

CHUNK_SAMPLES <- 25 # samples per chunk when melting wide assay matrices

# -------------------------------------------------------------------------
# CLI args (only input and output file locations)
# -------------------------------------------------------------------------

parse_args <- function(argv) {
  out <- list(out_dir = NULL, mae_path = NULL)
  for (arg in argv) {
    if (startsWith(arg, "--out-dir=")) {
      out$out_dir <- sub("--out-dir=", "", arg)
    } else if (startsWith(arg, "--mae-path=")) {
      out$mae_path <- sub("--mae-path=", "", arg)
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!is.null(args$out_dir)) {
  OUT_DIR <- args$out_dir
}
if (!is.null(args$mae_path)) {
  MAE_PATH <- args$mae_path
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# General helpers (mirrors clean_* helpers in preclinical_seed_ctrpv2.py)
# -------------------------------------------------------------------------

clean_na <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | x == "NS" | x == "N/A" | x == "Not Reported" | is.na(x)] <- NA_character_
  x
}

clean_gene_id <- function(x) {
  x <- clean_na(x)
  ifelse(is.na(x), NA_character_, sub("\\..*$", "", x))
}

clean_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x) | !is.finite(x), NA_integer_, as.integer(round(x)))
}

clean_float <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x) | !is.finite(x), NA_real_, x)
}

clean_bool_from_01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA, x != 0)
}

first_token <- function(x, sep = ";") {
  x <- clean_na(x)
  ifelse(is.na(x), NA_character_, vapply(strsplit(x, sep, fixed = TRUE), `[`, character(1), 1))
}

# Writes a data.table to `path`, truncating on first write and appending on
# subsequent calls, so large assays can be streamed in column chunks.
write_chunk <- function(dt, path, first_write) {
  fwrite(dt, path, append = !first_write, col.names = first_write)
  invisible(NULL)
}

# -------------------------------------------------------------------------
# Load the MAE
# -------------------------------------------------------------------------

if (!file.exists(MAE_PATH)) {
  stop(
    "Missing required raw MAE file: ", MAE_PATH, ". ",
    "Place TCGA_SARC_mae_final.qs there (or pass --mae-path=...)."
  )
}

cat("Loading MAE from", MAE_PATH, "\n")
mae <- qread(MAE_PATH)

col_data <- as.data.frame(colData(mae))
sample_map_dt <- as.data.table(as.data.frame(sampleMap(mae)))
setnames(sample_map_dt, c("assay", "primary", "colname"))

cat("Patients in colData:", nrow(col_data), "\n")
cat("Distinct sample barcodes across all assays:", length(unique(sample_map_dt$colname)), "\n")

# -------------------------------------------------------------------------
# clinical_sample.csv
# -------------------------------------------------------------------------

extract_clinical_sample <- function(col_data, sample_map_dt, out_dir) {
  sample_dt <- unique(sample_map_dt[, .(id = colname, patientid = primary)])

  col_dt <- as.data.table(col_data, keep.rownames = "patientid")
  demo_cols <- c(
    "patientid", "race_standardized", "ethnicity_standardized", "sex",
    "age", "histo", "cancer_type", "tissue_or_organ_of_origin"
  )
  missing_demo <- setdiff(demo_cols, colnames(col_dt))
  if (length(missing_demo) > 0) {
    stop(
      "colData(mae) is missing expected columns: ",
      paste(missing_demo, collapse = ", ")
    )
  }
  col_dt <- col_dt[, ..demo_cols]

  out <- merge(sample_dt, col_dt, by = "patientid", all.x = TRUE)

  out <- data.table(
    id = clean_na(out$id),
    race = clean_na(out$race_standardized),
    ethnicity = clean_na(out$ethnicity_standardized),
    sex = clean_na(out$sex),
    age = clean_int(out$age),
    histology = clean_na(out$histo),
    tissue = clean_na(out$cancer_type),
    tissue_origin = clean_na(out$tissue_or_organ_of_origin)
  )

  out <- out[!is.na(id)]
  out <- unique(out, by = "id")
  setorder(out, id)

  fwrite(out, file.path(out_dir, "clinical_sample.csv"))
  cat("Wrote clinical_sample.csv:", nrow(out), "rows\n")
  invisible(out)
}

# -------------------------------------------------------------------------
# pre_clinical_gene.csv + Hugo-symbol -> Ensembl-ID lookup
# -------------------------------------------------------------------------

build_gene_table_and_symbol_lookup <- function(mae, out_dir) {
  rna_row <- as.data.table(rowData(experiments(mae)[["RNA"]]))
  cnv_row <- as.data.table(rowData(experiments(mae)[["CNV"]]))

  gene_dt <- rbindlist(
    list(
      rna_row[, .(id = clean_gene_id(gene_id), name = clean_na(gene_name))],
      cnv_row[, .(id = clean_gene_id(gene_id), name = clean_na(gene_name))]
    ),
    fill = TRUE
  )
  gene_dt <- gene_dt[!is.na(id)]
  gene_dt <- unique(gene_dt, by = "id")
  setorder(gene_dt, id)

  fwrite(gene_dt, file.path(out_dir, "pre_clinical_gene.csv"))
  cat("Wrote pre_clinical_gene.csv:", nrow(gene_dt), "rows\n")

  # Hugo symbol -> Ensembl ID, used to resolve Mutation rownames (Hugo_Symbol)
  # to gene_id. Kept separate from gene_dt because gene symbols are not
  # guaranteed unique; first match wins.
  symbol_lookup <- unique(gene_dt[!is.na(name), .(name, id)], by = "name")
  setnames(symbol_lookup, c("symbol", "gene_id"))

  list(gene_dt = gene_dt, symbol_lookup = symbol_lookup)
}

# -------------------------------------------------------------------------
# Generic: wide gene x sample matrix -> long CSV, chunked over samples
# -------------------------------------------------------------------------

# Generic: wide gene x sample matrix -> long CSV, chunked over samples
# -------------------------------------------------------------------------

extract_gene_matrix_long <- function(
  mat,
  out_path,
  value_name,
  value_clean_fn,
  chunk_samples
) {
  n_genes <- nrow(mat)
  n_samples <- ncol(mat)

  total_rows <- 0
  first_write <- TRUE
  sample_cols <- colnames(mat)[seq_len(n_samples)]

  for (start in seq(1, n_samples, by = chunk_samples)) {
    end <- min(start + chunk_samples - 1, n_samples)
    cols <- sample_cols[start:end]

    sub_mat <- mat[seq_len(n_genes), cols, drop = FALSE]

    # as.table()/as.data.table() melts a matrix into (rowname, colname, value)
    # triples -- the first two columns are the actual dimname labels, not
    # positional indices.
    long_dt <- as.data.table(as.table(sub_mat))
    setnames(long_dt, c("gene_id", "sample_id", "value"))
    long_dt[, gene_id := clean_gene_id(as.character(gene_id))]
    long_dt[, sample_id := clean_na(as.character(sample_id))]
    long_dt[, value := value_clean_fn(value)]
    long_dt <- long_dt[!is.na(gene_id) & !is.na(sample_id) & !is.na(value)]
    setnames(long_dt, "value", value_name)

    if (nrow(long_dt) > 0) {
      write_chunk(long_dt, out_path, first_write)
      first_write <- FALSE
      total_rows <- total_rows + nrow(long_dt)
    }
  }

  if (first_write) {
    # Nothing written (e.g. all-NA matrix); still emit a header-only file.
    fwrite(
      data.table(sample_id = character(), gene_id = character())[, (value_name) := numeric()],
      out_path
    )
  }

  cat("Wrote", out_path, ":", total_rows, "rows\n")
  invisible(total_rows)
}

# -------------------------------------------------------------------------
# clinical_rna.csv
# -------------------------------------------------------------------------

extract_clinical_rna <- function(mae, out_dir) {
  rna_se <- experiments(mae)[["RNA"]]
  mat <- assay(rna_se, "tpm_unstrand")

  extract_gene_matrix_long(
    mat = mat,
    out_path = file.path(out_dir, "clinical_rna.csv"),
    value_name = "value",
    value_clean_fn = clean_float,
    chunk_samples = CHUNK_SAMPLES
  )
}

# -------------------------------------------------------------------------
# clinical_cnv.csv
# -------------------------------------------------------------------------

extract_clinical_cnv <- function(mae, out_dir) {
  cnv_se <- experiments(mae)[["CNV"]]
  mat <- assay(cnv_se, "copy_number")

  extract_gene_matrix_long(
    mat = mat,
    out_path = file.path(out_dir, "clinical_cnv.csv"),
    value_name = "value",
    value_clean_fn = clean_int,
    chunk_samples = CHUNK_SAMPLES
  )
}

# -------------------------------------------------------------------------
# clinical_mutation.csv (mutation + oncoprint share dims; process together)
# -------------------------------------------------------------------------

extract_clinical_mutation <- function(mae, symbol_lookup, out_dir) {
  mut_se <- experiments(mae)[["Mutation"]]
  mutation_mat <- assay(mut_se, "mutation")
  oncoprint_mat <- assay(mut_se, "oncoprint")

  n_genes <- nrow(mutation_mat)
  n_samples <- ncol(mutation_mat)

  hugo_symbols <- clean_na(rownames(mutation_mat)[seq_len(n_genes)])
  n_unmapped <- sum(!is.na(hugo_symbols) & !(hugo_symbols %in% symbol_lookup$symbol))
  if (n_unmapped > 0) {
    cat(
      "Mutation: skipping", n_unmapped,
      "gene rows whose Hugo symbol did not resolve to an Ensembl ID (no pre_clinical_gene match)\n"
    )
  }

  out_path <- file.path(out_dir, "clinical_mutation.csv")
  sample_cols <- colnames(mutation_mat)[seq_len(n_samples)]
  total_rows <- 0
  first_write <- TRUE

  for (start in seq(1, n_samples, by = CHUNK_SAMPLES)) {
    end <- min(start + CHUNK_SAMPLES - 1, n_samples)
    cols <- sample_cols[start:end]

    sub_mut <- mutation_mat[seq_len(n_genes), cols, drop = FALSE]
    sub_onc <- oncoprint_mat[seq_len(n_genes), cols, drop = FALSE]

    mut_long <- as.data.table(as.table(sub_mut))
    setnames(mut_long, c("hugo_symbol", "sample_id", "mutation"))
    onc_long <- as.data.table(as.table(sub_onc))
    setnames(onc_long, c("hugo_symbol", "sample_id", "oncoprint"))

    long_dt <- cbind(mut_long, oncoprint = onc_long$oncoprint)
    long_dt[, hugo_symbol := clean_na(as.character(hugo_symbol))]
    long_dt[, sample_id := clean_na(as.character(sample_id))]
    long_dt[, mutation := clean_bool_from_01(mutation)]
    long_dt[, oncoprint := clean_na(as.character(oncoprint))]
    long_dt[, gene_id := symbol_lookup$gene_id[match(hugo_symbol, symbol_lookup$symbol)]]

    long_dt <- long_dt[!is.na(gene_id) & !is.na(sample_id) & (!is.na(mutation) | !is.na(oncoprint))]
    long_dt <- long_dt[, .(sample_id, gene_id, mutation, oncoprint)]

    if (nrow(long_dt) > 0) {
      write_chunk(long_dt, out_path, first_write)
      first_write <- FALSE
      total_rows <- total_rows + nrow(long_dt)
    }
  }

  if (first_write) {
    fwrite(
      data.table(sample_id = character(), gene_id = character(), mutation = logical(), oncoprint = character()),
      out_path
    )
  }

  cat("Wrote clinical_mutation.csv:", total_rows, "rows\n")
}

# -------------------------------------------------------------------------
# clinical_mirna.csv
# -------------------------------------------------------------------------

extract_clinical_mirna <- function(mae, out_dir) {
  mirna_se <- experiments(mae)[["miRNA"]]
  mat <- assay(mirna_se, "rpm")

  n_genes <- nrow(mat)
  n_samples <- ncol(mat)

  out_path <- file.path(out_dir, "clinical_mirna.csv")
  sample_cols <- colnames(mat)[seq_len(n_samples)]
  total_rows <- 0
  first_write <- TRUE

  for (start in seq(1, n_samples, by = CHUNK_SAMPLES)) {
    end <- min(start + CHUNK_SAMPLES - 1, n_samples)
    cols <- sample_cols[start:end]

    sub_mat <- mat[seq_len(n_genes), cols, drop = FALSE]
    long_dt <- as.data.table(as.table(sub_mat))
    setnames(long_dt, c("id", "sample_id", "value"))
    long_dt[, id := clean_na(as.character(id))]
    long_dt[, sample_id := clean_na(as.character(sample_id))]
    long_dt[, value := clean_float(value)]
    long_dt[, gene_id := NA_character_] # no reliable mature-miRNA -> Ensembl mapping
    long_dt <- long_dt[!is.na(id) & !is.na(sample_id) & !is.na(value)]
    long_dt <- long_dt[, .(id, sample_id, gene_id, value)]

    if (nrow(long_dt) > 0) {
      write_chunk(long_dt, out_path, first_write)
      first_write <- FALSE
      total_rows <- total_rows + nrow(long_dt)
    }
  }

  if (first_write) {
    fwrite(
      data.table(id = character(), sample_id = character(), gene_id = character(), value = numeric()),
      out_path
    )
  }

  cat("Wrote clinical_mirna.csv:", total_rows, "rows\n")
}

# -------------------------------------------------------------------------
# clinical_antigen.csv + clinical_rppa.csv
# -------------------------------------------------------------------------

extract_clinical_antigen_and_rppa <- function(mae, gene_dt, out_dir) {
  rppa_se <- experiments(mae)[["RPPA"]]
  row_dt <- as.data.table(rowData(rppa_se))

  antigen_dt <- data.table(
    id = clean_na(row_dt$AGID),
    catalogue_number = clean_na(row_dt$catalog_number),
    peptide_target = clean_na(row_dt$peptide_target)
  )
  # Best-effort exact match of the RPPA peptide-target label to a known gene
  # symbol; left NA when there is no exact match (antigen names commonly use
  # non-standard suffixes, e.g. phospho-site tags such as "_pS65").
  gene_by_symbol <- unique(gene_dt[!is.na(name), .(name = toupper(name), id)], by = "name")
  antigen_dt[, peptide_target_gene := gene_by_symbol$id[match(toupper(peptide_target), gene_by_symbol$name)]]
  antigen_dt <- antigen_dt[!is.na(id)]
  antigen_dt <- unique(antigen_dt, by = "id")
  setorder(antigen_dt, id)

  fwrite(antigen_dt, file.path(out_dir, "clinical_antigen.csv"))
  cat("Wrote clinical_antigen.csv:", nrow(antigen_dt), "rows\n")

  mat <- assay(rppa_se, "expression")

  n_genes <- nrow(mat)
  n_samples <- ncol(mat)

  out_path <- file.path(out_dir, "clinical_rppa.csv")
  sample_cols <- colnames(mat)[seq_len(n_samples)]
  total_rows <- 0
  first_write <- TRUE

  for (start in seq(1, n_samples, by = CHUNK_SAMPLES)) {
    end <- min(start + CHUNK_SAMPLES - 1, n_samples)
    cols <- sample_cols[start:end]

    sub_mat <- mat[seq_len(n_genes), cols, drop = FALSE]
    long_dt <- as.data.table(as.table(sub_mat))
    setnames(long_dt, c("antigen_id", "sample_id", "value"))
    long_dt[, antigen_id := clean_na(as.character(antigen_id))]
    long_dt[, sample_id := clean_na(as.character(sample_id))]
    long_dt[, value := clean_float(value)]
    long_dt <- long_dt[!is.na(antigen_id) & !is.na(sample_id) & !is.na(value)]
    long_dt <- long_dt[, .(sample_id, antigen_id, value)]

    if (nrow(long_dt) > 0) {
      write_chunk(long_dt, out_path, first_write)
      first_write <- FALSE
      total_rows <- total_rows + nrow(long_dt)
    }
  }

  if (first_write) {
    fwrite(data.table(sample_id = character(), antigen_id = character(), value = numeric()), out_path)
  }

  cat("Wrote clinical_rppa.csv:", total_rows, "rows\n")
}

# -------------------------------------------------------------------------
# clinical_probe.csv + clinical_methylation.csv
# -------------------------------------------------------------------------

extract_clinical_probe_and_methylation <- function(mae, out_dir) {
  methyl_se <- experiments(mae)[["Methylation"]]
  row_dt <- as.data.table(rowData(methyl_se))

  probe_dt <- data.table(
    id = clean_na(rownames(methyl_se)),
    name = first_token(row_dt$gene, sep = ";")
  )
  probe_dt <- probe_dt[!is.na(id)]
  probe_dt <- unique(probe_dt, by = "id")
  setorder(probe_dt, id)

  fwrite(probe_dt, file.path(out_dir, "clinical_probe.csv"))
  cat("Wrote clinical_probe.csv:", nrow(probe_dt), "rows\n")

  mat <- assay(methyl_se, 1) # unnamed single assay of beta values

  n_genes <- nrow(mat)
  n_samples <- ncol(mat)

  out_path <- file.path(out_dir, "clinical_methylation.csv")
  sample_cols <- colnames(mat)[seq_len(n_samples)]
  total_rows <- 0
  first_write <- TRUE

  for (start in seq(1, n_samples, by = CHUNK_SAMPLES)) {
    end <- min(start + CHUNK_SAMPLES - 1, n_samples)
    cols <- sample_cols[start:end]

    sub_mat <- mat[seq_len(n_genes), cols, drop = FALSE]
    long_dt <- as.data.table(as.table(sub_mat))
    setnames(long_dt, c("probe_id", "sample_id", "value"))
    long_dt[, probe_id := clean_na(as.character(probe_id))]
    long_dt[, sample_id := clean_na(as.character(sample_id))]
    long_dt[, value := clean_float(value)]
    long_dt <- long_dt[!is.na(probe_id) & !is.na(sample_id) & !is.na(value)]
    long_dt <- long_dt[, .(sample_id, probe_id, value)]

    if (nrow(long_dt) > 0) {
      write_chunk(long_dt, out_path, first_write)
      first_write <- FALSE
      total_rows <- total_rows + nrow(long_dt)
    }
  }

  if (first_write) {
    fwrite(data.table(sample_id = character(), probe_id = character(), value = numeric()), out_path)
  }

  cat("Wrote clinical_methylation.csv:", total_rows, "rows\n")
}

# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------

invisible(extract_clinical_sample(col_data, sample_map_dt, OUT_DIR))
gene_info <- build_gene_table_and_symbol_lookup(mae, OUT_DIR)

extract_clinical_rna(mae, OUT_DIR)
extract_clinical_cnv(mae, OUT_DIR)
extract_clinical_mutation(mae, gene_info$symbol_lookup, OUT_DIR)
extract_clinical_mirna(mae, OUT_DIR)
extract_clinical_antigen_and_rppa(mae, gene_info$gene_dt, OUT_DIR)
extract_clinical_probe_and_methylation(mae, OUT_DIR)

cat("\nFinished extracting TCGA_SARC clinical CSVs into:", OUT_DIR, "\n")
