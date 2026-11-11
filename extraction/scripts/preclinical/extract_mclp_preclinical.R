suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(PharmacoGx)
  library(miRBaseConverter)
  library(reshape2)
})

#MCLP RDS PATH
MCLP_PATH <- "extraction/data/raw/preclinical/MCLP_RPPA_SE.rds"

#QC Cell Line CSV PATH
STS_CELL_LINE_PATH <- "extraction/data/raw/preclinical/All_Psets_sarcoma_cell_line_QC.csv"

#Output Path
OUTPUT_PATH <- "extraction/data/proc/preclinical/MCLP"

dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)

#extract and setup cell line information from the MCLP pset object
filter_cell_lines <- function(pset) {

  qc_cell_lines <- read.csv(STS_CELL_LINE_PATH, stringsAsFactors = FALSE)
  qc_cell_lines <- qc_cell_lines[qc_cell_lines$mod_tissueid == "Soft Tissue", ]
  
  if ("sample" %in% slotNames(pset)) {
    sample_info <- pset@sample
  } else {
    sample_info <- pset@colData
  }

  cols_to_check <- c('cell_line', 'sampleid')
  cell_line_column <- NULL

  for (col in cols_to_check) {
    if (col %in% colnames(sample_info)) {
      cell_line_column <- col
      break
    }
  }

  if (is.null(cell_line_column)) {
    stop("No cell line column found in the pset object.")
  }

  cell_line_table <- qc_cell_lines[tolower(qc_cell_lines$cell_line) %in% tolower(sample_info[[cell_line_column]]), c('cell_line', 'tissueid', 'mod_tissueid', 'Cellosaurus.Accession.id', 'cellosaurus.category', 'cellosaurus.sex', 'cellosaurus.age')]

  colnames(cell_line_table) <- c('cell_line_name', 'tissueid', 'mod_tissueid', 'accession', 'category', 'sex', 'age')

  return(cell_line_table[!duplicated(cell_line_table$cell_line_name), ])

    
}
#Construct sample table from the filtered cell line dataframe
create_sample_table <- function(cell_line_dataframe) {

    sample_table <- data.frame(
        id = cell_line_dataframe$cell_line_name,
        cell_line_name = cell_line_dataframe$cell_line_name,
        stringsAsFactors = FALSE
    )
    sample_table$id <- paste0("MCLP_", sample_table$id)
    return(sample_table)

}

#Extract and construct RPPA data from pSet object
extract_rppa_data <- function(pset, samples) {
    mclp_assay <- assay(pset)
    row_data <- rowData(pset)
    mclp_rppa <- mclp_assay[, colnames(mclp_assay) %in% samples]
    mclp_rppa <- melt(mclp_rppa, varnames = c("feature_id", "cell_line"), value.name = "value")
    mclp_rppa <- mclp_rppa[!is.na(mclp_rppa$value), ]
    mclp_rppa <- merge(mclp_rppa, row_data, by.x = "feature_id", by.y = "row.names", all.x = TRUE)
    mclp_rppa <- mclp_rppa[c("feature_id", "cell_line", "gene_name", "gene_id", "value")]
    colnames(mclp_rppa) <- c("feature_id", "sample_id", "gene_symbol_primary", "gene_id", "value")

    # Keep only the first gene symbol when multiple genes are listed
    mclp_rppa$gene_symbol_primary <- trimws(
    sub("///.*$", "", mclp_rppa$gene_symbol_primary)
    )

# Keep only the first gene ID, then remove the Ensembl version suffix
    mclp_rppa$gene_id <- trimws(
    sub("\\..*$", "",
        sub("///.*$", "", mclp_rppa$gene_id))
    )
    #Apply MCLP prefix to sample_id
    mclp_rppa$sample_id <- paste0("MCLP_", mclp_rppa$sample_id)
    return(mclp_rppa)
}

#Read RDS file
mclp_preclinical <- readRDS(MCLP_PATH)

#Update Pset Object if need be
mclp_preclinical <- updateObject(mclp_preclinical)
#Create filtered cell line table
filtered_cell_line_table <- filter_cell_lines(mclp_preclinical)

#Create sample table
sample_table <- create_sample_table(filtered_cell_line_table)



#Grab and format RPPA data from the MCLP pset object
rppa_table <- extract_rppa_data(mclp_preclinical, sample_table$cell_line_name)

#Write RPPA data to CSV
write.csv(rppa_table, paste0(OUTPUT_PATH, "/pre_clinical_rppa.csv"), row.names = FALSE)

#Write sample table to CSV
write.csv(sample_table, paste0(OUTPUT_PATH, "/pre_clinical_sample.csv"), row.names = FALSE)

#Write cell line table to CSV
write.csv(filtered_cell_line_table, paste0(OUTPUT_PATH, "/pre_clinical_cell_line.csv"), row.names = FALSE)
