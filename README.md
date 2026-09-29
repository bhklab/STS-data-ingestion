## Preclinical Data Pipeline

The preclinical pipeline extracts soft-tissue cell-line data from PharmacoSet/RDS sources, normalizes the extracted data into a common schema, and seeds the processed data into MySQL.

The datasets currently supported are:

- PRISM
- CTRPv2
- CCLE
- GCSI
- GDSCv1
- GDSCv2

The general workflow is:

```text
Raw RDS / metadata files
        ↓
Dataset-specific R extraction
        ↓
Processed preclinical CSVs
        ↓
Dataset-specific Python seeders
        ↓
MySQL preclinical tables
        ↓
Gene-name normalization
        ↓
Drug + cell-line annotation enrichment
```

### 1. Install the environment

From the repository root:

```bash
pixi install
```

All commands below should also be run from the repository root.

---

### 2. Required raw preclinical files

Raw preclinical inputs are stored under:

```text
extraction/data/raw/preclinical/
```

The complete extraction currently expects the following RDS files (ensure same naming convention):

```text
CCLE_2015.rds
CCLE_2019.rds
GCSI_2019.rds
PSet_PRISM.rds
Pset_CTRPv2.rds
PSet_GDSCv1.rds
PSet_GDSCv2.rds
```

The following metadata/mapping files are also required:

```text
All_PSets_sarcoma_cell_line_QC.csv
combined_datasets.csv
biomart_genes.csv
```

`All_PSets_sarcoma_cell_line_QC.csv` is used to define the final soft-tissue cell-line set. Dataset-specific extraction scripts match their samples to rows where:

```text
mod_tissueid == "Soft Tissue"
```

`combined_datasets.csv` contains metadata used to populate the `datasets` table.

`biomart_genes.csv` provides Ensembl gene ID / gene-name mappings and is also used when normalizing the final `pre_clinical_gene` table.

A quick input check can be run with:

```bash
ls -lh extraction/data/raw/preclinical/
```

---

### 3. Optional: clean the existing preclinical database

For a complete clean reload, remove the existing preclinical tables first.

```sql
USE sts_portal_beta;

SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS pre_clinical_cell_line_info;
DROP TABLE IF EXISTS drugs;

DROP TABLE IF EXISTS pre_clinical_mirna;
DROP TABLE IF EXISTS pre_clinical_rppa;
DROP TABLE IF EXISTS pre_clinical_methylation_tss_1kb;
DROP TABLE IF EXISTS pre_clinical_massspec_intensity;

DROP TABLE IF EXISTS pre_clinical_mutation;
DROP TABLE IF EXISTS pre_clinical_rna_seq;
DROP TABLE IF EXISTS pre_clinical_microarray;
DROP TABLE IF EXISTS pre_clinical_copy_number_variation;

DROP TABLE IF EXISTS pre_clinical_treatment_response;
DROP TABLE IF EXISTS pre_clinical_sample;
DROP TABLE IF EXISTS pre_clinical_gene;
DROP TABLE IF EXISTS pre_clinical_cell_line;

DROP TABLE IF EXISTS pre_clinical_dataset;
DROP TABLE IF EXISTS datasets;

SET FOREIGN_KEY_CHECKS = 1;
```

The dataset seeders recreate the required tables through the SQLAlchemy models.

For a clean database reload, do **not** use `--replace`.

---

### 4. Run all preclinical extractions

The extraction scripts read the raw RDS files and produce normalized CSVs under:

```text
extraction/data/proc/preclinical/
```

On macOS, the following invocation also prevents R from using packages from the user's global R library and increases the available R vector heap:

```bash
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_prism_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ctrpv2_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ccle_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gcsi_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gdscv1_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gdscv2_preclinical.R
```

The processed output directories are:

```text
extraction/data/proc/preclinical/PRISM/
extraction/data/proc/preclinical/CTRPv2/
extraction/data/proc/preclinical/CCLE/
extraction/data/proc/preclinical/GCSI/
extraction/data/proc/preclinical/GDSCv1/
extraction/data/proc/preclinical/GDSCv2/
```

Depending on the dataset, generated files include:

```text
pre_clinical_cell_line.csv
pre_clinical_sample.csv
pre_clinical_treatment_response.csv
pre_clinical_gene.csv
pre_clinical_rna_seq.csv
pre_clinical_microarray.csv
pre_clinical_copy_number_variation.csv
pre_clinical_mutation.csv
```

CCLE additionally generates:

```text
pre_clinical_mirna.csv
pre_clinical_rppa.csv
pre_clinical_methylation_tss_1kb.csv
pre_clinical_massspec_intensity.csv
```

The CCLE feature profiles are sourced from:

```text
mirna.mimat_expression
proteomics.rppa
methylation.tss_1kb
proteomics.massspec_intensity
```

CCLE uses both the 2019 and 2015 PSets. The 2019 PSet provides the primary molecular profiles, while the 2015 PSet is used for the legacy microarray data and treatment response data.

> **CCLE extractor note:** the finalized extractor containing miRNA, RPPA, methylation, and mass-spec support should be the version stored as `extract_ccle_preclinical.R`. If the finalized file is still named `extract_ccle_preclinical_with_rppa_methylation_massspec.R`, run that filename instead.

---

### 5. Seed all datasets into MySQL

After all extraction scripts complete successfully, seed each dataset:

```bash
pixi run python -m seeding.scripts.preclinical_seed_prism
pixi run python -m seeding.scripts.preclinical_seed_ctrpv2
pixi run python -m seeding.scripts.preclinical_seed_ccle
pixi run python -m seeding.scripts.preclinical_seed_gcsi
pixi run python -m seeding.scripts.preclinical_seed_gdscv1
pixi run python -m seeding.scripts.preclinical_seed_gdscv2
```

The seeders populate the shared preclinical schema, including:

```text
datasets
pre_clinical_cell_line
pre_clinical_sample
pre_clinical_treatment_response
pre_clinical_gene

pre_clinical_rna_seq
pre_clinical_microarray
pre_clinical_copy_number_variation
pre_clinical_mutation
```

CCLE also populates:

```text
pre_clinical_mirna
pre_clinical_rppa
pre_clinical_methylation_tss_1kb
pre_clinical_massspec_intensity
```

Mutation loading behavior is:

```text
CCLE     -> loaded
GDSCv1   -> loaded
GDSCv2   -> loaded
GCSI     -> intentionally skipped
```

Mutation values are standardized to binary values:

```text
wt                  -> 0
other mutation call -> 1
```

RNA-seq values are normalized by the dataset seeders to the shared:

```text
log2(TPM + 1)
```

representation where applicable.

For the CCLE miRNA, RPPA, methylation, and mass-spec layers, rows that cannot be mapped to a valid Ensembl `gene_id` are intentionally skipped rather than inserted with a null gene ID.

---

### 6. Normalize gene names

Each molecular dataset may contribute Ensembl IDs with incomplete, inconsistent, or missing gene names.

After **all six dataset seeders have finished**, normalize the shared `pre_clinical_gene` table using the local BioMart snapshot:

```bash
pixi run python -m seeding.scripts.update_gene_table_names
```

The script uses:

```text
extraction/data/raw/preclinical/biomart_genes.csv
```

and expects the BioMart columns:

```text
Gene stable ID
Gene name
```

Ensembl version/postfix suffixes are normalized before matching. For example:

```text
ENSG00000002586.20_PAR_Y -> ENSG00000002586
ENSG00000123456.7        -> ENSG00000123456
```

The update is deliberately run only after all dataset seeders so that `pre_clinical_gene` already contains the union of gene IDs contributed by CCLE, GCSI, GDSCv1, and GDSCv2.

Audit files are produced for updated genes, already-correct names, and IDs not found in the BioMart snapshot.

---

### 7. Seed drug annotations

Once all treatment-response rows are available, populate the `drugs` table from AnnotationDB:

```bash
pixi run python -m seeding.scripts.preclinical_seed_drugs_from_annotationdb --batch-size 250
```

The script collects the unique non-null compound CIDs from `pre_clinical_treatment_response`, queries AnnotationDB in batches, and upserts the returned annotations into `drugs`.

AnnotationDB currently limits this seeder to a maximum batch size of 250.

---

### 8. Seed cell-line annotations

After all datasets have populated `pre_clinical_cell_line`, enrich the unique Cellosaurus accessions:

```bash
pixi run python -m seeding.scripts.preclinical_seed_cell_line_info_from_annotationdb --batch-size 250
```

This step:

```text
pre_clinical_cell_line.accession
        ↓
AnnotationDB
        ↓
cell-line metadata
        ↓
OncoTree disease resolution
        ↓
pre_clinical_cell_line_info
```

The seeder also writes audit CSVs for Cellosaurus accessions not returned by AnnotationDB and disease descriptions that could not be resolved through OncoTree.

---

### 9. Validate the load

After the complete pipeline finishes, basic row counts can be checked with:

```sql
USE sts_portal_beta;

SELECT * FROM datasets;

SELECT COUNT(*) AS cell_lines
FROM pre_clinical_cell_line;

SELECT COUNT(*) AS samples
FROM pre_clinical_sample;

SELECT COUNT(*) AS treatment_responses
FROM pre_clinical_treatment_response;

SELECT COUNT(*) AS genes
FROM pre_clinical_gene;

SELECT COUNT(*) AS rna_seq
FROM pre_clinical_rna_seq;

SELECT COUNT(*) AS microarray
FROM pre_clinical_microarray;

SELECT COUNT(*) AS cnv
FROM pre_clinical_copy_number_variation;

SELECT COUNT(*) AS mutations
FROM pre_clinical_mutation;

SELECT COUNT(*) AS mirna
FROM pre_clinical_mirna;

SELECT COUNT(*) AS rppa
FROM pre_clinical_rppa;

SELECT COUNT(*) AS methylation_tss_1kb
FROM pre_clinical_methylation_tss_1kb;

SELECT COUNT(*) AS massspec_intensity
FROM pre_clinical_massspec_intensity;

SELECT COUNT(*) AS drugs
FROM drugs;

SELECT COUNT(*) AS cell_line_info
FROM pre_clinical_cell_line_info;
```

The datasets table should contain:

```text
PRISM
CTRPv2
CCLE
GCSI
GDSCv1
GDSCv2
```

---

### 10. Complete end-to-end command sequence

For a normal clean run after the old database tables have been dropped, the complete command sequence is:

```bash
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_prism_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ctrpv2_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ccle_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gcsi_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gdscv1_preclinical.R
R_LIBS_USER=NULL R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gdscv2_preclinical.R

pixi run python -m seeding.scripts.preclinical_seed_prism
pixi run python -m seeding.scripts.preclinical_seed_ctrpv2
pixi run python -m seeding.scripts.preclinical_seed_ccle
pixi run python -m seeding.scripts.preclinical_seed_gcsi
pixi run python -m seeding.scripts.preclinical_seed_gdscv1
pixi run python -m seeding.scripts.preclinical_seed_gdscv2

pixi run python -m seeding.scripts.update_gene_table_names

pixi run python -m seeding.scripts.preclinical_seed_drugs_from_annotationdb --batch-size 250
pixi run python -m seeding.scripts.preclinical_seed_cell_line_info_from_annotationdb --batch-size 250
```

The ordering of the final three operations is important:

```text
dataset seeders
    ↓
update_gene_table_names
    ↓
drug enrichment
    ↓
cell-line enrichment
```

`update_gene_table_names` should not be run before the molecular dataset seeders because additional datasets can introduce additional Ensembl IDs into `pre_clinical_gene`.

---

### 11. Reloading a single dataset

The dataset seeders support `--replace` when an individual dataset needs to be regenerated without rebuilding the entire preclinical database.

For example:

```bash
pixi run python -m seeding.scripts.preclinical_seed_ccle --replace
```

For a **full clean rebuild**, dropping the preclinical tables and running the normal seeders without `--replace` is preferred.
