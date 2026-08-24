# Preclinical Data Pipeline: End-to-End Rerun Guide

This README documents the full end-to-end process for regenerating extracted preclinical CSVs and reseeding the `sts_portal_beta` MySQL database.

The workflow is:

1. Confirm required raw input files are present.
2. Drop existing preclinical tables.
3. Run extraction scripts to generate processed CSVs.
4. Seed each dataset into MySQL.
5. Seed drug annotations from AnnotationDB.
6. Seed cell-line annotations from AnnotationDB and OncoTree.
7. Run basic validation checks.

---

## 1. Required raw input files

All raw input files must be placed in:

```bash
extraction/data/raw/preclinical/
```

Required PSet files:

```text
CCLE_2019.rds
CCLE_2015.rds
GCSI_2019.rds
PSet_PRISM.rds
Pset_CTRPv2.rds
PSet_GDSCv2.rds
```

Source notes:

```text
CCLE_2015.rds  -> ORCESTRA
CCLE_2019.rds  -> ORCESTRA
PSet_GDSCv2.rds -> ORCESTRA
GCSI_2019.rds  -> Azure
PSet_PRISM.rds -> Azure
Pset_CTRPv2.rds -> Azure
```

Required CSV files:

```text
All_PSets_sarcoma_cell_line_QC.csv
combined_datasets.csv
```

Purpose of the CSVs:

```text
All_PSets_sarcoma_cell_line_QC.csv
  Used by all extraction scripts for the final mod_tissueid == "Soft Tissue" filtering.

combined_datasets.csv
  Used by all seeding scripts to populate the datasets table.
```

Recommended check:

```bash
ls -lh extraction/data/raw/preclinical/
```

You should see the six `.rds` files and the two `.csv` files listed above.

---

## 2. Drop existing preclinical tables

Before doing a full clean reload, drop the existing preclinical tables from `sts_portal_beta`.

Open MySQL:

```bash
mysql -u YOUR_USER -p sts_portal_beta
```

Then run:

```sql
USE sts_portal_beta;

SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS pre_clinical_cell_line_info;
DROP TABLE IF EXISTS drugs;
DROP TABLE IF EXISTS pre_clinical_microarray;
DROP TABLE IF EXISTS pre_clinical_copy_number_variation;
DROP TABLE IF EXISTS pre_clinical_mutation;
DROP TABLE IF EXISTS pre_clinical_rna_seq;
DROP TABLE IF EXISTS pre_clinical_treatment_response;
DROP TABLE IF EXISTS pre_clinical_sample;
DROP TABLE IF EXISTS pre_clinical_gene;
DROP TABLE IF EXISTS pre_clinical_cell_line;
DROP TABLE IF EXISTS pre_clinical_dataset;
DROP TABLE IF EXISTS datasets;

SET FOREIGN_KEY_CHECKS = 1;
```

This removes all preclinical data and schema objects so the rerun starts cleanly.

---

## 3. Run extraction scripts

Run from the repository root.

```bash
R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_prism_preclinical.R

R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ctrpv2_preclinical.R

R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ccle_preclinical.R

R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gcsi_preclinical.R

R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gdscv2_preclinical.R
```

The extraction scripts read raw PSets from:

```bash
extraction/data/raw/preclinical/
```

and write processed CSVs into:

```text
extraction/data/proc/preclinical/PRISM/
extraction/data/proc/preclinical/CTRPv2/
extraction/data/proc/preclinical/CCLE/
extraction/data/proc/preclinical/GCSI/
extraction/data/proc/preclinical/GDSCv2/
```

Expected processed outputs vary by dataset, but include files such as:

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

Notes:

- PRISM and CTRPv2 only produce sample/cell-line/treatment-response style outputs.
- Mutation CSVs may be generated for molecular datasets, but mutation insertion is intentionally skipped during seeding.
- Final cell-line selection is based on `All_PSets_sarcoma_cell_line_QC.csv` where `mod_tissueid == "Soft Tissue"`.
- Original dataset-specific soft-tissue criteria are used only for audit/comparison outputs, not for final extraction.

---

## 4. Seed extracted datasets into MySQL

Run from the repository root after all extraction scripts complete successfully.

```bash
pixi run python -m seeding.scripts.preclinical_seed_prism

pixi run python -m seeding.scripts.preclinical_seed_ctrpv2

pixi run python -m seeding.scripts.preclinical_seed_ccle

pixi run python -m seeding.scripts.preclinical_seed_gcsi

pixi run python -m seeding.scripts.preclinical_seed_gdscv2
```

These scripts insert the processed CSV data into MySQL.

The seeders populate:

```text
datasets
pre_clinical_cell_line
pre_clinical_sample
pre_clinical_treatment_response
pre_clinical_gene
pre_clinical_rna_seq
pre_clinical_microarray
pre_clinical_copy_number_variation
```

Mutation insertion is intentionally skipped.

---

## 5. RNA-seq standardization

RNA-seq values are standardized during seeding to:

```text
log2(TPM + 1)
```

Dataset-specific behavior:

```text
CCLE:
  Source RNA-seq is already TPM.
  Seeder transforms TPM -> log2(TPM + 1).

GCSI:
  Source RNA-seq is log2(TPM + 0.001).
  Seeder transforms log2(TPM + 0.001) -> TPM -> log2(TPM + 1).

GDSCv2:
  Source RNA-seq is log2(TPM + 0.001).
  Seeder transforms log2(TPM + 0.001) -> TPM -> log2(TPM + 1).
```

---

## 6. CNV handling

CNV values are stored as gene-level log2 copy-number values.

Dataset-specific behavior:

```text
CCLE:
  Source CNV is linear copy-ratio.
  Seeder transforms value -> log2(value).
  Non-positive values are treated as missing and skipped.

GCSI:
  Source CNV already appears log2-scaled.
  Seeder stores CNV as-is.

GDSCv2:
  Source CNV already appears log2-scaled.
  Seeder stores CNV as-is.
```

---

## 7. Gene ID cleaning

Gene Ensembl IDs are cleaned by removing any suffix beginning with a period.

Example:

```text
ENSG00000002586.20_PAR_Y -> ENSG00000002586
```

The gene table uses global Ensembl IDs. If a gene ID already exists, seeders skip inserting that gene again and do not replace the existing gene name.

---

## 8. Seed the drugs table from treatment-response CIDs

Run after all treatment-response rows have been seeded.

```bash
pixi run python -m seeding.scripts.preclinical_seed_drugs_from_annotationdb --batch-size 250
```

Optional test run:

```bash
pixi run python -m seeding.scripts.preclinical_seed_drugs_from_annotationdb \
  --limit 10 \
  --batch-size 10
```

This script:

1. Queries unique non-null `cid` values from `pre_clinical_treatment_response`.
2. Calls AnnotationDB using batches of up to 250 compounds.
3. Inserts returned compound annotations into the `drugs` table.
4. Stores only the first mechanism object from the AnnotationDB `mechanisms` list.
5. Writes audit CSVs for returned and missing CIDs.

AnnotationDB route pattern:

```text
https://annotationdb.bhklab.ca/compound/many?compound=<cid>&compound=<cid>&format=json&bioassay=false&mechanism=true&toxicity=false&golden_bioassay=false
```

Audit outputs:

```text
extraction/data/proc/preclinical/drugs/drug_cids_returned_annotationdb.csv
extraction/data/proc/preclinical/drugs/drug_cids_not_found_annotationdb.csv
```

---

## 9. Seed cell-line info from Cellosaurus accessions

Run after all `pre_clinical_cell_line` rows have been seeded.

```bash
pixi run python -m seeding.scripts.preclinical_seed_cell_line_info_from_annotationdb --batch-size 250
```

Optional test run:

```bash
pixi run python -m seeding.scripts.preclinical_seed_cell_line_info_from_annotationdb \
  --limit 10 \
  --batch-size 10
```

This script:

1. Queries unique non-null Cellosaurus accessions from `pre_clinical_cell_line.accession`.
2. Calls AnnotationDB cell-line endpoint in batches.
3. Inserts returned annotations into `pre_clinical_cell_line_info`.
4. Flattens disease objects into pipe-separated fields.
5. Queries OncoTree using disease descriptions to resolve disease hierarchy.

AnnotationDB route pattern:

```text
https://annotationdb.bhklab.ca/cell_line/many?cell_lines=CVCL_1058,CVCL_1205&format=json
```

Disease flattening:

```text
disease_ids = C8971|Orphanet_99757
disease_descriptions = Embryonal rhabdomyosarcoma|Embryonal rhabdomyosarcoma
```

OncoTree behavior:

```text
For each disease description:
  1. Query OncoTree by disease name.
  2. If no hit / 404, move to the next disease description.
  3. If a hit is found, recursively query by parent code until level == 2.
  4. Store:
       first_level  = parent from level-2 response, with underscores replaced by spaces
       second_level = name from level-2 response
```

If none of the disease descriptions resolve to OncoTree, then:

```text
first_level = NULL
second_level = NULL
```

Audit outputs:

```text
extraction/data/proc/preclinical/cell_line_info/cell_line_accessions_returned_annotationdb.csv
extraction/data/proc/preclinical/cell_line_info/cell_line_accessions_not_found_annotationdb.csv
extraction/data/proc/preclinical/cell_line_info/cell_line_info_oncotree_unresolved.csv
```

---

## 10. Full clean run command sequence

After dropping the tables in MySQL, run:

```bash
R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_prism_preclinical.R
R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ctrpv2_preclinical.R
R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_ccle_preclinical.R
R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gcsi_preclinical.R
R_MAX_VSIZE=64Gb pixi run Rscript extraction/scripts/preclinical/extract_gdscv2_preclinical.R

pixi run python -m seeding.scripts.preclinical_seed_prism
pixi run python -m seeding.scripts.preclinical_seed_ctrpv2
pixi run python -m seeding.scripts.preclinical_seed_ccle
pixi run python -m seeding.scripts.preclinical_seed_gcsi
pixi run python -m seeding.scripts.preclinical_seed_gdscv2

pixi run python -m seeding.scripts.preclinical_seed_drugs_from_annotationdb --batch-size 250

pixi run python -m seeding.scripts.preclinical_seed_cell_line_info_from_annotationdb --batch-size 250
```

---

## 11. Quick validation queries

Run these in MySQL after seeding.

```sql
USE sts_portal_beta;

SELECT * FROM datasets;

SELECT COUNT(*) AS cell_line_count FROM pre_clinical_cell_line;
SELECT COUNT(*) AS sample_count FROM pre_clinical_sample;
SELECT COUNT(*) AS treatment_response_count FROM pre_clinical_treatment_response;
SELECT COUNT(*) AS gene_count FROM pre_clinical_gene;
SELECT COUNT(*) AS rnaseq_count FROM pre_clinical_rna_seq;
SELECT COUNT(*) AS microarray_count FROM pre_clinical_microarray;
SELECT COUNT(*) AS cnv_count FROM pre_clinical_copy_number_variation;
SELECT COUNT(*) AS drug_count FROM drugs;
SELECT COUNT(*) AS cell_line_info_count FROM pre_clinical_cell_line_info;
```

Check RNA-seq range:

```sql
SELECT
  MIN(value) AS min_rnaseq,
  AVG(value) AS mean_rnaseq,
  MAX(value) AS max_rnaseq
FROM pre_clinical_rna_seq;
```

RNA-seq should be non-negative because it is stored as `log2(TPM + 1)`.

Check CNV range:

```sql
SELECT
  MIN(value) AS min_cnv,
  AVG(value) AS mean_cnv,
  MAX(value) AS max_cnv
FROM pre_clinical_copy_number_variation;
```

CNV can be negative because it is stored as log2 copy-number value.

Check missing drug annotations:

```bash
cat extraction/data/proc/preclinical/drugs/drug_cids_not_found_annotationdb.csv
```

Check unresolved OncoTree mappings:

```bash
cat extraction/data/proc/preclinical/cell_line_info/cell_line_info_oncotree_unresolved.csv
```

---

## 12. Notes on `--replace`

For a full clean reload, do not use `--replace`; drop the tables and run the seeders normally.

Use `--replace` only when reloading one dataset into an existing database.

Example:

```bash
pixi run python -m seeding.scripts.preclinical_seed_ccle --replace
```

This deletes and reinserts only the CCLE dataset rows. It does not reset every preclinical table.

---

## 13. Common troubleshooting

### Missing raw PSet file

If an extractor fails because a PSet is missing, confirm the exact file name and location:

```bash
ls -lh extraction/data/raw/preclinical/
```

### MySQL foreign key errors

For a full rerun, make sure the tables were dropped in dependency order with:

```sql
SET FOREIGN_KEY_CHECKS = 0;
```

and then re-enabled with:

```sql
SET FOREIGN_KEY_CHECKS = 1;
```

### AnnotationDB missing CIDs or accessions

Missing API results are written to audit CSVs. They do not necessarily indicate a pipeline failure.

### OncoTree 404 responses

A 404 from OncoTree name search means no tumor type matched that disease description. The seeder should move to the next disease description. If none resolve, `first_level` and `second_level` are stored as `NULL`.
