# Extraction

Suite of tools for extracting various data layers from clinical and preclinical MAE's to be seeded into the STS MySQL database.

## Directory Structure

### data

Container for all raw data objects, extracted CSVs, and database backups:
- `/data/raw`: Contains raw (MAEs, PSets) and semi-processed data files.
- `/data/proc`: Contains the extracted CSV files ready to be piped into seeding scripts.
- `/data/backups`: Contains timestamped CSV database dumps and `backup_tables.py` utility.

### scripts

Scripts used for extracting and/or formatting extracted data into the csvs needed for MySQL ingestion