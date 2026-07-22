# Extraction

Suite of tools for extracting various data layers from clinical and preclinical MAE's to be seeded into the STS MySQL database.

## Directory Structure

### data

Container for all raw data objects and their extracted csvs. `/data/raw` contains raw (MAE's) and semi processed data files. `/data/proc` contains the extracted csvs files that are ready to be piped into the seeding scripts in 

### scripts

Scripts used for extracting and/or formatting extracted data into the csvs needed for MySQL ingestion