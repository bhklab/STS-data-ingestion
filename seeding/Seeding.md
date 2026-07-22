# Seeding

This directory contains all of the scripts needed for creating tables for our MySQL database and also seeding them from the csvs extracted from the MAE's. 

## Scripts

`/scripts/clinical_seeding_coordinator.py` is used to create and seed the clinical tables, `/scripts/preclinical_seeding_coordinator.py` is used to create and seed the pre clinical tables.

`/scripts/seeding_coordinator_engine.py` is referenced by both the clinical and preclinical seeding scripts to create a SQL Alchemy engine attached to our STS MySQL database to allow for easy creation/manipulation of tables.

## Models

There is only a singular relevant model file for this seeding coordinator and that is the `/models/tables.py` file. This file contains all of the tables SQL Alchemy table definitions for the STS database. 

## data

Contains the raw csv files needed for seeding out clinical and preclinical MySQL tables.

## How to run seeding files

Due to the reliance of relatively located files, the seeding coordinators must be run as modules. For example the `preclinical_seeding_coordinator.py` file can be run by executing the following command:

```bash
pixi run python -m seeding.scripts.preclinical_seeding_coordinator
```
