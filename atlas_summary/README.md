# Atlas Summary

Code and data for loading the CRISPR-Cas Atlas into a PostgreSQL database and performing downstream analysis of Cas13b (Type VI-B) CRISPR systems.

## Overview

This folder handles the ETL pipeline from the CRISPR-Cas Atlas JSONL source files into a structured relational database, and contains notebooks for summarising and visualising the resulting data. Analyses include operon composition, spacer properties, Cas13b co-occurrence with other Cas genes, and protospacer location identification.

## Contents

- **`notebooks/`** — Jupyter notebooks for analysis and visualisation:
  - `atlas_summary.ipynb` — General CRISPR-Cas Atlas summary statistics.
  - `cas13b_vi_summary.ipynb` — Focused analysis of Type VI-B (Cas13b) systems.
  - `class2_cas_summary.ipynb` — Broader Class 2 Cas system comparisons.
  - `protospacer_location.ipynb` — Analysis of protospacer positions relative to target sequences.
- **`src/`** — Supporting scripts:
  - `load_data2sql.py` — ETL script to parse JSONL files and load operon data into the PostgreSQL database.
  - `extract.py` — Functions for querying and extracting data from the database.
  - `get_cas13b_cooccurrence.py` — Analysis of Cas13b co-occurrence with other Cas proteins.
  - `upload_parallel.py` — Parallel upload utility for large datasets.
  - `migrate_sql_tables.py` — Database migration utilities.
  - `identify_protospacer_loc/` — Module for identifying protospacer locations in phage genomes.
    - `main.py` — Orchestrates the pipeline: extracts Cas13b spacers (25–35 nt) from the database, runs parallel BLASTN searches against the Phagescope phage database, and compiles results with protospacer coordinates and terminase features.
    - `utils.py` — Utility functions for database queries, BLASTN execution, parsing BLAST output into structured records, fetching phage terminase features, and reorienting genomic coordinates.
- **`sql/`** — SQL assets including schema definitions, queries, views, indexes, and stored functions.

