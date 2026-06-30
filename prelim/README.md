# Prelim

Preliminary CRISPR-Cas13b spacer analysis on a small dataset to establish the analysis pipeline and explore initial genomic signals from spacer sequences.

## Overview

This folder contains the end-to-end pipeline for fetching genome data from NCBI, identifying Cas13b loci, running CRISPR array prediction, and performing exploratory analysis of spacer nucleotide content and properties. Results from this stage feed into the broader analyses in `atlas_summary/` and `phage/`.

## Contents

- **`premlinary_analysis.ipynb`** — Main exploratory notebook. Analyses spacer sequences from a preliminary dataset, including nucleotide content, length distributions, and sequence logo generation.
- **`pipeline.py`** — Core pipeline script orchestrating:
  - `retrieve_data_from_db` — Query the `cas13_bacterial_db` PostgreSQL database for candidate genomes.
  - `fetch_data` — Download genome FASTA, GFF, and GenBank files from NCBI FTP for each biosample.
  - `cas13b_search` — Identify Cas13b loci in downloaded genomes using annotation data.
  - `run_crispr_prediction` — Run CRISPR array prediction (via `CrisprDetection`) on genome sequences.
  - `upload_result` — Upload CRISPR prediction results back to the database.
- **`runCRISPRidentify.py`** — Standalone script to run CRISPR array prediction on a set of biosample accessions using `CrisprDetection`.
- **`analysis/`** — Processed CRISPR analysis outputs:
  - `preliminary_crisprs.csv` / `all_spacers.csv` — CRISPR array and spacer data for the full preliminary dataset.
  - `SAMN03382583_crisprs.csv` / `SAMN03382583_spacers.csv` — Array and spacer data for a specific biosample.
  - `repeat_sequences.csv` — Repeat sequences identified in CRISPR arrays.
  - `spacer_nt_content.csv` — Nucleotide content calculations for spacer sequences.
- **`data/`** — Input metadata:
  - `ncbi_assembly.csv` — NCBI assembly metadata for the preliminary dataset.
  - `prelim_dataset_final.csv` — Final curated list of samples used in the analysis.
  - `P5-125_seqs.csv` — Sequence data for the P5-125 reference strain.
- **`mydb.db`** — Local SQLite database used during development.
