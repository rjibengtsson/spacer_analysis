# Modules

Shared Python library used across the project. Each module provides reusable utilities for a specific domain of the analysis pipeline.

## Contents

- **`analysis_utils.py`** — Core analysis utilities including spacer/gRNA sequence processing, nucleotide content calculation, per-position nucleotide frequency counting, Fisher's exact test helpers, and dataclasses for CRISPR array candidates (`ArrayCandidate`) and phage elements (`PhageElement`).
- **`atlas_handle.py`** — Functions for parsing and handling data from the CRISPR-Cas Atlas JSONL files, including CRISPR subtype counting and spacer retrieval.
- **`blast_utils.py`** — Utilities for running BLASTN searches, formatting spacer sequences as FASTA input, and processing BLAST output; supports parallel execution.
- **`genome_utils.py`** — Functions for fetching genome FASTA files from NCBI, parsing GFF annotations, and working with contig sequences (`CasContig` dataclass).
- **`annotation.py`** — Genome annotation utilities including GFF parsing via `gffutils` and downloading genome files from NCBI FTP paths.
- **`database_utils.py`** — SQLAlchemy-based utilities for interacting with the PostgreSQL database, including dataclasses for database entities (e.g., `Allthebacteria`).
- **`crispr_detecion.py`** — Wrapper around CRISPR detection tools (`CrisprDetection` class) for predicting CRISPR arrays in genome sequences via subprocess calls.
- **`visualization.py`** — Plotting utilities (currently in development).
- **`test_atlas_handle.py`** — Unit tests for `atlas_handle.py`.
- **`test_files/`** — Test input files for unit tests (e.g., spacer lists).
