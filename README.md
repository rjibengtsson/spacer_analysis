# Spacer Sequence Analysis

Analysis of CRISPR spacer sequences from Cas13b (Type VI-B) systems, investigating molecular signatures that underlie functional guide RNA potency.

## Overview

This repository brings together a pipeline for discovering, characterising, and interpreting Cas13b spacer sequences derived from large-scale bacterial genomic data. Starting from the CRISPR-Cas Atlas and public NCBI assemblies, the workflow covers:

1. **Genome retrieval and CRISPR detection** — downloading bacterial genomes and identifying CRISPR arrays using CRISPRidentify and CRISPRCasFinder.
2. **Database ingestion** — loading parsed operon and spacer data into a PostgreSQL database for scalable querying.
3. **Spacer–phage mapping** — BLASTN searches of spacer sequences against the Phagescope phage database to identify protospacer targets in phage genomes.
4. **Sequence analysis** — nucleotide content, length distributions, positional composition, and sequence logo generation for spacer and protospacer sequences.
5. **Phylogenetics** — clustering and multiple sequence alignment of Cas13b protein sequences to explore evolutionary diversity.
6. **Visualisation and reporting** — figures and summary tables for conference presentations and ongoing analysis.

## Repository Structure

| Folder | Description | Status |
|---|---|---|
| `prelim/` | Preliminary pipeline: genome fetch, Cas13b locus identification, CRISPR prediction, initial spacer analysis | Complete |
| `atlas_summary/` | ETL from CRISPR-Cas Atlas into PostgreSQL; Cas13b spacer BLASTN; protospacer location analysis | Active |
| `phage/` | Phage genome analysis: spacer-phage BLASTN, protospacer nucleotide characterisation | Complete |
| `phylogeny/` | Cas13b protein clustering, multiple sequence alignment, HEPN domain identification | Complete |
| `abacbs2025/` | Figures and analysis for the ABACBS 2025 conference poster | Complete |
| `modules/` | Shared Python library (BLAST utilities, database access, genome utilities, CRISPR detection wrappers) | - |
| `tools/` | Third-party tools: BLAST, CD-HIT, CRISPRidentify, CRISPRCasFinder, bedtools | - |

See each subfolder's own README for details on the analysis and how to reproduce it.

## Setup

Dependencies are managed with conda. To create the environment:

```bash
conda env create -f environment.yaml
conda activate modules
```

Key dependencies: Python 3.10, BLAST 2.14, Biopython, pandas, numpy, plotnine, logomaker, SQLAlchemy, psycopg2.

A PostgreSQL database (`crisprcas_atlas`) is required for the `atlas_summary/` and `prelim/` pipelines. Connection credentials are stored in `modules/db_credentials.json` (not tracked in version control).
