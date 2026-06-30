# Phage

Analysis of spacer sequences against phage genomes to identify protospacer targets and characterise the nucleotide composition of integrated versus non-integrated phage sequences.

## Overview

Spacer sequences from CRISPR arrays are searched against the Phagescope phage database using BLASTN. The resulting hits are used to identify protospacer regions in phage genomes and compare nucleotide properties between spacer-targeted regions and the broader phage genome background.

## Contents

- **`phage_genome_analysis.ipynb`** — Main analysis notebook. Covers BLASTN searches of spacer sequences against the Phagescope database, protospacer identification, nucleotide content analysis, and visualisation of sequence logos and distributions.
- **`spacer_phageblastn_results.csv`** — BLASTN results for all spacers searched against the Phagescope database.
- **`test_env/`** — Development and testing scripts:
  - `blastn_phagescope.py` — Script to run BLASTN searches of filtered spacer sequences against the Phagescope database in parallel.
  - `extract_orf.py` — Script to extract protospacer–phage genome intersections for a target species (e.g., *Prevotella pectinovora*), linking spacer hits to phage genome ORF coordinates.
  - `phage_genomes.fasta` / `phage_sequences_*.fasta` / `Prevotella_pectinovora_*.fasta` — Phage genome and sequence files used for testing.
  - `protospacer_phage_intersection_Prevotella_pectinovora.csv` — Test output intersecting protospacer hits with phage genome annotations.
