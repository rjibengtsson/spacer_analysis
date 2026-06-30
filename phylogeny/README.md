# Phylogeny

Phylogenetic analysis and characterisation of Cas13b protein sequences from the CRISPR-Cas Atlas.

## Overview

This folder contains sequence data, clustering results, and analysis code for exploring the diversity of Cas13b orthologs. Key steps include parsing CD-HIT clusters of Cas13b proteins, retrieving spacer/guide sequences for representative clusters, multiple sequence alignment, and identification of conserved functional motifs (HEPN domain).

## Contents

- **`crisprcasatlas.ipynb`** — Main analysis notebook. Parses CD-HIT clustering output, retrieves Cas13b sequences and guide RNAs from the CRISPR-Cas Atlas, and summarises sequence diversity across clusters.
- **`pspCas13b_cluster.list`** — List of 128 accessions belonging to the PspCas13b cluster from CD-HIT.
- **`pspCas13b_cluster_aligned.fasta`** — Multiple sequence alignment of the 128-sequence PspCas13b cluster, for use in phylogenetic tree construction.
- **`PspCas13b.fasta`** — Protein sequence for PspCas13b (*Prevotella* sp.).
- **`BzCas13b.fasta`** — Protein sequence for BzCas13b (*Bergeyella zoohelcum*).
- **`prevotella_pectinovora_cas13b.fasta`** — Cas13b protein sequence from *Prevotella pectinovora*.
- **`P5-125/`** — GenBank file (`GCF_000833995.1`) for the reference genome associated with PspCas13b.
- **`scripts/`**
  - `HEPN_finder.py` — Scans Cas13b protein sequences for HEPN domain motifs (pattern `R.{4}H`) using regex and reports match positions.
