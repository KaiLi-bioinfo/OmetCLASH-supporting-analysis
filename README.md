# OmetCLASH: Supporting Analysis of MIWI Small RNA-seq and LACE-seq Data

This repository contains scripts for processing MIWI-associated small RNA-seq and LACE-seq data.

Small RNA-seq

Reads are adapter-trimmed and collapsed into unique sequences with read counts. Known small RNAs (rRNA, tRNA, snoRNA, snRNA and miRNA) are identified using Bowtie2 and RNAcentral references. Sequences are then aligned to a mouse genome reference. The included summary script reports read lengths and nucleotide preferences (1U).

LACE-seq

Reads are trimmed to remove adapters, low-quality bases and homopolymer tails. After filtering rRNA-matching reads, the remaining reads are aligned to the mouse GRCm39 genome using Bowtie2. Peaks are called with Piranha, annotated with HOMER, and visualized using CPM-normalized coverage tracks from deepTools.

Scripts

* small_rna/run.sh — small RNA-seq preprocessing and alignment
* small_rna/summary.py — basic small RNA summaries
* lace_seq/run.sh — LACE-seq preprocessing, alignment and peak calling
* config.example.env — reference paths and settings
* environment.yml — software dependencies

Note: These scripts were reorganized from the analysis steps provided for the project. They are not a complete copy of the original analysis pipeline; additional piRNA QC and IgG-background filtering steps are not included. Reference settings should be checked against the data being analyzed.
