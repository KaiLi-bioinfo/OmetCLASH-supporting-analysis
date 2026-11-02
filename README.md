# MIWI RIP small RNA-seq and LACE-seq analysis

Reproducible, configurable analysis workflows adapted from the analysis commands shared in the manuscript preparation conversation.

> **Important:** This repository is a clean-room reconstruction from pasted code excerpts, not an exact archive of the original executed scripts. Validate parameters and rerun on a test sample before citing it as the deposited analysis code.

## Contents
- `small_rna/run.sh`: adapter trimming, sequence collapsing, RNAcentral screening and genome alignment.
- `small_rna/summary.py`: length, 1U/10A, and RPM summaries from sequence tables.
- `lace_seq/run.sh`: LACE-seq R1 trimming, rRNA depletion, genome mapping, peak calling, coverage tracks and optional annotation.
- `config.example.env`: paths and analysis settings.
- `environment.yml`: software dependencies.

## Setup
```bash
conda env create -f environment.yml
conda activate miwi-analysis
cp config.example.env config.env
# Edit config.env to point to local FASTA/Bowtie2 indices
```

## Small RNA-seq
```bash
bash small_rna/run.sh config.env sample_id reads.fastq.gz output/small_rna
python small_rna/summary.py output/small_rna/sample_id.seq_annot.tsv.gz output/small_rna/sample_id.metrics.tsv
```

## LACE-seq
```bash
bash lace_seq/run.sh config.env sample_id reads.R1.fastq.gz output/lace_seq
```

## Analysis provenance and caveats
- The pasted original small-RNA preprocessing used an **mm10 + ERCC92** index and sent **only sequences not matched to RNAcentral** to the genome aligner. The later requested manuscript description (GRCm39, re-align *all* sequences) is **not identical**. The reconstructed workflow exposes `SMRNA_GENOME_MODE=all|unmatched` and uses a configurable genome index; do not claim the manuscript workflow was executed without checking original output provenance.
- The original LACE-seq script mapped to **GRCm39**, removed human/mouse rRNA matches, used Bowtie2 `-k 1` (not unique mapping), Piranha 50-nt bins, zero-truncated negative binomial, P < 1e-4, and bamCoverage CPM, 25-bp bins.
- `small_rna/summary.py` is a compact QC summary, **not** a replacement for the original extensive `02_pirna_qc.py` or `03_base_percentage.py` figures/ping-pong analysis. Upload the originals to include those analyses verbatim after code review.
- LACE-seq IgG subtraction is **not included**: the supplied excerpts showed an `rmIgG` filename but no actual subtraction algorithm. Do not describe background subtraction as performed by this repository.
- Bowtie2's `-k 1` does not guarantee uniquely mapped reads. `NM<=1` is an optional post-alignment filter for the small-RNA table.

## Data
No FASTQ files, genome indices, personal paths or unpublished data are bundled. Supply your own reference indices and document the versions used.

## Citation
Cite cutadapt, Bowtie2, SAMtools, BEDTools, Piranha, HOMER and deepTools as appropriate. See `docs/software.md`.

## License
MIT (see LICENSE).
