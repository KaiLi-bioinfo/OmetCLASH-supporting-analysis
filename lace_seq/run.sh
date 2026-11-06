#!/usr/bin/env bash
#LACE-seq data processing

#Only Read 1 (R1) was used for analysis because the LACE-seq

#inserts were short and Read 2 (R2) consisted predominantly

#f poly(G) sequences (GGGGGGGGG), providing no useful

#equence information for genomic alignment.
#
set -euo pipefail
[[ $# -eq 4 ]] || { echo "Usage: $0 config.env SAMPLE reads.R1.fastq.gz OUTDIR" >&2; exit 2; }
source "$1"; sample=$2; input=$3; out=$4; mkdir -p "$out"; : "${LACE_RRNA_INDEX:?}" "${LACE_GENOME_INDEX:?}"
threads=${THREADS:-4}; p="$out/$sample"
cutadapt -j "$threads" -q 30,0 -a ATCTCGTATGCCGTCTTCTGCTT -m 18 --max-n 0.25 --trim-n -o "$p.step1.fq.gz" "$input" > "$p.step1.log"
cutadapt -j "$threads" -u 4 -a 'A{9}' -a 'G{9}' -m 17 -n 2 -o "$p.step2.fq.gz" "$p.step1.fq.gz" > "$p.step2.log"
cutadapt -j "$threads" -a 'A{9}' -a 'G{9}' -a 'T{9}' -m 17 -n 2 -o "$p.step3.fq.gz" "$p.step2.fq.gz" > "$p.step3.log"
cutadapt -j "$threads" -q 30,0 -a GCCGTCTTCT -m 18 --max-n 0.25 --trim-n -o "$p.clean.fq.gz" "$p.step3.fq.gz" > "$p.step4.log"
bowtie2 -p "$threads" --very-sensitive -x "$LACE_RRNA_INDEX" -U "$p.clean.fq.gz" --un-gz "$p.no_rRNA.fq.gz" -S /dev/null 2> "$p.rrna.log"
bowtie2 -p "$threads" -k 1 -x "$LACE_GENOME_INDEX" -U "$p.no_rRNA.fq.gz" 2> "$p.map.log" | samtools sort -@ "$threads" -o "$p.sorted.bam" -
samtools index "$p.sorted.bam"
bedtools bamtobed -i "$p.sorted.bam" > "$p.bed"
Piranha "$p.bed" -s -p 0.0001 -b 50 -d ZeroTruncatedNegativeBinomial -o "$p.peaks.bed"
bamCoverage -b "$p.sorted.bam" -o "$p.CPM.bw" --normalizeUsing CPM --binSize 25 -p "$threads"
if [[ -n ${HOMER_GENOME:-} ]] && command -v annotatePeaks.pl >/dev/null 2>&1; then
 annotatePeaks.pl "$p.peaks.bed" "$HOMER_GENOME" -annStats "$p.annStats.txt" -cpu "$threads" > "$p.peaks.annotation.tsv"
fi
rm -f "$p.step1.fq.gz" "$p.step2.fq.gz" "$p.step3.fq.gz"
echo "Outputs: $p.peaks.bed $p.CPM.bw"
