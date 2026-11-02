#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 4 ]] || { echo "Usage: $0 config.env SAMPLE reads.fastq.gz OUTDIR" >&2; exit 2; }
source "$1"; sample=$2; input=$3; out=$4; mkdir -p "$out"
: "${SMRNA_INDEX:?}" "${SMRNA_GENOME_INDEX:?}" "${SMRNA_ADAPTER:?}"
threads=${THREADS:-4}; mode=${SMRNA_GENOME_MODE:-unmatched}
[[ $mode == all || $mode == unmatched ]] || { echo 'SMRNA_GENOME_MODE must be all or unmatched' >&2; exit 2; }
p="$out/$sample"
cutadapt -j "$threads" -a "$SMRNA_ADAPTER" -q 20 -m 18 -M 45 --discard-untrimmed -o "$p.trim.fq.gz" "$input" > "$p.cutadapt.log"
gzip -cd "$p.trim.fq.gz" | awk 'NR%4==2' | LC_ALL=C sort | uniq -c | awk '{print ">"$2"_"$1"\n"$2}' > "$p.collapsed.fa"
bowtie2 -p "$threads" -f -x "$SMRNA_INDEX" -U "$p.collapsed.fa" --un "$p.unmatched.fa" -S "$p.smallrna.sam" 2> "$p.smallrna.log"
if [[ $mode == all ]]; then genome_input="$p.collapsed.fa"; else genome_input="$p.unmatched.fa"; fi
bowtie2 -p "$threads" -f --very-sensitive -L 16 -x "$SMRNA_GENOME_INDEX" -U "$genome_input" -S "$p.genome.sam" 2> "$p.genome.log"
# Classify small-RNA matches separately; 'all' mode can produce duplicate sequence rows with both annotation and genomic hit.
{
 printf 'sequence\tcount\tclass\tchrom\tpos5\tstrand\tmapq\n'
 samtools view -F 0x904 "$p.smallrna.sam" | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print a[1],a[2],"smRNA",$3,".",".",$5}'
 samtools view -F 0x900 "$p.genome.sam" | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); seq=a[1];cnt=a[2];if (int($2/4)%2) {print seq,cnt,"unmapped",".",".",".",0;next} nm=-1;for(i=12;i<=NF;i++)if($i~/^NM:i:/)nm=substr($i,6)+0; cls=(nm==0?"genome_0mm":nm==1?"genome_1mm":"genome_2plus"); rev=int($2/16)%2; cig=$6;rl=0;while(match(cig,/^[0-9]+[MIDNSHP=X]/)){len=substr(cig,1,RLENGTH-1)+0;op=substr(cig,RLENGTH,1);if(op~/[MDN=X]/)rl+=len;cig=substr(cig,RLENGTH+1)};pos=rev?$4+rl-1:$4;print seq,cnt,cls,$3,pos,(rev?"-":"+"),$5}'
} | gzip -c > "$p.seq_annot.tsv.gz"
echo "Output: $p.seq_annot.tsv.gz"
