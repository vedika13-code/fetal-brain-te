#!/usr/bin/env bash
# Step 3: splice-aware alignment to the whole genome with HISAT2, keeping up to 100 alignments per read
set -euo pipefail
source config.sh
IDX=$REF/hg38/genome
mkdir -p $WORK/bam $RESULTS/align

for S in $SAMPLES; do
  if [ -s $WORK/bam/$S.bam ]; then echo "== $S already aligned, skipping =="; continue; fi
  echo "== HISAT2 on $S (started $(date +%H:%M)) =="
  hisat2 -p $THREADS -x $IDX \
         -1 $WORK/trimmed/${S}_R1.fq.gz -2 $WORK/trimmed/${S}_R2.fq.gz \
         -k 100 --no-unal --mm \
         --summary-file $RESULTS/align/${S}.hisat2_summary.txt \
    | samtools view -b -o $WORK/bam/$S.bam -
  echo "   finished $(date +%H:%M)"
done

for S in $SAMPLES; do echo "---- $S ----"; cat $RESULTS/align/${S}.hisat2_summary.txt; done
echo "Contigs in the index that are alternate haplotypes (_alt): $(samtools view -H $WORK/bam/$(echo $SAMPLES | cut -d' ' -f1).bam | grep -c '_alt' || true)"
