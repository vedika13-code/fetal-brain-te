#!/usr/bin/env bash
# Step 5: gene counts and TE counts with featureCounts (reverse-stranded library: -s 2)
# TE counts 4 ways: {subfamily, locus} x {unique reads only, multimappers included}.
# --fraction is used in BOTH modes so reads overlapping two TEs are split the same way;
# the only difference between uniq and multi is then -M (multimappers).
set -euo pipefail
source config.sh
mkdir -p $WORK/counts $RESULTS/counts
STRAND_N=2

for S in $SAMPLES; do
  B=$WORK/bam/$S.main.bam
  echo "== $S: genes =="
  featureCounts -T $THREADS -p --countReadPairs -s $STRAND_N -a $REF/genes.exons.gtf \
    -o $WORK/counts/$S.genes.txt $B > /dev/null 2>&1
  for LEVEL in subfamily locus; do
    G=gene_id; [ $LEVEL = locus ] && G=transcript_id
    echo "== $S: TE $LEVEL =="
    featureCounts -T $THREADS -p --countReadPairs -s $STRAND_N -a $REF/te.gtf -g $G -O --fraction \
      -o $WORK/counts/$S.te_${LEVEL}_uniq.txt $B > /dev/null 2>&1
    featureCounts -T $THREADS -p --countReadPairs -s $STRAND_N -a $REF/te.gtf -g $G -O --fraction -M \
      -o $WORK/counts/$S.te_${LEVEL}_multi.txt $B > /dev/null 2>&1
  done
done
cp $WORK/counts/*.summary $RESULTS/counts/
