#!/usr/bin/env bash
# Step 4: remove alignments to _alt/_random/chrUn contigs and recompute NH per read
set -euo pipefail
source config.sh
mkdir -p $RESULTS/align
head -1 > $RESULTS/align/04_alt_effect.tsv < /dev/null || true
first=1
for S in $SAMPLES; do
  echo "== $S =="
  python scripts/04_primary_only.py $WORK/bam/$S.bam $WORK/bam/$S.main.bam $WORK/bam/$S.alt_stats.tsv $S
  if [ $first = 1 ]; then cat $WORK/bam/$S.alt_stats.tsv > $RESULTS/align/04_alt_effect.tsv; first=0
  else tail -1 $WORK/bam/$S.alt_stats.tsv >> $RESULTS/align/04_alt_effect.tsv; fi
done
cat $RESULTS/align/04_alt_effect.tsv
