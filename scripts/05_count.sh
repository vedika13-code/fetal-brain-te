#!/usr/bin/env bash
# Step 5: strandedness, gene counts, and TE counts with featureCounts
# TE counts are made 4 ways: {subfamily, locus} x {unique reads only, multimappers split fractionally}
set -euo pipefail
source config.sh
mkdir -p $WORK/counts $RESULTS/counts
STRAND_N=2    # reverse-stranded library (see 05_strandedness.tsv): featureCounts -s 2

for S in $SAMPLES; do
  B=$WORK/bam/$S.main.bam
  echo "== $S: genes =="
  featureCounts -T $THREADS -p --countReadPairs -s $STRAND_N -a $REF/genes.exons.gtf \
    -o $WORK/counts/$S.genes.txt $B > /dev/null 2>&1

  for LEVEL in subfamily locus; do
    G=gene_id; [ $LEVEL = locus ] && G=transcript_id
    echo "== $S: TE $LEVEL, unique only =="
    featureCounts -T $THREADS -p --countReadPairs -s $STRAND_N -a $REF/te.gtf -g $G -O \
      -o $WORK/counts/$S.te_${LEVEL}_uniq.txt $B > /dev/null 2>&1
    echo "== $S: TE $LEVEL, multimappers split fractionally =="
    featureCounts -T $THREADS -p --countReadPairs -s $STRAND_N -a $REF/te.gtf -g $G -O \
      -M --fraction \
      -o $WORK/counts/$S.te_${LEVEL}_multi.txt $B > /dev/null 2>&1
  done
done

# summary: how many read pairs each run assigned
echo -e "sample\tgenes\tte_subfamily_uniq\tte_subfamily_multi" > $RESULTS/counts/05_assigned.tsv
for S in $SAMPLES; do
  a() { awk '$1=="Assigned"{print $2}' $WORK/counts/$S.$1.txt.summary; }
  echo -e "$S\t$(a genes)\t$(a te_subfamily_uniq)\t$(a te_subfamily_multi)" >> $RESULTS/counts/05_assigned.tsv
done
cp $WORK/counts/*.summary $RESULTS/counts/
cat $RESULTS/counts/05_assigned.tsv
