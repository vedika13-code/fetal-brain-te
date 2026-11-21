#!/usr/bin/env bash
# Step 6: prepare annotation tables (subfamily age, locus context) and run the analysis
set -euo pipefail
source config.sh
C=$WORK/counts
A=$RESULTS/analysis
T=$WORK/analysis_tmp
mkdir -p $A $T

# helper: pull one attribute (e.g. gene_id) out of a GTF attribute column
ATTR='function attr(s, key) { if (match(s, key " \"[^\"]+\"")) return substr(s, RSTART+length(key)+2, RLENGTH-length(key)-3); return "" }'

# 1. per-subfamily: family, class, mean divergence from consensus (age proxy), number of copies
awk -F'\t' "$ATTR"'
  { g=attr($9,"gene_id"); n[g]++; d[g]+=attr($9,"milliDiv"); fam[g]=attr($9,"family_id"); cl[g]=attr($9,"class_id") }
  END { print "subfamily\tfamily\tclass\tmean_div\tcopies"; for (g in n) printf "%s\t%s\t%s\t%.1f\t%d\n", g, fam[g], cl[g], d[g]/n[g], n[g] }' \
  $REF/te.gtf > $A/subfamily_info.tsv

# 2. expressed loci: >=5 read pairs (multimappers included) in at least one sample
FILES=""; for S in $SAMPLES; do FILES="$FILES $C/$S.te_locus_multi.txt"; done
awk -F'\t' '$1!~/^#/ && $1!="Geneid" && $7+0>=5 {print $1}' $FILES | sort -u > $T/expressed_ids.txt
echo "Expressed TE loci (>=5 read pairs in a sample): $(wc -l < $T/expressed_ids.txt)"

# 3. BED files: expressed TE loci, gene bodies (with names), exons
awk -F'\t' "$ATTR"' NR==FNR {keep[$1]=1; next}
  { t=attr($9,"transcript_id"); if (t in keep) print $1"\t"$4-1"\t"$5"\t"t"\t.\t"$7 }' $T/expressed_ids.txt $REF/te.gtf > $T/loci.bed
awk -F'\t' "$ATTR"' $3=="gene" {print $1"\t"$4-1"\t"$5"\t"attr($9,"gene_name")"\t.\t"$7}' $REF/genes.gtf > $T/genes.bed
awk -F'\t' '$3=="exon" {print $1"\t"$4-1"\t"$5"\t.\t.\t"$7}' $REF/genes.gtf > $T/exons.bed

# 4. context: overlap with exons / gene bodies on the same (-s) or opposite (-S) strand
bedtools intersect -a $T/loci.bed -b $T/exons.bed -s -u | cut -f4 > $T/exonic_sense.txt
bedtools intersect -a $T/loci.bed -b $T/exons.bed -S -u | cut -f4 > $T/exonic_antisense.txt
bedtools intersect -a $T/loci.bed -b $T/genes.bed -s -u | cut -f4 > $T/genic_sense.txt
bedtools intersect -a $T/loci.bed -b $T/genes.bed -S -u | cut -f4 > $T/genic_antisense.txt
bedtools intersect -a $T/loci.bed -b $T/genes.bed -wa -wb | awk '{print $4"\t"$10}' | sort -u \
  | awk '{h[$1] = ($1 in h) ? h[$1]","$2 : $2} END {for (k in h) print k"\t"h[k]}' > $T/hosts.txt

# priority: exonic sense > exonic antisense > intronic sense > intronic antisense > intergenic
awk -F'\t' -v T=$T '
  BEGIN { while ((getline l < (T"/exonic_sense.txt")) > 0) es[l]=1; while ((getline l < (T"/exonic_antisense.txt")) > 0) ea[l]=1
          while ((getline l < (T"/genic_sense.txt")) > 0) gs[l]=1;   while ((getline l < (T"/genic_antisense.txt")) > 0) ga[l]=1
          while ((getline l < (T"/hosts.txt")) > 0) { split(l, a, "\t"); h[a[1]]=a[2] }
          print "locus\tcontext\thost_genes" }
  { id=$4; c = (id in es) ? "exonic_sense" : (id in ea) ? "exonic_antisense" : (id in gs) ? "intronic_sense" : (id in ga) ? "intronic_antisense" : "intergenic"
    print id"\t"c"\t"h[id] }' $T/loci.bed > $A/locus_context.tsv

# 5. figures + summary
python scripts/06_analysis.py "$SAMPLES" $C $A $REF
