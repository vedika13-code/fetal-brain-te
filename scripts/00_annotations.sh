#!/usr/bin/env bash
# Step 0: gene annotation (GENCODE) and a TE annotation built from UCSC RepeatMasker (hg38)
set -euo pipefail
source config.sh
mkdir -p $REF $RESULTS

# 1. Genes: GENCODE v44, main chromosomes only
curl -fL https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz \
  | gunzip -c | awk '/^#/ || $1 ~ /^chr([0-9]+|X|Y)$/' > $REF/genes.gtf

# 2. TEs: RepeatMasker table -> GTF. One line per TE copy.
#    rmsk columns: 3=milliDiv  6=chrom 7=start 8=end 10=strand 11=repName 12=repClass 13=repFamily
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz | gunzip -c \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      $6 ~ /^chr([0-9]+|X|Y)$/ && $12 ~ /^(LINE|SINE|LTR|DNA|Retroposon)$/ {
        n[$11]++
        attr = "gene_id \"" $11 "\"; transcript_id \"" $11 "_dup" n[$11] "\"; family_id \"" $13 "\"; class_id \"" $12 "\"; milliDiv \"" $3 "\";"
        print $6, "rmsk", "exon", $7+1, $8, ".", $10, ".", attr
      }' > $REF/te.gtf

# 3. Summary of what's in it
{
  echo -e "class\tcopies"
  awk -F'class_id "' '{split($2,a,"\""); print a[1]}' $REF/te.gtf | sort | uniq -c | awk '{print $2"\t"$1}'
  echo -e "total\t$(wc -l < $REF/te.gtf)"
  echo -e "genes_gtf_lines\t$(wc -l < $REF/genes.gtf)"
} > $RESULTS/00_annotation_summary.tsv
cat $RESULTS/00_annotation_summary.tsv
echo "Example lines:"; grep 'gene_id "L1HS"' $REF/te.gtf | head -2
