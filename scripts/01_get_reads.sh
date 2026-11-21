#!/usr/bin/env bash
# Step 1: stream the first N reads of each sample from ENA (no full download)
set -euo pipefail
source config.sh
mkdir -p $DATA/fastq $RESULTS
echo -e "sample\tena_sample\trun\tlayout\treads_in_run\treads_kept" > $RESULTS/01_samples.tsv

for S in $SAMPLES; do
  ACC=${ENA[$S]}
  # ask ENA which sequencing run belongs to this sample, and where its FASTQs are
  info=$(curl -fsL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$ACC&result=read_run&fields=run_accession,library_layout,read_count,fastq_ftp&format=tsv" | sed -n 2p)
  run=$(echo "$info" | cut -f1); layout=$(echo "$info" | cut -f2); total=$(echo "$info" | cut -f3)
  urls=$(echo "$info" | cut -f4 | tr ';' '\n')
  echo "== $S ($ACC): run $run, $layout, $total reads in full run =="

  set +o pipefail   # 'head' stops the download early on purpose; don't treat that as an error
  if [ "$layout" = "PAIRED" ]; then
    u1=$(echo "$urls" | grep '_1.fastq.gz'); u2=$(echo "$urls" | grep '_2.fastq.gz')
    curl -sL "https://$u1" | gunzip -c 2>/dev/null | head -n $((4*MAX_READS)) | gzip > $DATA/fastq/${S}_R1.fq.gz
    curl -sL "https://$u2" | gunzip -c 2>/dev/null | head -n $((4*MAX_READS)) | gzip > $DATA/fastq/${S}_R2.fq.gz
  else
    u1=$(echo "$urls" | head -1)
    curl -sL "https://$u1" | gunzip -c 2>/dev/null | head -n $((4*MAX_READS)) | gzip > $DATA/fastq/${S}_R1.fq.gz
  fi
  set -o pipefail

  kept=$(( $(zcat $DATA/fastq/${S}_R1.fq.gz | wc -l) / 4 ))
  echo -e "$S\t$ACC\t$run\t$layout\t$total\t$kept" >> $RESULTS/01_samples.tsv
done
cat $RESULTS/01_samples.tsv
echo "Read length (first read):"; zcat $DATA/fastq/early_CS20_R1.fq.gz | sed -n 2p | awk '{print length($0)}'
