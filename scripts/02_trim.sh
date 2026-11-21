#!/usr/bin/env bash
# Step 2: adapter + quality trimming and QC in one pass with fastp
set -euo pipefail
source config.sh
mkdir -p $WORK/trimmed $RESULTS/qc

for S in $SAMPLES; do
  echo "== fastp on $S =="
  fastp -i $DATA/fastq/${S}_R1.fq.gz -I $DATA/fastq/${S}_R2.fq.gz \
        -o $WORK/trimmed/${S}_R1.fq.gz -O $WORK/trimmed/${S}_R2.fq.gz \
        --detect_adapter_for_pe --cut_right --cut_right_mean_quality 20 --length_required 36 \
        -w $THREADS -j $RESULTS/qc/${S}.fastp.json -h $RESULTS/qc/${S}.fastp.html 2> $WORK/trimmed/${S}.fastp.log
done

# summary table from fastp's JSON reports
python - <<'PY'
import json
samples = open("config.sh").read().split('SAMPLES="')[1].split('"')[0].split()
print("sample\treads_before\treads_after\tpct_kept\tQ30_after\tGC\tduplication\tadapter_trimmed_reads")
rows = []
for s in samples:
    d = json.load(open(f"results/qc/{s}.fastp.json"))
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    rows.append(f'{s}\t{b["total_reads"]}\t{a["total_reads"]}\t{100*a["total_reads"]/b["total_reads"]:.1f}%\t'
                f'{100*a["q30_rate"]:.1f}%\t{100*a["gc_content"]:.1f}%\t{100*d["duplication"]["rate"]:.1f}%\t'
                f'{d["adapter_cutting"]["adapter_trimmed_reads"]}')
print("\n".join(rows))
open("results/qc/02_fastp_summary.tsv", "w").write(
    "sample\treads_before\treads_after\tpct_kept\tQ30_after\tGC\tduplication\tadapter_trimmed_reads\n" + "\n".join(rows) + "\n")
PY
