#!/usr/bin/env python3
"""
Keep only alignments on the main chromosomes (chr1-22, X, Y, M) and fix the NH tag.

Why: the UCSC hg38 index contains alternate-haplotype contigs (chr*_alt). A read from a
region that exists both on a main chromosome and on an _alt contig gets two equally good
alignments, so it *looks* like a multimapper even though it comes from one locus. Our TE
annotation covers only the main chromosomes, so we drop _alt/_random/chrUn alignments and
recompute NH (number of alignments per read) from what is left.

Input must be grouped by read name (HISAT2's unsorted output is).
Usage: primary_only.py in.bam out.bam stats.tsv sample_name
"""
import re
import sys

import pysam

MAIN = re.compile(r"^chr([0-9]+|X|Y|M)$")


def fix_mate(records):
    """records: alignments of ONE mate on main chromosomes. Set NH, keep exactly one primary."""
    n = len(records)
    for r in records:
        r.set_tag("NH", n)
    if n and all(r.is_secondary for r in records):
        records[0].flag &= ~0x100  # promote one alignment to primary
    return records


def main(inp, out, stats_path, sample):
    bam = pysam.AlignmentFile(inp, "rb")
    keep = {i for i, name in enumerate(bam.references) if MAIN.match(name)}
    outb = pysam.AlignmentFile(out, "wb", template=bam)
    s = dict(reads=0, multi_before=0, multi_after=0, became_unique=0, lost_all=0)

    def flush(group):
        if not group:
            return
        s["reads"] += 1
        was_multi = any(r.get_tag("NH") > 1 for r in group if r.has_tag("NH"))
        kept = [r for r in group if not r.is_unmapped and r.reference_id in keep]
        mates = {}
        for r in kept:
            mates.setdefault(1 if r.is_read1 else 2 if r.is_read2 else 0, []).append(r)
        for recs in mates.values():
            fix_mate(recs)
        is_multi = any(len(v) > 1 for v in mates.values())
        s["multi_before"] += was_multi
        s["multi_after"] += is_multi
        if was_multi and kept and not is_multi:
            s["became_unique"] += 1
        if not kept:
            s["lost_all"] += 1
        for r in kept:
            outb.write(r)

    group, name = [], None
    for r in bam.fetch(until_eof=True):
        if r.query_name != name:
            flush(group)
            group, name = [], r.query_name
        group.append(r)
    flush(group)
    outb.close()

    with open(stats_path, "w") as fh:
        fh.write("sample\tread_pairs\tmultimapped_before\tmultimapped_after\tbecame_unique\tonly_on_alt_contigs\n")
        fh.write(f"{sample}\t{s['reads']}\t{s['multi_before']}\t{s['multi_after']}\t{s['became_unique']}\t{s['lost_all']}\n")
    print(open(stats_path).read().strip())


if __name__ == "__main__":
    main(*sys.argv[1:5])
