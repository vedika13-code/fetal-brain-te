#!/usr/bin/env python3
"""
Step 6 analysis: turn featureCounts TE tables into figures and tables.

Questions:
  1. Which TE classes do the TE reads come from?                        (fig1)
  2. What does ignoring multimappers cost, and for which subfamilies?   (fig2)
     Young subfamilies (low divergence from consensus) should lose most.
  3. Where do expressed TE loci sit relative to genes, and on which strand? (fig3)
     Stranded data: a TE read is counted only on the TE's own strand.
       - inside a gene, SAME strand as the gene  -> could be the gene's pre-mRNA (read-through)
       - inside a gene, OPPOSITE strand          -> cannot be the gene's own transcript
       - intergenic                              -> no host gene to explain it
  4. (optional) Expressed TE loci inside SFARI autism genes.

n = 1 sample per stage: descriptive only, no statistical tests.
Usage: 06_analysis.py <samples space-separated> <counts_dir> <analysis_dir> <ref_dir>
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
GREY, INK = "#8a8a85", "#2b2b29"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e5e0", "axes.axisbelow": True,
                     "legend.frameon": False, "figure.dpi": 150})
CLASS_COL = {"LINE": C[0], "SINE": C[1], "LTR": C[2], "DNA": C[3], "Retroposon": C[4]}
LABEL = ["L1HS", "L1PA2", "L1PA3", "AluYa5", "AluYb8", "AluY", "SVA_D", "SVA_F", "HERVK-int", "LTR5_Hs",
         "L1MC4", "AluSx", "MIR", "L2a", "HERVH-int"]


def read_fc(path):
    """featureCounts table -> Series(id -> count). Reads only the id and count columns."""
    df = pd.read_csv(path, sep="\t", comment="#", usecols=[0, 6])
    df.columns = ["id", "count"]
    return df.set_index("id")["count"].astype(float)


def assigned(summary):
    s = pd.read_csv(summary, sep="\t", index_col=0)
    return float(s.loc["Assigned"].iloc[0])


def main(samples, cdir, adir, rdir):
    os.makedirs(adir, exist_ok=True)
    info = pd.read_csv(f"{adir}/subfamily_info.tsv", sep="\t", index_col=0)  # subfamily -> family, class, mean_div, copies

    # ---------------- subfamily level ----------------
    sub = info.copy()
    lib = {}
    for s in samples:
        lib[s] = assigned(f"{cdir}/{s}.genes.txt.summary")
        sub[f"uniq_{s}"] = read_fc(f"{cdir}/{s}.te_subfamily_uniq.txt").reindex(sub.index).fillna(0)
        sub[f"multi_{s}"] = read_fc(f"{cdir}/{s}.te_subfamily_multi.txt").reindex(sub.index).fillna(0)
        sub[f"cpm_{s}"] = sub[f"multi_{s}"] / lib[s] * 1e6
    sub["multi_total"] = sum(sub[f"multi_{s}"] for s in samples)
    sub["uniq_total"] = sum(sub[f"uniq_{s}"] for s in samples)
    sub["pct_recovered_by_unique"] = 100 * sub["uniq_total"] / sub["multi_total"].replace(0, np.nan)
    sub.sort_values("multi_total", ascending=False).round(2).to_csv(f"{adir}/subfamily_counts.tsv", sep="\t")

    # fig1: class composition
    cls = pd.DataFrame({s: sub.groupby("class")[f"multi_{s}"].sum() for s in samples})
    cls = 100 * cls / cls.sum()
    fig, ax = plt.subplots(figsize=(6.5, 3))
    x = np.arange(len(cls)); w = 0.8 / len(samples)
    for i, s in enumerate(samples):
        b = ax.bar(x + i * w - 0.4 + w / 2, cls[s], w - 0.03, color=[C[0], C[1]][i % 2], label=s)
        ax.bar_label(b, fmt="%.0f", fontsize=7, padding=1, color=INK)
    ax.set_xticks(x, cls.index); ax.set_ylabel("% of TE-assigned reads")
    ax.set_title("1. TE reads by class (multimappers included)", loc="left", fontsize=10)
    ax.legend(); ax.grid(axis="x", visible=False)
    fig.tight_layout(); fig.savefig(f"{adir}/fig1_te_classes.png"); plt.close(fig)

    # fig2: % of signal recovered by unique-only counting vs divergence (age proxy)
    d = sub[sub["multi_total"] >= 50].copy()
    d["div_pct"] = d["mean_div"] / 10.0
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for c, col in CLASS_COL.items():
        dd = d[d["class"] == c]
        ax.scatter(dd["div_pct"], dd["pct_recovered_by_unique"], s=12 + 30 * np.log10(dd["multi_total"]),
                   color=col, alpha=0.75, label=c, edgecolor="white", linewidth=0.6)
    for name in LABEL:
        if name in d.index:
            ax.annotate(name, (d.loc[name, "div_pct"], d.loc[name, "pct_recovered_by_unique"]),
                        fontsize=7, color=INK, xytext=(4, 2), textcoords="offset points")
    ax.set_xlabel("mean divergence from subfamily consensus (%)  ->  older")
    ax.set_ylabel("% of signal kept by unique-only counting")
    ax.set_ylim(0, 105)
    ax.set_title("2. Young TEs lose the most reads when multimappers are discarded", loc="left", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(f"{adir}/fig2_unique_vs_age.png"); plt.close(fig)

    # ---------------- locus level ----------------
    ctx = pd.read_csv(f"{adir}/locus_context.tsv", sep="\t", index_col=0)  # locus -> context, host_genes
    loc = pd.DataFrame({s: read_fc(f"{cdir}/{s}.te_locus_multi.txt").reindex(ctx.index).fillna(0) for s in samples})
    locu = pd.DataFrame({s: read_fc(f"{cdir}/{s}.te_locus_uniq.txt").reindex(ctx.index).fillna(0) for s in samples})
    loc = loc.join(locu, rsuffix="_uniq").join(ctx)
    loc.round(2).to_csv(f"{adir}/expressed_loci.tsv", sep="\t")

    order = ["exonic_sense", "exonic_antisense", "intronic_sense", "intronic_antisense", "intergenic"]
    cnt = pd.DataFrame({s: loc.loc[loc[s] >= 5, "context"].value_counts() for s in samples}).reindex(order).fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    x = np.arange(len(cnt)); w = 0.8 / len(samples)
    for i, s in enumerate(samples):
        b = ax.bar(x + i * w - 0.4 + w / 2, cnt[s], w - 0.03, color=[C[0], C[1]][i % 2], label=s)
        ax.bar_label(b, fontsize=7, padding=1, color=INK)
    ax.set_xticks(x, [o.replace("_", "\n") for o in order]); ax.set_ylabel("TE loci with >=5 reads")
    ax.set_title("3. Where expressed TE loci sit, relative to genes and their strand", loc="left", fontsize=10)
    ax.legend(); ax.grid(axis="x", visible=False)
    fig.tight_layout(); fig.savefig(f"{adir}/fig3_locus_context.png"); plt.close(fig)

    # ---------------- SFARI (optional) ----------------
    sf_path = f"{rdir}/sfari_genes.txt"
    sf_lines = []
    if os.path.exists(sf_path):
        sfari = {g.strip() for g in open(sf_path) if g.strip()}
        hosts = loc["host_genes"].fillna("").str.split(",")
        in_sf = hosts.apply(lambda gs: sorted(set(gs) & sfari))
        sfl = loc[in_sf.str.len() > 0].copy()
        sfl["sfari_genes"] = in_sf[in_sf.str.len() > 0].str.join(",")
        sfl = sfl[(sfl[samples] >= 5).any(axis=1)]
        sfl.sort_values(samples[-1], ascending=False).round(2).to_csv(f"{adir}/sfari_expressed_loci.tsv", sep="\t")
        n_genes = len(set(",".join(sfl["sfari_genes"]).split(","))) if len(sfl) else 0
        sf_lines = [f"- Expressed TE loci (>=5 reads in a sample) inside SFARI genes: {len(sfl)}, in {n_genes} genes",
                    "- By context: " + ", ".join(f"{k}={v}" for k, v in sfl["context"].value_counts().items()),
                    "- Top genes: " + ", ".join(sfl["sfari_genes"].str.split(",").explode().value_counts().head(10).index)]

    # ---------------- summary ----------------
    with open(f"{adir}/summary.md", "w") as fh:
        fh.write("# TE expression in prenatal brain: summary\n\nn = 1 per stage; descriptive only.\n\n")
        fh.write("| sample | gene read pairs | TE (unique only) | TE (with multimappers) | gain from multimappers |\n|---|---|---|---|---|\n")
        for s in samples:
            u, m = sub[f"uniq_{s}"].sum(), sub[f"multi_{s}"].sum()
            fh.write(f"| {s} | {lib[s]:,.0f} | {u:,.0f} | {m:,.0f} | +{100*(m-u)/u:.0f}% |\n")
        fh.write("\n## TE reads by class (% of TE reads, multimappers included)\n\n" + cls.round(1).to_markdown() + "\n")
        young = sub.loc[[n for n in ["L1HS", "L1PA2", "AluYa5", "AluYb8", "SVA_F"] if n in sub.index],
                        ["class", "mean_div", "uniq_total", "multi_total", "pct_recovered_by_unique"]]
        fh.write("\n## Young subfamilies: how much unique-only counting keeps\n\n" + young.round(1).to_markdown() + "\n")
        fh.write("\n## Expressed TE loci (>=5 reads) by genomic context\n\n" + cnt.to_markdown() + "\n")
        top = sub.sort_values(f"cpm_{samples[0]}", ascending=False).head(15)[["class"] + [f"cpm_{s}" for s in samples]]
        fh.write("\n## Top 15 subfamilies (CPM, relative to gene-assigned read pairs)\n\n" + top.round(1).to_markdown() + "\n")
        if sf_lines:
            fh.write("\n## SFARI autism genes\n\n" + "\n".join(sf_lines) + "\n")
    print(open(f"{adir}/summary.md").read())


if __name__ == "__main__":
    main(sys.argv[1].split(), sys.argv[2], sys.argv[3], sys.argv[4])
