#!/usr/bin/env python3
"""06_consensus.py -- consensus of the MHCflurry and MHCnuggets HLA-I arms.

Run ONCE both predictor tables exist. Merges them per peptide, tags each peptide by
agreement/confidence, and writes tables + figures.

Inputs (each: peptide, sample, binding_allele, <predictor>_binder_call [SB|WB|NB]):
  --flurry <mhcflurry pred tsv>   --nuggets <mhcnuggets pred tsv>

Outputs (--outdir):
  consensus.tsv                     per peptide: both calls, both alleles, agree, confidence tier
  consensus_by_allele.tsv           per-allele agreement stats
  consensus_concordance.png         binder-overlap Venn + call confusion heatmap + per-allele
                                    concordance heatmap, with % agreement/disagreement + footnote
  consensus_counts_by_allele.png    grouped bars: #SB/#WB/#NB per allele per predictor (colour=predictor)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLASSES = ["SB", "WB", "NB"]
FLURRY_C, NUGGETS_C = "#1f77b4", "#e0821e"   # predictor colours (used everywhere)


def _find(df, *names, contains=None):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n in low:
            return low[n]
    if contains:
        for c in df.columns:
            if contains in c.lower():
                return c
    return None


def _load(path, tag):
    d = pd.read_csv(path, sep="\t", dtype=str)
    pep = _find(d, "peptide")
    call = _find(d, contains="binder_call")
    allele = _find(d, "binding_allele")
    samp = _find(d, "sample")
    if not (pep and call and allele):
        sys.exit(f"{path}: need peptide + *_binder_call + binding_allele; got {list(d.columns)}")
    out = pd.DataFrame({"peptide": d[pep].astype(str),
                        f"{tag}_call": d[call].astype(str).str.upper().str.strip(),
                        f"{tag}_allele": d[allele].astype(str)})
    out["sample"] = d[samp].astype(str) if samp else "NA"
    return out.dropna(subset=["peptide"]).drop_duplicates(["peptide", "sample"])


def _confidence(fc, nc):
    fb, nb = fc in ("SB", "WB"), nc in ("SB", "WB")
    return "both_binder" if fb and nb else ("flurry_only" if fb else ("nuggets_only" if nb else "neither"))


def draw_venn(ax, f_only, both, n_only, neither):
    """Dependency-free 2-set Venn of binder (SB/WB) peptides."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(plt.Circle((0.40, 0.52), 0.30, color=FLURRY_C, alpha=0.45))
    ax.add_patch(plt.Circle((0.62, 0.52), 0.30, color=NUGGETS_C, alpha=0.45))
    ax.text(0.26, 0.52, f"{f_only}", ha="center", va="center", fontsize=15, fontweight="bold")
    ax.text(0.51, 0.52, f"{both}", ha="center", va="center", fontsize=16, fontweight="bold", color="#222")
    ax.text(0.76, 0.52, f"{n_only}", ha="center", va="center", fontsize=15, fontweight="bold")
    ax.text(0.40, 0.88, "MHCflurry\nbinders", ha="center", va="center", fontsize=10, color=FLURRY_C, fontweight="bold")
    ax.text(0.62, 0.16, "MHCnuggets\nbinders", ha="center", va="center", fontsize=10, color=NUGGETS_C, fontweight="bold")
    ax.text(0.51, 0.30, "both = high-confidence", ha="center", fontsize=8, style="italic", color="#333")
    ax.text(0.99, 0.02, f"neither (both NB): {neither}", ha="right", fontsize=9, color="#555")
    ax.set_title("Binder overlap (SB/WB)", fontsize=12)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flurry", required=True)
    ap.add_argument("--nuggets", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--allele-basis", choices=["flurry", "nuggets"], default="flurry",
                    help="which predictor's binding_allele groups the per-allele concordance stats")
    ap.add_argument("--min-allele-n", type=int, default=20,
                    help="min peptides for an allele to be ranked in the most/least-agreed footnote")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    f, n = _load(a.flurry, "flurry"), _load(a.nuggets, "nuggets")
    m = f.merge(n, on=["peptide", "sample"], how="inner")
    if m.empty:
        sys.exit("no shared (peptide, sample) between the two tables")
    m = m[m["flurry_call"].isin(CLASSES) & m["nuggets_call"].isin(CLASSES)].copy()
    m["agree"] = m["flurry_call"] == m["nuggets_call"]
    m["confidence"] = [_confidence(fc, nc) for fc, nc in zip(m["flurry_call"], m["nuggets_call"])]
    m["concordance"] = np.where(m["agree"], "both_" + m["flurry_call"],
                                m["flurry_call"] + "_vs_" + m["nuggets_call"])
    m["allele"] = m[f"{a.allele_basis}_allele"]

    N = len(m)
    agree_pct = 100.0 * m["agree"].mean()
    dis_pct = 100.0 - agree_pct
    conf = pd.crosstab(m["flurry_call"], m["nuggets_call"]).reindex(index=CLASSES, columns=CLASSES, fill_value=0)

    # per-allele concordance (grouped by the chosen allele basis)
    per = []
    for al, g in m.groupby("allele"):
        per.append({"allele": al, "n": len(g),
                    "both_SB": int(((g.flurry_call == "SB") & (g.nuggets_call == "SB")).sum()),
                    "both_WB": int(((g.flurry_call == "WB") & (g.nuggets_call == "WB")).sum()),
                    "both_NB": int(((g.flurry_call == "NB") & (g.nuggets_call == "NB")).sum()),
                    "disagree": int((~g.agree).sum()),
                    "agree_pct": round(100.0 * g.agree.mean(), 1)})
    per = pd.DataFrame(per).sort_values("n", ascending=False).reset_index(drop=True)

    elig = per[per.n >= a.min_allele_n]
    if len(elig):
        top = elig.loc[elig.agree_pct.idxmax()]
        bot = elig.loc[elig.agree_pct.idxmin()]
        foot = (f"Most concordant allele: {top.allele} ({top.agree_pct:.0f}% agree, n={int(top.n)})   |   "
                f"Least concordant: {bot.allele} ({bot.agree_pct:.0f}% agree, n={int(bot.n)})")
    else:
        foot = "(too few peptides per allele to rank agreement)"

    # confidence tiers + venn quantities
    tiers = m["confidence"].value_counts().to_dict()
    both_b = tiers.get("both_binder", 0)
    f_only = tiers.get("flurry_only", 0)
    n_only = tiers.get("nuggets_only", 0)
    neither = tiers.get("neither", 0)

    # ---- tables ----
    m.to_csv(os.path.join(a.outdir, "consensus.tsv"), sep="\t", index=False)
    per.to_csv(os.path.join(a.outdir, "consensus_by_allele.tsv"), sep="\t", index=False)

    # ================= FIG 1: concordance (Venn + confusion + per-allele) =================
    fig = plt.figure(figsize=(19, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.35], wspace=0.28)
    axV = fig.add_subplot(gs[0]); draw_venn(axV, f_only, both_b, n_only, neither)

    axC = fig.add_subplot(gs[1])
    im = axC.imshow(conf.values, cmap="Blues")
    axC.set_xticks(range(3)); axC.set_xticklabels(CLASSES)
    axC.set_yticks(range(3)); axC.set_yticklabels(CLASSES)
    axC.set_xlabel("MHCnuggets call"); axC.set_ylabel("MHCflurry call")
    vmax = conf.values.max() or 1
    for i in range(3):
        for j in range(3):
            v = int(conf.values[i, j])
            axC.text(j, i, str(v), ha="center", va="center",
                     color="white" if v > vmax * 0.5 else "black",
                     fontweight="bold" if i == j else "normal")
            if i == j:
                axC.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#27ae60", lw=3))
    axC.set_title(f"Call confusion  (n={N})\nagreement {agree_pct:.0f}%   ·   disagreement {dis_pct:.0f}%")

    axA = fig.add_subplot(gs[2])
    pa = per.set_index("allele")[["both_SB", "both_WB", "both_NB", "disagree"]]
    axA.imshow(pa.values, cmap="Blues", aspect="auto")
    axA.set_xticks(range(4)); axA.set_xticklabels(["both SB", "both WB", "both NB", "disagree"])
    axA.set_yticks(range(len(pa))); axA.set_yticklabels(pa.index)
    pamax = pa.values.max() or 1
    for i in range(len(pa)):
        for j in range(4):
            axA.text(j, i, str(int(pa.values[i, j])), ha="center", va="center",
                     color="white" if pa.values[i, j] > pamax * 0.5 else "black")
        axA.text(3.75, i, f"{per.set_index('allele').loc[pa.index[i], 'agree_pct']:.0f}%",
                 va="center", ha="left", fontsize=9, color="#27ae60", fontweight="bold")
    axA.set_xlim(-0.5, 4.4)
    axA.set_title("Per-allele concordance (peptide count)\n"
                  f"grouped by {a.allele_basis} best allele; green = % agree")

    fig.suptitle("MHCflurry vs MHCnuggets — consensus concordance", fontweight="bold", fontsize=15)
    fig.text(0.5, 0.015, foot, ha="center", fontsize=10.5, style="italic")
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    fig.savefig(os.path.join(a.outdir, "consensus_concordance.png"), dpi=150)
    plt.close(fig)

    # ================= FIG 2: grouped bars, #SB/#WB/#NB per allele per predictor =================
    fcnt = m.groupby(["flurry_allele", "flurry_call"]).size()
    ncnt = m.groupby(["nuggets_allele", "nuggets_call"]).size()
    alleles = sorted(set(m["flurry_allele"]) | set(m["nuggets_allele"]))
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), sharey=True)
    x = np.arange(len(alleles)); w = 0.38
    for ax, cls in zip(axes, CLASSES):
        fv = [int(fcnt.get((al, cls), 0)) for al in alleles]
        nv = [int(ncnt.get((al, cls), 0)) for al in alleles]
        ax.bar(x - w / 2, fv, w, color=FLURRY_C, label="MHCflurry")
        ax.bar(x + w / 2, nv, w, color=NUGGETS_C, label="MHCnuggets")
        ax.set_xticks(x); ax.set_xticklabels(alleles, rotation=45, ha="right", fontsize=8)
        ax.set_title({"SB": "Strong binders (SB)", "WB": "Weak binders (WB)", "NB": "Non-binders (NB)"}[cls])
        ax.set_xlabel("allele")
    axes[0].set_ylabel("number of unique peptides")
    axes[0].legend(title="Predictor", loc="upper left")
    fig.suptitle("Binder counts per allele per predictor", fontweight="bold", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(a.outdir, "consensus_counts_by_allele.png"), dpi=150)
    plt.close(fig)

    print(f"[consensus] n={N}  agreement={agree_pct:.1f}%  disagreement={dis_pct:.1f}%")
    print(f"[consensus] confidence tiers: both_binder={both_b}  flurry_only={f_only}  "
          f"nuggets_only={n_only}  neither={neither}")
    print(f"[consensus] {foot}")
    print("[consensus] wrote consensus.tsv, consensus_by_allele.tsv, "
          "consensus_concordance.png, consensus_counts_by_allele.png")


if __name__ == "__main__":
    main()
