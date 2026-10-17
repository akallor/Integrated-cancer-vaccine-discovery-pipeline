#!/usr/bin/env python3
"""07_id_diagnostics.py -- Module 2 identification diagnostics for one HLA class.

Compares Sage (pre-rescoring) vs MS2Rescore (post-rescoring) PSM tables and writes a 3-panel figure:
  (a) length distribution of the rescored UNMODIFIED peptides (canonical vs non-canonical, stacked)
  (b) #PSMs and #peptides identified, pre vs post rescoring (grouped bars)
  (c) canonical vs non-canonical peptides, pre vs post rescoring (grouped bars)

Both tables are filtered to TARGET PSMs at the chosen FDR (default 1%). Canonical vs non-canonical
comes from the |grp=CANON|/|grp=NONCANON| tag carried in the 'proteins' column (from BUILD_FASTA).
Column names are auto-detected so it works with either Sage or mokapot output layouts.
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRE_C, POST_C = "#7f8c8d", "#2c7fb8"          # pre / post
CANON_C, NONCANON_C = "#95a5a6", "#e0821e"    # canonical / non-canonical


def _find(df, *names, contains=None):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    if contains:
        for c in df.columns:
            if contains in c.lower():
                return c
    return None


def bare(pep):
    """Strip flanks (K.PEP.R), modifications ([..]/(..)), charge (/2) -> plain residues."""
    s = str(pep)
    if s.count(".") >= 2:
        parts = s.split(".")
        if len(parts) > 2:
            s = ".".join(parts[1:-1])
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r"\([^)]*\)", "", s)
    s = re.sub(r"/\d+$", "", s)
    return re.sub(r"[^A-Za-z]", "", s).upper()


def is_modified(pep):
    s = str(pep)
    return bool(re.search(r"[\[\(]", s) or re.search(r"[a-z]", s))


def category(proteins):
    p = str(proteins).upper()
    if "GRP=CANON" in p:
        return "canonical"
    if "GRP=NONCANON" in p:
        return "non-canonical"
    return "unknown"


def load_psms(path, kind, fdr):
    d = pd.read_csv(path, sep="\t", dtype=str)
    pep = _find(d, "peptide")
    prot = _find(d, "proteins", "protein", contains="protein")   # matches proteins / protein_list
    if pep is None:
        sys.exit(f"{path}: no 'peptide' column (cols={list(d.columns)})")
    if prot is None:
        looks_recovery = any(_find(d, c) for c in ("recovered", "searchable", "peptide_bare", "note", "class"))
        want = "Sage results.sage.tsv" if kind == "sage" else "MS2Rescore rescored_psms.tsv"
        hint = ("  This looks like the group_fdr recovery list, not the rescored PSM table."
                if looks_recovery else "")
        sys.exit(f"{path}: no 'proteins' column. The --{kind} input must be a PSM table with a "
                 f"'proteins' column carrying the |grp=CANON/NONCANON| tag (i.e. {want})." + hint)
    q = (_find(d, "spectrum_q", "psm_q") if kind == "sage"
         else (_find(d, contains="q-value") or _find(d, contains="q_value") or _find(d, contains="qvalue")))
    if q is None:
        q = _find(d, contains="_q")
    if q is None:
        print(f"[WARN] {path}: no q-value column found; keeping all target rows (no FDR filter).", file=sys.stderr)
    tgt, decoy, lab = _find(d, "is_target"), _find(d, "is_decoy"), _find(d, "label")

    d = d.rename(columns={pep: "_pep"})
    d["_prot"] = d[prot] if prot else ""
    d["_q"] = pd.to_numeric(d[q], errors="coerce") if q else 0.0
    if tgt:
        d["_target"] = d[tgt].astype(str).str.lower().isin(["true", "1", "target"])
    elif decoy:
        d["_target"] = ~d[decoy].astype(str).str.lower().isin(["true", "1", "-1", "decoy"])
    elif lab:
        lv = pd.to_numeric(d[lab], errors="coerce")
        d["_target"] = (lv > 0) if lv.notna().any() else ~d[lab].astype(str).str.lower().isin(["false", "decoy", "-1"])
    else:
        d["_target"] = True
    d = d[d["_target"] & (d["_q"] <= fdr)].copy()
    d["_bare"] = d["_pep"].map(bare)
    d["_mod"] = d["_pep"].map(is_modified)
    d["_cat"] = d["_prot"].map(category)
    d["_len"] = d["_bare"].str.len()
    return d


def counts(d):
    return dict(psms=len(d), peptides=int(d["_bare"].nunique()),
                canonical=int(d[d["_cat"] == "canonical"]["_bare"].nunique()),
                noncanonical=int(d[d["_cat"] == "non-canonical"]["_bare"].nunique()))


def _pair_bars(ax, labels, pre, post, ylab, title):
    x = np.arange(len(labels)); w = 0.38
    ax.bar(x - w / 2, pre, w, color=PRE_C, label="pre (Sage)")
    ax.bar(x + w / 2, post, w, color=POST_C, label="post (MS\u00b2Rescore)")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(ylab); ax.set_title(title); ax.legend()
    for i, (a, b) in enumerate(zip(pre, post)):
        ax.text(i - w / 2, a, str(a), ha="center", va="bottom", fontsize=8)
        ax.text(i + w / 2, b, str(b), ha="center", va="bottom", fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sage", required=True, help="Sage results.sage.tsv (pre-rescoring)")
    ap.add_argument("--rescored", required=True, help="MS2Rescore rescored_psms.tsv (post-rescoring)")
    ap.add_argument("--hla-class", default="I")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--fdr", type=float, default=0.01)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    cls = a.hla_class

    pre = load_psms(a.sage, "sage", a.fdr)
    post = load_psms(a.rescored, "rescore", a.fdr)
    cpre, cpost = counts(pre), counts(post)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.4))

    # (a) length distribution — rescored UNMODIFIED peptides, canonical vs non-canonical stacked
    um = post[~post["_mod"]].drop_duplicates("_bare")
    lens = sorted(int(x) for x in um["_len"].unique())
    if lens:
        can = [int(((um["_cat"] == "canonical") & (um["_len"] == L)).sum()) for L in lens]
        non = [int(((um["_cat"] == "non-canonical") & (um["_len"] == L)).sum()) for L in lens]
        axes[0].bar(lens, can, color=CANON_C, label="canonical")
        axes[0].bar(lens, non, bottom=can, color=NONCANON_C, label="non-canonical")
        axes[0].set_xticks(lens)
        axes[0].legend()
    axes[0].set_xlabel("peptide length (aa)"); axes[0].set_ylabel("unique peptides")
    axes[0].set_title(f"(a) Length \u2014 rescored unmodified peptides\nHLA-{cls}  (n={len(um)})")

    # (b) PSMs & peptides pre vs post
    _pair_bars(axes[1], ["PSMs", "peptides"],
               [cpre["psms"], cpre["peptides"]], [cpost["psms"], cpost["peptides"]],
               "count", f"(b) IDs pre vs post rescoring\nHLA-{cls}  ({a.fdr*100:.0f}% FDR)")

    # (c) canonical vs non-canonical pre vs post
    _pair_bars(axes[2], ["canonical", "non-canonical"],
               [cpre["canonical"], cpre["noncanonical"]], [cpost["canonical"], cpost["noncanonical"]],
               "unique peptides", f"(c) Canonical vs non-canonical\nHLA-{cls}  (pre vs post)")

    fig.suptitle(f"Module 2 \u2014 identification diagnostics (HLA-{cls})", fontweight="bold", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = os.path.join(a.outdir, f"id_diagnostics_HLA-{cls}.png")
    fig.savefig(out, dpi=150); plt.close(fig)

    pd.DataFrame([dict(stage="pre",  hla_class=cls, **cpre),
                  dict(stage="post", hla_class=cls, **cpost)]).to_csv(
        os.path.join(a.outdir, f"id_diagnostics_HLA-{cls}.tsv"), sep="\t", index=False)
    print(f"[id-diag] HLA-{cls}  pre={cpre}  post={cpost}")
    gain_p = 100 * (cpost['peptides'] - cpre['peptides']) / cpre['peptides'] if cpre['peptides'] else 0
    print(f"[id-diag] rescoring peptide gain: {gain_p:+.1f}%")
    print(f"[id-diag] wrote {out} and id_diagnostics_HLA-{cls}.tsv")


if __name__ == "__main__":
    main()
