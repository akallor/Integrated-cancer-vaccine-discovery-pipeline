#!/usr/bin/env python3
"""
run_explain.py -- CLI driver + plotting for the MHCflurry explainability layer.

Example (EC2, real model):
  python run_explain.py --input results/mhcflurry_hlaI.tsv --outdir explain_out \
                        --backend mhcflurry --seed 7

Example (validation, no MHCflurry needed):
  python run_explain.py --input results/mhcflurry_hlaI.tsv --outdir explain_out \
                        --backend surrogate --seed 7
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

import mhcflurry_explain as mx

STRATA = [("CANON", "SB"), ("CANON", "WB"), ("NONCANON", "SB"), ("NONCANON", "WB")]
STRAT_LABEL = {("CANON", "SB"): "Canonical - strong",
               ("CANON", "WB"): "Canonical - weak",
               ("NONCANON", "SB"): "Non-canonical - strong",
               ("NONCANON", "WB"): "Non-canonical - weak"}

# muted, print-friendly palette
POS_C, NEG_C = "#2166ac", "#b2182b"


def _pep_positions(pep):
    """Position labels P1..P(L-1),POmega with residue letters."""
    L = len(pep)
    labs = []
    for i, a in enumerate(pep):
        if i == L - 1:
            labs.append(f"P\u03A9\n{a}")
        else:
            labs.append(f"P{i+1}\n{a}")
    return labs


# --------------------------------------------------------------------------------------
def plot_shap_grid(records, outpath):
    """4 strata (rows) x 2 peptides (cols) grid of positional-SHAP bar charts."""
    fig, axes = plt.subplots(4, 2, figsize=(11, 13))
    fig.suptitle("Positional SHAP — residue attribution to MHCflurry binding score\n"
                 "(blue = pushes toward strong binding, red = toward weak)",
                 fontsize=12, fontweight="bold", y=0.995)
    by_stratum = {s: [] for s in STRATA}
    for r in records:
        by_stratum[(r["group"], r["call"])].append(r)

    for row_i, stratum in enumerate(STRATA):
        recs = by_stratum[stratum]
        for col_i in range(2):
            ax = axes[row_i, col_i]
            if col_i >= len(recs):
                ax.axis("off"); continue
            r = recs[col_i]
            sv, pep = r["shap"], r["peptide"]
            colors = [POS_C if v >= 0 else NEG_C for v in sv]
            ax.bar(range(len(sv)), sv, color=colors, edgecolor="black", linewidth=0.4)
            ax.axhline(0, color="black", linewidth=0.6)
            ax.set_xticks(range(len(sv)))
            ax.set_xticklabels(_pep_positions(pep), fontsize=7)
            ax.set_title(f"{STRAT_LABEL[stratum]}  |  {pep}  |  {r['allele']}\n"
                         f"aff={r['affinity']:.0f} nM  score={r['wt_score']:.2f}",
                         fontsize=8)
            ax.set_ylabel("SHAP", fontsize=8)
            ax.tick_params(labelsize=7)
            m = max(0.02, np.abs(sv).max() * 1.25)
            ax.set_ylim(-m, m)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_ism_grid(records, outpath):
    """8 ISM delta heatmaps (L x 20). WT residue marked with a dot."""
    n = len(records)
    ncol = 2
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 3.0 * nrow))
    axes = np.atleast_2d(axes)
    fig.suptitle("In-silico saturation mutagenesis — model 'weights matrix' (Δ binding score vs WT)",
                 fontsize=12, fontweight="bold", y=1.0)
    for k, r in enumerate(records):
        ax = axes[k // ncol, k % ncol]
        delta, pep = r["ism_delta"], r["peptide"]
        vmax = max(0.05, np.abs(delta).max())
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
        im = ax.imshow(delta.T, aspect="auto", cmap="RdBu_r", norm=norm)
        ax.set_yticks(range(20)); ax.set_yticklabels(list(mx.AA), fontsize=6)
        ax.set_xticks(range(len(pep)))
        ax.set_xticklabels([f"P{i+1}" if i < len(pep)-1 else "PΩ" for i in range(len(pep))],
                           fontsize=7)
        for i, a in enumerate(pep):  # mark WT residue
            ax.plot(i, mx.AA_INDEX[a], "o", ms=4, mfc="none", mec="black", mew=1.1)
        ax.set_title(f"{r['group']}/{r['call']}  {pep}  {r['allele']}", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    for k in range(n, nrow * ncol):
        axes[k // ncol, k % ncol].axis("off")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_motif_comparison(df, allele_col, call_col, records, outdir):
    """For each (allele, length) among selected peptides: model motif vs empirical motif."""
    try:
        import logomaker
    except Exception:
        logomaker = None

    seen = {}
    for r in records:
        seen.setdefault((r["allele"], len(r["peptide"])), []).append(r)

    outpaths = []
    for (allele, L), recs in seen.items():
        # model motif = residue frequencies among peptides the MODEL scores as top binders
        backend = recs[0]["_backend"]
        model_pwm, n_model = mx.model_motif_by_sampling(backend, allele, L, seed=7)
        emp_pwm, n_emp = mx.empirical_pwm(df, allele_col, call_col, allele, L)

        model_bits = mx.pwm_to_bits(model_pwm)
        emp_bits = mx.pwm_to_bits(emp_pwm)
        cols = list(mx.AA)

        fig, axes = plt.subplots(2, 1, figsize=(1.05 * L + 2, 5.2))
        safe = allele.replace('*', '').replace(':', '').replace('/', '')
        if logomaker is not None:
            for ax, bits, ttl in [
                (axes[0], model_bits, f"MHCflurry model motif (top model-scored {L}-mers, n={n_model}) — {allele}"),
                (axes[1], emp_bits, f"Observed motif — {allele} strong binders in sample (n={n_emp})")]:
                dfm = pd.DataFrame(bits, columns=cols)  # RangeIndex 0..L-1 for logomaker
                logomaker.Logo(dfm, ax=ax, color_scheme="chemistry")
                ax.set_title(ttl, fontsize=9)
                ax.set_ylabel("bits", fontsize=8)
                ax.set_ylim(0, max(0.5, np.log2(20)))
                ax.set_xticks(range(L))
                ax.set_xticklabels([f"P{i+1}" if i < L-1 else "PΩ" for i in range(L)],
                                   fontsize=8)
        else:  # fallback: heatmaps if logomaker missing
            for ax, pwm, ttl in [(axes[0], model_pwm, f"Model motif {allele} {L}-mer"),
                                 (axes[1], emp_pwm, f"Observed motif {allele} {L}-mer (n={n_emp})")]:
                im = ax.imshow(pwm.T, aspect="auto", cmap="viridis")
                ax.set_yticks(range(20)); ax.set_yticklabels(cols, fontsize=6)
                ax.set_title(ttl, fontsize=9)
                fig.colorbar(im, ax=ax, fraction=0.025)

        # per-position cosine similarity between model and observed preference
        cos = []
        for i in range(L):
            a, b = model_pwm[i], emp_pwm[i]
            c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
            cos.append(c)
        fig.suptitle(f"Motif preference: model vs observed  |  mean per-position cosine = {np.mean(cos):.2f}",
                     fontsize=10, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        op = os.path.join(outdir, f"motif_{safe}_{L}mer.png")
        fig.savefig(op, dpi=160, bbox_inches="tight")
        plt.close(fig)
        outpaths.append((op, allele, L, float(np.mean(cos)), n_emp))
    return outpaths


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="MHCflurry Module-3 explainability layer")
    ap.add_argument("--input", required=True, help="mhcflurry_hlaI.tsv")
    ap.add_argument("--outdir", default="explain_out")
    ap.add_argument("--backend", choices=["mhcflurry", "surrogate"], default="mhcflurry")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--n-per-stratum", type=int, default=2)
    ap.add_argument("--group-col", default="category")
    ap.add_argument("--allele-col", default="binding_allele")
    ap.add_argument("--call-col", default="mhcflurry_binder_call")
    ap.add_argument("--aff-col", default="mhcflurry_affinity")
    ap.add_argument("--background", choices=["dataset", "random"], default="dataset")
    ap.add_argument("--bg-size", type=int, default=120)
    ap.add_argument("--shap-nsamples", default="auto")
    args = ap.parse_args()

    # resolve shap nsamples ("auto" -> integer heuristic; else int)
    if str(args.shap_nsamples).lower() == "auto":
        args.shap_nsamples = None  # set per-peptide below
    else:
        args.shap_nsamples = int(args.shap_nsamples)

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input, sep="\t")
    df["peptide"] = df["peptide"].astype(str).str.strip().str.upper()
    print(f"[load] {len(df)} rows, {df[args.allele_col].nunique()} alleles")

    backend = mx.load_backend(args.backend)
    print(f"[backend] {args.backend}")

    chosen = mx.select_peptides(df, args.group_col, args.allele_col, args.call_col,
                                args.aff_col, n_per_stratum=args.n_per_stratum, seed=args.seed)
    print(f"[select] {len(chosen)} peptides:")
    for s in chosen:
        print(f"    {s.group:9s} {s.call:2s}  {s.peptide:12s} {s.allele:14s} {s.affinity:8.1f} nM")

    records = []
    for s in chosen:
        L = len(s.peptide)
        if not set(s.peptide) <= set(mx.AA):
            print(f"    [skip] {s.peptide}: non-canonical residues", file=sys.stderr); continue
        bg = mx.build_background_pool(df, args.allele_col, s.allele, L,
                                      n=args.bg_size, seed=args.seed,
                                      mode=args.background, call_col=args.call_col)
        nsamp = args.shap_nsamples if args.shap_nsamples is not None else (2 * L + 512)
        sv, base, wt_score = mx.positional_shap(backend, s.peptide, s.allele, bg,
                                                nsamples=nsamp)
        M, delta, wt = mx.ism_matrix(backend, s.peptide, s.allele)
        records.append(dict(peptide=s.peptide, allele=s.allele, group=s.group, call=s.call,
                            affinity=s.affinity, shap=sv, base=base, wt_score=wt_score,
                            ism_full=M, ism_delta=delta, _backend=backend))

    # ---- figures ----
    plot_shap_grid(records, os.path.join(args.outdir, "positional_shap_grid.png"))
    plot_ism_grid(records, os.path.join(args.outdir, "ism_weights_grid.png"))
    motifs = plot_motif_comparison(df, args.allele_col, args.call_col, records, args.outdir)
    print(f"[motif] wrote {len(motifs)} model-vs-observed motif panels")

    # ---- tables ----
    rows = []
    for r in records:
        for i, (a, v) in enumerate(zip(r["peptide"], r["shap"])):
            rows.append(dict(peptide=r["peptide"], allele=r["allele"], group=r["group"],
                             call=r["call"], position=i + 1, residue=a, shap=v,
                             wt_score=r["wt_score"], affinity_nM=r["affinity"]))
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "positional_shap.tsv"),
                              sep="\t", index=False)
    for r in records:
        pd.DataFrame(r["ism_delta"], columns=list(mx.AA),
                     index=[f"P{i+1}" for i in range(len(r['peptide']))]).to_csv(
            os.path.join(args.outdir, f"ism_delta_{r['peptide']}_{r['allele'].replace('*','').replace(':','')}.csv"))
    pd.DataFrame([dict(peptide=r["peptide"], allele=r["allele"], group=r["group"],
                       call=r["call"], affinity_nM=r["affinity"], wt_score=r["wt_score"])
                  for r in records]).to_csv(
        os.path.join(args.outdir, "selected_peptides.tsv"), sep="\t", index=False)

    print(f"[done] outputs in {args.outdir}/")
    return records, motifs


if __name__ == "__main__":
    main()
