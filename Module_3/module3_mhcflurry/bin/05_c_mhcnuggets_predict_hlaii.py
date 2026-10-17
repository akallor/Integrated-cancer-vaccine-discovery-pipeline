#!/usr/bin/env python3
"""Module 3c -- MHCnuggets (class II) binding predictor against a COMMON-ALLELE PANEL.

>>> IMPORTANT DISCLAIMER <<<
This study (Chong 2020 / PXD013649) provides NO sample-specific HLA-II (DR/DQ/DP) typing -- the
HLA-II immunopeptidome was captured with pan-HLA-II antibodies. Therefore this module predicts
against a GENERIC panel of common HLA-II alleles (passed via --alleles), NOT the samples' own
alleles. Results are a population-level "could this bind a common class-II allele" signal, useful
as a coarse feature ONLY. For accurate, sample-matched HLA-II binding you MUST supply real DR/DQ/DP
typing (e.g. re-type from RNA-seq/exome with a class-II-capable HLA typer). Class-II binding
prediction is also weaker than class I across all tools -- lean on the MS presentation evidence.

Reuses mod3_common for sample/tissue/origin/category. Because the allele panel is sample-
independent, prediction is done once over all unique peptides per panel allele, then mapped back
to (sample, tissue, peptide) rows.

Output columns:
  peptide, sample, tissue, category, noncanonical_origin, binding_allele,
  binding_allele_source, mhcnuggets_ic50, mhcnuggets_binder_call

Run in the mhcnuggets env:
  python 03c_mhcnuggets_hlaII.py --peptides results/group_fdr_recovery_hlaII.tsv \
      --haplotypes PXD013649_hlaI_haplotypes.tsv \
      --alleles common_hlaII_panel.txt --out results/mhcnuggets_hlaII.tsv
"""
import sys, os, argparse, tempfile
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mod3_common as C

SB_NM, WB_NM = 50.0, 500.0

# Default top-10 common HLA-II alleles (global DRB1 -- most frequent & best-supported in MHCnuggets).
# Override with --alleles. DQ/DP can be added but require alpha-beta pairing (see README).
DEFAULT_PANEL = [
    "HLA-DRB1*07:01", "HLA-DRB1*15:01", "HLA-DRB1*03:01", "HLA-DRB1*01:01", "HLA-DRB1*11:01",
    "HLA-DRB1*13:01", "HLA-DRB1*04:01", "HLA-DRB1*08:01", "HLA-DRB1*13:02", "HLA-DRB1*11:04",
]


def to_mhcnuggets_allele(a):
    return str(a).replace("*", "").strip()


def load_panel(spec):
    """--alleles is either a path (one allele per line) or a comma-separated string. Falls back to default."""
    if not spec:
        return list(DEFAULT_PANEL)
    if os.path.exists(spec):
        with open(spec) as fh:
            return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    return [a.strip() for a in spec.split(",") if a.strip()]


def run_mhcnuggets(peptides, allele_mn, class_, **predict_kwargs):
    from mhcnuggets.src.predict import predict
    pf = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    pf.write("\n".join(peptides) + "\n"); pf.close()
    out_path = pf.name + ".out.csv"
    try:
        predict(class_=class_, peptides_path=pf.name, mhc=allele_mn, output=out_path, **predict_kwargs)
        res = pd.read_csv(out_path)
        ic_col = next((c for c in res.columns if "ic50" in c.lower()), None)
        pep_col = next((c for c in res.columns if "pep" in c.lower()), res.columns[0])
        if ic_col is None:
            print(f"[WARN] no ic50 column for {allele_mn}; cols={list(res.columns)}", file=sys.stderr)
            return {}
        return dict(zip(res[pep_col].astype(str), pd.to_numeric(res[ic_col], errors="coerce")))
    finally:
        for p in (pf.name, out_path):
            if os.path.exists(p):
                os.remove(p)


def rank_binder_call(rank, sb_rank, wb_rank):
    try:
        r = float(rank)
    except (TypeError, ValueError):
        return "NA"
    return "SB" if r <= sb_rank else ("WB" if r <= wb_rank else "NB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", default="results/group_fdr_recovery_hlaII.tsv")
    ap.add_argument("--haplotypes", default="PXD013649_hlaI_haplotypes.tsv",
                    help="used ONLY for the sample list / sample-name matching, not for alleles")
    ap.add_argument("--alleles", default=None,
                    help="common HLA-II panel: a file (one allele/line) or comma-separated list. "
                         "Format like HLA-DRB1*15:01. Defaults to a top-10 common DRB1 panel.")
    ap.add_argument("--out", default="results/mhcnuggets_hlaII.tsv")
    ap.add_argument("--pep-col", default="peptide")
    ap.add_argument("--run-col", default="run")
    ap.add_argument("--group-col", default="group")
    ap.add_argument("--annotation", default="build/peptides.tsv")
    ap.add_argument("--tissue-annotation", default=None)
    ap.add_argument("--tumor-token", default="tum")
    ap.add_argument("--normal-token", default="norm")
    ap.add_argument("--min-len", type=int, default=9, help="class-II peptides are longer; default 9")
    ap.add_argument("--max-len", type=int, default=30)
    ap.add_argument("--bg-per-len", type=int, default=300,
                    help="random background peptides per length for rank calibration")
    ap.add_argument("--sb-rank", type=float, default=0.5)
    ap.add_argument("--wb-rank", type=float, default=2.0)
    ap.add_argument("--max-bg-binder-rate", type=float, default=None,
                    help="drop panel alleles whose random-background binder rate (%%IC50<500) exceeds this")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    panel = load_panel(args.alleles)
    print("=" * 78, file=sys.stderr)
    print("DISCLAIMER: predicting HLA-II binding against a GENERIC common-allele panel, NOT the",
          file=sys.stderr)
    print("samples' own DR/DQ/DP types (this dataset has no HLA-II typing). Coarse feature only;",
          file=sys.stderr)
    print("supply real sample-specific HLA-II typing for accurate binding calls.", file=sys.stderr)
    print(f"Panel ({len(panel)} alleles): {panel}", file=sys.stderr)
    print("=" * 78, file=sys.stderr)

    haplo = C.load_haplotypes(args.haplotypes)
    df = pd.read_csv(args.peptides, sep="\t", dtype=str)
    uniq = C.prepare(df, haplo, pep_col=args.pep_col, run_col=args.run_col,
                     group_col=args.group_col, annotation=args.annotation,
                     tissue_annotation=args.tissue_annotation,
                     tumor_token=args.tumor_token, normal_token=args.normal_token,
                     min_len=args.min_len, max_len=args.max_len)

    all_peps = sorted(uniq["pep_bare"].unique())
    print(f"\n=== Plan: {len(all_peps)} unique peptides x {len(panel)} panel alleles ===")
    if args.dry_run:
        print("[DRY RUN] stopping before MHCnuggets.")
        return
    try:
        import mhcnuggets  # noqa
    except ImportError:
        sys.exit("mhcnuggets not installed. In the mhcnuggets env: pip install mhcnuggets")

    from mhcnuggets.src.predict import predict as _pred_fn
    predict_kwargs = C.mhcnuggets_ba_kwargs(_pred_fn)
    lengths = sorted({len(p) for p in all_peps})

    # allele panel is sample-independent -> per allele: background rank + observed IC50 (once)
    per_allele_ic, per_allele_rank, bg_binder_rate = {}, {}, {}
    for a in panel:
        mn = to_mhcnuggets_allele(a)
        try:
            bg = C.random_peptides(lengths, args.bg_per_len)
            bg_ic = run_mhcnuggets(bg, mn, "II", **predict_kwargs)
            per_allele_rank[a] = C.rank_from_background(bg_ic)
            per_allele_ic[a] = run_mhcnuggets(all_peps, mn, "II", **predict_kwargs)
            vals = sorted(float(v) for v in bg_ic.values() if v == v)
            rate = 100.0 * sum(v < 500 for v in vals) / len(vals) if vals else float("nan")
            bg_binder_rate[a] = rate
            med = vals[len(vals) // 2] if vals else float("nan")
            flag = "  <-- NON-SPECIFIC: model calls random peptides as binders; ranks/calls unreliable" if rate > 10 else ""
            print(f"[QC] {a} ({mn}): background %IC50<500nM={rate:.1f}  median={med:.0f}{flag}", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] panel allele {a} ({mn}) failed: {e}", file=sys.stderr)

    best_rows = []
    for pep in all_peps:
        best = None   # (allele, ic50, rank)
        for a in per_allele_ic:
            if args.max_bg_binder_rate is not None and bg_binder_rate.get(a, 0.0) > args.max_bg_binder_rate:
                continue
            ic = per_allele_ic[a].get(pep); rk = per_allele_rank.get(a)
            if ic is None or ic != ic or rk is None:
                continue
            r = rk(ic, len(pep))
            if best is None or r < best[2]:
                best = (a, float(ic), r)
        if best is not None:
            best_rows.append((pep, best[0], best[1], best[2]))
    if not best_rows:
        sys.exit("No MHCnuggets class-II predictions produced.")
    pred = pd.DataFrame(best_rows, columns=["pep_bare", "binding_allele", "mhcnuggets_ic50", "mhcnuggets_rank"])

    m = uniq.merge(pred, on="pep_bare", how="left")
    out = pd.DataFrame({
        "peptide": m["pep_bare"],
        "sample": m["sample"],
        "tissue": m["tissue"],
        "category": m["category"],
        "noncanonical_origin": m["noncanonical_origin"],
        "binding_allele": m["binding_allele"],
        "binding_allele_source": "common_panel_NOT_sample_specific",   # unmissable in the data
        "mhcnuggets_ic50": m["mhcnuggets_ic50"],
        "mhcnuggets_rank": m["mhcnuggets_rank"],
    })
    out["mhcnuggets_binder_call"] = out["mhcnuggets_rank"].map(
        lambda r: rank_binder_call(r, args.sb_rank, args.wb_rank))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {args.out}: {len(out)} (peptide, sample, tissue) rows.")
    print("  REMINDER: binding_allele is from a GENERIC panel, not sample-specific typing.")


if __name__ == "__main__":
    main()
