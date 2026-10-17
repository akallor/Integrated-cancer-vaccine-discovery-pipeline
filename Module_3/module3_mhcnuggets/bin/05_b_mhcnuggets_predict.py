#!/usr/bin/env python3
"""Module 3b -- MHCnuggets (class I) binding predictor, percentile-rank calibrated.

Parallel HLA-I arm to the MHCflurry script. MHCnuggets outputs IC50 (nM) per allele, but the
raw IC50 is NOT comparable across alleles: different alleles use different model heads
(BA vs BA_to_HLAp) on different output scales, so a naive "best (min-IC50) allele" collapses
onto whichever allele's model runs low, and an absolute 50/500 nM call disagrees wildly with
MHCflurry's percentile call.

Fix (this version): for every allele we score a length-matched background of random peptides
with the SAME per-allele model, convert each observed IC50 to an allele-normalised PERCENTILE
RANK, pick the best allele by RANK (not raw IC50), and call SB/WB/NB on the rank
(SB <= 0.5%, WB <= 2% -- the MHCflurry convention). Rank is scale-independent, so it removes
both the cross-allele scale bias and the incomparability with the MHCflurry arm.

Output columns:
  peptide, sample, tissue, category, noncanonical_origin, binding_allele,
  mhcnuggets_ic50, mhcnuggets_rank, mhcnuggets_binder_call
"""
import sys, os, argparse, tempfile
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mod3_common as C


def to_mhcnuggets_allele(a):
    """'HLA-A*02:01' -> MHCnuggets 'HLA-A02:01' (asterisk removed)."""
    return str(a).replace("*", "").strip()


def run_mhcnuggets(peptides, allele_mn, class_, **predict_kwargs):
    """Predict IC50 (nM) for one allele over a list of peptides. Returns {peptide: ic50}.
    predict_kwargs (e.g. ba_models=True) force a single model type where the installed
    MHCnuggets supports it -- see mod3_common.mhcnuggets_ba_kwargs()."""
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


def predict_binding(uniq, haplo, scorer, *, bg_per_len=500, sb_rank=0.5, wb_rank=2.0,
                    predict_kwargs=None, log=sys.stderr):
    """Core prediction, decoupled from MHCnuggets for testability.

    scorer(peptides, allele_mn, class_, **kw) -> {peptide: ic50}
    Returns a DataFrame: sample, pep_bare, binding_allele, mhcnuggets_ic50, mhcnuggets_rank
    (best allele chosen by percentile rank, length-matched, allele-normalised).
    """
    predict_kwargs = predict_kwargs or {}
    lengths = sorted({len(p) for p in uniq["pep_bare"].unique()})
    ranker_cache, ic_cache = {}, {}   # allele_mn -> ranker / {pep: ic50}   (background is sample-independent)
    best_rows = []
    for sample in sorted(uniq["sample"].unique()):
        peps = sorted(uniq.loc[uniq["sample"] == sample, "pep_bare"].unique())
        alleles = haplo.get(sample, [])
        if not peps or not alleles:
            continue
        for a in alleles:
            mn = to_mhcnuggets_allele(a)
            if mn not in ranker_cache:
                bg = C.random_peptides(lengths, bg_per_len)
                try:
                    bg_ic = scorer(bg, mn, "I", **predict_kwargs)
                    ranker_cache[mn] = C.rank_from_background(bg_ic)
                except Exception as e:
                    print(f"[WARN] background scoring failed for {mn}: {e}", file=log)
                    ranker_cache[mn] = None
            key = (sample, mn)
            if key not in ic_cache:
                try:
                    ic_cache[key] = scorer(peps, mn, "I", **predict_kwargs)
                except Exception as e:
                    print(f"[WARN] {sample}/{a} ({mn}) MHCnuggets failed: {e}", file=log)
                    ic_cache[key] = {}
        # best allele per peptide -> lowest RANK (allele-normalised), not lowest IC50
        for pep in peps:
            best = None   # (allele, ic50, rank)
            for a in alleles:
                mn = to_mhcnuggets_allele(a)
                ic = ic_cache.get((sample, mn), {}).get(pep)
                rk = ranker_cache.get(mn)
                if ic is None or ic != ic or rk is None:
                    continue
                r = rk(ic, len(pep))
                if best is None or r < best[2]:
                    best = (a, float(ic), r)
            if best is not None:
                best_rows.append((sample, pep, best[0], best[1], best[2]))
        print(f"[INFO] {sample}: ranked {len(peps)} peptides x {len(alleles)} alleles "
              f"(bg={bg_per_len}/len, lengths={lengths})", file=log)

    pred = pd.DataFrame(best_rows, columns=["sample", "pep_bare", "binding_allele",
                                            "mhcnuggets_ic50", "mhcnuggets_rank"])
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", default="results/group_fdr_recovery_hlaI.tsv")
    ap.add_argument("--haplotypes", default="PXD013649_hlaI_haplotypes.tsv")
    ap.add_argument("--out", default="results/mhcnuggets_hlaI.tsv")
    ap.add_argument("--pep-col", default="peptide")
    ap.add_argument("--run-col", default="run")
    ap.add_argument("--group-col", default="group")
    ap.add_argument("--annotation", default="build/peptides.tsv")
    ap.add_argument("--tissue-annotation", default=None)
    ap.add_argument("--tumor-token", default="tum")
    ap.add_argument("--normal-token", default="norm")
    ap.add_argument("--bg-per-len", type=int, default=500,
                    help="random background peptides per length for rank calibration (default 500)")
    ap.add_argument("--sb-rank", type=float, default=0.5, help="SB percentile-rank cutoff (default 0.5%%)")
    ap.add_argument("--wb-rank", type=float, default=2.0, help="WB percentile-rank cutoff (default 2%%)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    haplo = C.load_haplotypes(args.haplotypes)
    df = pd.read_csv(args.peptides, sep="\t", dtype=str)
    uniq = C.prepare(df, haplo, pep_col=args.pep_col, run_col=args.run_col,
                     group_col=args.group_col, annotation=args.annotation,
                     tissue_annotation=args.tissue_annotation,
                     tumor_token=args.tumor_token, normal_token=args.normal_token,
                     min_len=8, max_len=15)

    print("\n=== Per-sample plan (MHCnuggets class I, percentile-rank) ===")
    for s in sorted(uniq["sample"].unique()):
        npep = uniq.loc[uniq["sample"] == s, "pep_bare"].nunique()
        print(f"  {s:12s} peptides={npep:6d}  alleles={haplo.get(s, [])}")
    if args.dry_run:
        print("\n[DRY RUN] stopping before MHCnuggets.")
        return

    try:
        import mhcnuggets  # noqa
        from mhcnuggets.src.predict import predict as _pred_fn
    except ImportError:
        sys.exit("mhcnuggets not installed. In the mhcnuggets env: pip install mhcnuggets")

    # force one model head across alleles where supported (else rank still compensates, but warn)
    predict_kwargs = C.mhcnuggets_ba_kwargs(_pred_fn)

    pred = predict_binding(uniq, haplo, run_mhcnuggets, bg_per_len=args.bg_per_len,
                           sb_rank=args.sb_rank, wb_rank=args.wb_rank, predict_kwargs=predict_kwargs)
    if pred.empty:
        sys.exit("No MHCnuggets predictions produced.")

    m = uniq.merge(pred, on=["sample", "pep_bare"], how="left")
    out = pd.DataFrame({
        "peptide": m["pep_bare"], "sample": m["sample"], "tissue": m["tissue"],
        "category": m["category"], "noncanonical_origin": m["noncanonical_origin"],
        "binding_allele": m["binding_allele"], "mhcnuggets_ic50": m["mhcnuggets_ic50"],
        "mhcnuggets_rank": m["mhcnuggets_rank"],
    })
    out["mhcnuggets_binder_call"] = out["mhcnuggets_rank"].map(
        lambda r: rank_binder_call(r, args.sb_rank, args.wb_rank))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {args.out}: {len(out)} (peptide, sample, tissue) rows "
          f"(call on percentile rank: SB<={args.sb_rank}%, WB<={args.wb_rank}%).")
    for cat in ["canonical", "non-canonical"]:
        sub = out[out.category == cat]
        if len(sub):
            print(f"  {cat:14s} n={len(sub):6d}  calls: {sub['mhcnuggets_binder_call'].value_counts().to_dict()}")
    ba = out["binding_allele"].value_counts().to_dict()
    print(f"  best-allele (by rank): {ba}")
    print("  (Now comparable to the MHCflurry arm: both call binders on allele-normalised "
          "percentile rank; canonical MS ligands should be mostly SB/WB.)")


if __name__ == "__main__":
    main()
