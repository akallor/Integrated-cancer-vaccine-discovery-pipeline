#!/usr/bin/env python3
"""mod3_aggregate.py -- derive the diagnostic-plot summary tables from a combined
Module-3 prediction table.

>>> INFERRED SCHEMA <<<  No aggregation script was provided; the four static plot
scripts read fixed-name pre-aggregated TSVs. This produces exactly those files.
Predictor-agnostic: auto-detects the *_binder_call column, so it works for both
MHCflurry and MHCnuggets prediction outputs. Verify/replace with your own if you
already have one.
"""
import os, argparse
import pandas as pd

def pick(cols, *cands, contains=None):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    if contains:
        for c in cols:
            if contains in c.lower():
                return c
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    df = pd.read_csv(a.input, sep=None, engine="python")
    pep      = pick(df.columns, "peptide")
    sample   = pick(df.columns, "sample")
    tissue   = pick(df.columns, "tissue")
    category = pick(df.columns, "category")
    allele   = pick(df.columns, "binding_allele", "allele")
    call     = pick(df.columns, contains="binder_call")
    assert pep and sample and category and allele and call, \
        f"aggregate: missing required columns in {list(df.columns)}"

    df["_cat"]  = df[category].map(lambda v: "non_canonical" if str(v).strip().lower().startswith("non") else "canonical")
    df["_call"] = df[call].astype(str).str.upper().str.strip()

    # 1) sample x {SB,WB,NB} unique-peptide counts
    bc = (df.groupby([sample, "_call"])[pep].nunique().unstack(fill_value=0)
            .reindex(columns=["SB", "WB", "NB"], fill_value=0))
    bc.index.name = "sample"
    bc.reset_index().to_csv(f"{a.outdir}/df_binder_raw_counts.tsv", sep="\t", index=False)

    # 2) sample, peptide(total unique), per-allele unique counts
    per_allele = df.groupby([sample, allele])[pep].nunique().unstack(fill_value=0)
    total = df.groupby(sample)[pep].nunique().rename("peptide")
    ac = pd.concat([total, per_allele], axis=1).fillna(0).astype(int)
    ac.index.name = "sample"
    ac.reset_index().to_csv(f"{a.outdir}/df_allele_pep_counts.tsv", sep="\t", index=False)

    # 3) sample x {canonical, non_canonical}
    og = (df.groupby([sample, "_cat"])[pep].nunique().unstack(fill_value=0)
            .reindex(columns=["canonical", "non_canonical"], fill_value=0))
    og.index.name = "sample"
    og.reset_index().to_csv(f"{a.outdir}/df_origin_grp.tsv", sep="\t", index=False)

    # 4) tissue x {canonical, non_canonical}
    if tissue:
        tg = (df.groupby([tissue, "_cat"])[pep].nunique().unstack(fill_value=0)
                .reindex(columns=["canonical", "non_canonical"], fill_value=0))
        tg.index.name = "tissue"
        tg.reset_index().to_csv(f"{a.outdir}/df_tissue_grp.tsv", sep="\t", index=False)

    print("aggregated -> df_binder_raw_counts.tsv, df_allele_pep_counts.tsv, "
          "df_origin_grp.tsv, df_tissue_grp.tsv")

if __name__ == "__main__":
    main()
