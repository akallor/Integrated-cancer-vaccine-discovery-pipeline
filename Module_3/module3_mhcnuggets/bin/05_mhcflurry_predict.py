#!/usr/bin/env python3
"""Module 3 -- MHCflurry (class I) binding affinity + antigen processing + presentation.

Runs the MHCflurry Class1PresentationPredictor on ALL recovered HLA-I peptides (canonical AND
non-canonical -- do not pre-filter; canonical peptides are your calibration set, and
presented-but-predicted-non-binder non-canonical peptides are findings, not noise).

Per-peptide it emits the best allele across the ORIGIN SAMPLE's own HLA-I genotype, plus
affinity, affinity percentile, processing score, presentation score, presentation percentile.

Sample assignment: each PSM's `run` (filename) is matched to a Sample in the haplotype table.
This is deliberately defensive because the haplotype table's sample names DO NOT cleanly match
the filenames:
  - letter-O vs zero:      0D5P/0NVC/0MM745 in filenames  vs  ...O...  in the table
  - transposition typo:    table "05DP"  ==  sample 0D5P (Chong flagship)
  - hyphens:               table "C3N-02289"  vs filename "C3N02289"
Matching normalizes O->0 and strips non-alphanumerics, plus an explicit alias map for the
transposition. Unmatched / ambiguous runs are REPORTED, not silently dropped.

Usage:
  python 05_mhcflurry_predict.py \
      --peptides results/group_fdr_recovery_hlaI.tsv \
      --haplotypes PXD013649_hlaI_haplotypes.tsv \
      --out results/mhcflurry_hlaI.tsv \
      [--dry-run]        # verify sample/allele/peptide assignment WITHOUT running MHCflurry

Prereqs (in the mhcflurry env):
  pip install mhcflurry && mhcflurry-downloads fetch models_class1_presentation
"""
import sys, os, re, argparse
import pandas as pd

AA = set("ACDEFGHIKLMNPQRSTVWY")
VALID_PEP = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]{8,15}$')   # MHCflurry class I supports 8-15mers

# Explicit aliases for haplotype-table sample names whose spelling can't be normalized to the
# filename token automatically (transposition). Value = list of tokens that appear in filenames.
SAMPLE_ALIASES = {
    "05DP": ["0D5P", "OD5P"],          # table typo/transposition -> Chong flagship 0D5P
    "OMM745": ["0MM745", "OMM745"],
    "ONVC": ["0NVC", "ONVC"],
    "C3N-02289": ["C3N02289", "C3N-02289"],
    "C3N-02671": ["C3N02671", "C3N-02671"],
}


def norm(s):
    """Uppercase, collapse letter-O to zero, strip non-alphanumerics (removes hyphens/underscores)."""
    return re.sub(r'[^A-Z0-9]', '', str(s).upper().replace('O', '0'))


def strip_peptide(p):
    """Bare residues only: drop flanks (K.PEP.L), mod annotations [..]/(..), non-residue chars.
    NOTE: keeps I and L distinct (MHCflurry needs the true sequence -- do NOT I->L here)."""
    s = str(p).strip().upper()
    if s.count(".") >= 2:                       # K.PEPTIDE.L -> PEPTIDE
        parts = s.split(".")
        if len(parts) > 2:
            s = ".".join(parts[1:-1])
    s = re.sub(r"\[[^\]]*\]", "", s)            # [+15.9949]
    s = re.sub(r"\([^)]*\)", "", s)             # (ox)
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    return s


def load_haplotypes(path):
    """Return {sample: [unique allele strings]} from the Sample<TAB>Haplotype table."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    out = {}
    for _, r in df.iterrows():
        sample = str(r["Sample"]).strip()
        alleles = [a.strip() for a in str(r["Haplotype"]).split(";") if a.strip()]
        # dedupe while preserving order (homozygous samples list an allele twice)
        seen, uniq = set(), []
        for a in alleles:
            if a not in seen:
                seen.add(a); uniq.append(a)
        out[sample] = uniq
    return out


def build_sample_tokens(samples):
    """Precompute normalized match-tokens per sample (using aliases where defined)."""
    tokens = {}
    for s in samples:
        cands = SAMPLE_ALIASES.get(s, [s])
        tokens[s] = [norm(c) for c in cands]
    return tokens


def assign_sample(run, sample_tokens):
    """Return (sample, ambiguous). Match any normalized token as a substring of the normalized run."""
    nrun = norm(run)
    hits = [s for s, toks in sample_tokens.items() if any(t and t in nrun for t in toks)]
    if len(hits) == 1:
        return hits[0], False
    if len(hits) > 1:
        return hits[0], True     # ambiguous: report but take the first
    return None, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", default="results/group_fdr_recovery_hlaI.tsv")
    ap.add_argument("--haplotypes", default="PXD013649_hlaI_haplotypes.tsv")
    ap.add_argument("--out", default="results/mhcflurry_hlaI.tsv")
    ap.add_argument("--pep-col", default="peptide", help="column holding the peptide (default: peptide)")
    ap.add_argument("--run-col", default="run", help="column holding the run/filename (default: run)")
    ap.add_argument("--group-col", default="group", help="CANON/NONCANON column (default: group)")

    # --- (a) non-canonical origin annotation (peptide-keyed, intrinsic to the sequence) ---
    ap.add_argument("--annotation", default="build/peptides.tsv",
                    help="Peptide-origin annotation file (TSV). Preferred schema: two columns, "
                         "'peptide' and 'noncanonical_origin' (e.g. lncRNA/TE/uORF). Falls back to "
                         "a 'source' column (as in build/peptides.tsv). Matching is I/L-insensitive. "
                         "Origin is intrinsic to the peptide, so it is keyed on sequence, not sample. "
                         "Optional -- if absent, non-canonical peptides get 'unknown'.")

    # --- (b) tumor/normal tissue axis (SAMPLE/RUN-keyed, NOT peptide-keyed) ---
    ap.add_argument("--tissue-annotation", default=None,
                    help="Optional TSV mapping sample -> tissue. Schema: two columns, 'sample' and "
                         "'tissue' (values: tumor/normal). Overrides/supplies the run-substring "
                         "detection. Use when filenames don't encode tissue.")
    ap.add_argument("--tumor-token", default="tum",
                    help="case-insensitive substring in the run marking TUMOR (default: 'tum')")
    ap.add_argument("--normal-token", default="norm",
                    help="case-insensitive substring in the run marking NORMAL (default: 'norm')")

    ap.add_argument("--dry-run", action="store_true",
                    help="report sample/allele/peptide assignment and exit (no MHCflurry)")
    args = ap.parse_args()

    # ---- load inputs ----
    haplo = load_haplotypes(args.haplotypes)
    sample_tokens = build_sample_tokens(list(haplo))
    df = pd.read_csv(args.peptides, sep="\t", dtype=str)
    for c in (args.pep_col, args.run_col):
        if c not in df.columns:
            sys.exit(f"Column '{c}' not in {args.peptides}. Columns: {list(df.columns)}")

    n_in = len(df)

    # ---- assign sample from run ----
    assign = df[args.run_col].map(lambda r: assign_sample(r, sample_tokens))
    df["sample"] = [a[0] for a in assign]
    df["_ambiguous"] = [a[1] for a in assign]

    unmatched = df["sample"].isna()
    ambiguous = df["_ambiguous"]
    if unmatched.any():
        runs = sorted(df.loc[unmatched, args.run_col].dropna().unique())[:10]
        print(f"[WARN] {unmatched.sum()} PSMs matched NO sample. Example runs:", file=sys.stderr)
        for r in runs:
            print(f"        {r}", file=sys.stderr)
    if ambiguous.any():
        print(f"[WARN] {int(ambiguous.sum())} PSMs matched MULTIPLE samples (took first). "
              f"Check SAMPLE_ALIASES for token collisions.", file=sys.stderr)
    df = df[df["sample"].notna()].copy()

    # ---- (b) tissue axis: tumor vs normal, per RUN (a sample can be BOTH) ----
    # Primary signal = run substring (tum/norm); optional sample->tissue file overrides/fills.
    tumtok, normtok = args.tumor_token.lower(), args.normal_token.lower()
    def tissue_from_run(run):
        r = str(run).lower()
        has_t, has_n = (tumtok in r), (normtok in r)
        if has_t and not has_n: return "tumor"
        if has_n and not has_t: return "normal"
        if has_t and has_n:     return "ambiguous"   # both tokens present -> flag
        return "unknown"
    df["tissue"] = df[args.run_col].map(tissue_from_run)

    tissue_override = {}
    if args.tissue_annotation:
        if os.path.exists(args.tissue_annotation):
            ta = pd.read_csv(args.tissue_annotation, sep="\t", dtype=str)
            scol = next((c for c in ta.columns if c.lower() == "sample"), None)
            tcol = next((c for c in ta.columns if c.lower() == "tissue"), None)
            if scol and tcol:
                tissue_override = {str(r[scol]).strip(): str(r[tcol]).strip().lower()
                                   for _, r in ta.iterrows()}
                # apply override ONLY where the run didn't already give a clear tumor/normal
                mask = df["tissue"].isin(["unknown", "ambiguous"])
                df.loc[mask, "tissue"] = df.loc[mask, "sample"].map(
                    lambda s: tissue_override.get(s, "unknown"))
                print(f"[INFO] applied tissue overrides for {len(tissue_override)} samples "
                      f"from {args.tissue_annotation}", file=sys.stderr)
            else:
                print(f"[WARN] {args.tissue_annotation} lacks 'sample'/'tissue' columns "
                      f"({list(ta.columns)}); using run-substring tissue only.", file=sys.stderr)
        else:
            print(f"[WARN] tissue annotation {args.tissue_annotation} not found; "
                  f"using run-substring tissue only.", file=sys.stderr)
    tvc = df["tissue"].value_counts().to_dict()
    print(f"[INFO] tissue (per PSM): {tvc}", file=sys.stderr)
    if df["tissue"].isin(["unknown", "ambiguous"]).any():
        n_amb = int(df["tissue"].isin(["unknown","ambiguous"]).sum())
        print(f"[WARN] {n_amb} PSMs have unknown/ambiguous tissue "
              f"(neither/both of '{tumtok}'/'{normtok}' in run, no override). "
              f"They are kept with tissue label as-is.", file=sys.stderr)

    # ---- clean peptides + validity filter ----
    df["pep_bare"] = df[args.pep_col].map(strip_peptide)
    valid = df["pep_bare"].map(lambda p: bool(VALID_PEP.match(p)))
    n_bad = int((~valid).sum())
    if n_bad:
        egs = df.loc[~valid, args.pep_col].head(5).tolist()
        print(f"[INFO] Dropping {n_bad} peptides outside 8-15aa / non-standard residues. "
              f"e.g. {egs}", file=sys.stderr)
    df = df[valid].copy()

    # category tag
    if args.group_col in df.columns:
        df["category"] = df[args.group_col].map(
            lambda g: "non-canonical" if str(g).upper() == "NONCANON" else "canonical")
    else:
        df["category"] = "unknown"

    # ---- collapse to unique (sample, tissue, peptide) with a category ----
    #   Keying on tissue too is CRITICAL: the same peptide in a tumor run vs a normal run must
    #   stay as TWO rows -- that split IS the tumor-specificity signal.
    #   prefer non-canonical if a peptide is ever tagged non-canonical (conservative)
    def pick_cat(s):
        return "non-canonical" if (s == "non-canonical").any() else s.iloc[0]
    uniq = (df.groupby(["sample", "tissue", "pep_bare"], as_index=False)
              .agg(category=("category", pick_cat)))
    print(f"[INFO] {n_in} PSMs -> {len(uniq)} unique (sample, tissue, peptide) rows "
          f"across {uniq['sample'].nunique()} samples "
          f"({uniq.groupby('tissue').size().to_dict()} by tissue).", file=sys.stderr)

    # ---- (a) non-canonical origin subtype (lncRNA / TE / ...) from the annotation file ----
    #   canonical peptides get 'canonical_origin'; non-canonical get their annotated origin, else 'unknown'.
    #   Preferred schema: columns 'peptide' + 'noncanonical_origin'. Fallback: 'peptide' + 'source'.
    #   Match on I->L-normalized sequence (MS can't distinguish I/L; build/peptides.tsv is I->L-normalized).
    def il(p):
        return str(p).upper().replace("I", "L")
    origin_map = {}
    if args.annotation and os.path.exists(args.annotation):
        ann = pd.read_csv(args.annotation, sep="\t", dtype=str)
        pcol = next((c for c in ann.columns if c.lower() == "peptide"), None)
        ocol = (next((c for c in ann.columns if c.lower() == "noncanonical_origin"), None)
                or next((c for c in ann.columns if c.lower() == "source"), None))
        if pcol and ocol:
            for _, r in ann.iterrows():
                key = il(strip_peptide(r[pcol]))
                val = str(r[ocol]).strip()
                if key:
                    origin_map[key] = val
            print(f"[INFO] loaded {len(origin_map)} peptide->origin entries from {args.annotation} "
                  f"(peptide col='{pcol}', origin col='{ocol}')", file=sys.stderr)
        else:
            print(f"[WARN] {args.annotation} needs a 'peptide' column and a "
                  f"'noncanonical_origin' (or 'source') column; found {list(ann.columns)}. "
                  f"noncanonical_origin will be 'unknown'.", file=sys.stderr)
    else:
        print(f"[WARN] annotation file {args.annotation} not found; "
              f"noncanonical_origin will be 'unknown' for non-canonical peptides.", file=sys.stderr)

    def origin_of(row):
        if row["category"] == "canonical":
            return "canonical_origin"
        return origin_map.get(il(row["pep_bare"]), "unknown")
    uniq["noncanonical_origin"] = uniq.apply(origin_of, axis=1)

    # ---- per-sample plan (+ dry-run report) ----
    print("\n=== Per-sample plan ===")
    plan = {}
    for sample in sorted(uniq["sample"].unique()):
        peps = sorted(uniq.loc[uniq["sample"] == sample, "pep_bare"].unique())
        alleles = haplo.get(sample, [])
        plan[sample] = (peps, alleles)
        print(f"  {sample:12s} peptides={len(peps):6d}  alleles={alleles}")

    if args.dry_run:
        print("\n[DRY RUN] stopping before MHCflurry. Review sample/allele/peptide assignment above.")
        return

    # ---- MHCflurry ----
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        sys.exit("mhcflurry not installed. In the mhcflurry env:\n"
                 "  pip install mhcflurry && mhcflurry-downloads fetch models_class1_presentation")
    predictor = Class1PresentationPredictor.load()

    # supported-allele set (drop rare alleles MHCflurry can't handle, e.g. B*15:220, B*02:01)
    try:
        supported = set(predictor.supported_alleles)
    except Exception:
        try:
            supported = set(predictor.affinity_predictor.supported_alleles)
        except Exception:
            supported = None

    results = []
    for sample, (peps, alleles) in plan.items():
        use = alleles
        if supported is not None:
            use = [a for a in alleles if a in supported]
            dropped = [a for a in alleles if a not in supported]
            if dropped:
                print(f"[WARN] {sample}: MHCflurry does not support {dropped} -- excluded.",
                      file=sys.stderr)
        if not use:
            print(f"[WARN] {sample}: no supported alleles -- skipping {len(peps)} peptides.",
                  file=sys.stderr)
            continue
        if not peps:
            continue
        res = predictor.predict(
            peptides=peps,
            alleles=use,                      # flat list -> single genotype, best allele per peptide
            include_affinity_percentile=True,
            verbose=0,
        )
        res = res.rename(columns={"peptide": "pep_bare"})
        res["sample"] = sample
        results.append(res)

    if not results:
        sys.exit("No predictions produced (no supported alleles for any sample?).")
    pred = pd.concat(results, ignore_index=True)

    # ---- assemble the final tagged table ----
    pred = pred.merge(uniq, on=["sample", "pep_bare"], how="left")

    out = pd.DataFrame({
        "peptide": pred["pep_bare"],
        "sample": pred["sample"],
        "tissue": pred.get("tissue"),
        "category": pred["category"],
        "noncanonical_origin": pred.get("noncanonical_origin"),
        "binding_allele": pred.get("best_allele"),
        "mhcflurry_affinity": pred.get("affinity"),
        "mhcflurry_affinity_percentile": pred.get("affinity_percentile"),
        "mhcflurry_processing_score": pred.get("processing_score"),
        "mhcflurry_presentation_score": pred.get("presentation_score"),
        "mhcflurry_presentation_percentile": pred.get("presentation_percentile"),
    })

    # convenience binder call on the affinity percentile-rank (SB<=0.5%, WB<=2%) -- for downstream ranking
    def binder(pct):
        try:
            pct = float(pct)
        except (TypeError, ValueError):
            return "NA"
        return "SB" if pct <= 0.5 else ("WB" if pct <= 2.0 else "NB")
    out["mhcflurry_binder_call"] = out["mhcflurry_affinity_percentile"].map(binder)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    # ---- summary (canonical peptides are the calibration control: most should be binders) ----
    print(f"\nWrote {args.out}: {len(out)} (peptide, sample) predictions.")
    for cat in ["canonical", "non-canonical"]:
        sub = out[out.category == cat]
        if len(sub):
            calls = sub["mhcflurry_binder_call"].value_counts().to_dict()
            print(f"  {cat:14s} n={len(sub):6d}  binder calls: {calls}")
    nc = out[out.category == "non-canonical"]
    if len(nc):
        print(f"  non-canonical by origin: {nc['noncanonical_origin'].value_counts().to_dict()}")
    if "tissue" in out.columns:
        print(f"  rows by tissue: {out['tissue'].value_counts().to_dict()}")
        # tumor-specific view: non-canonical strong binders seen in tumor
        tsb = out[(out.category=="non-canonical") & (out.tissue=="tumor")
                  & (out.mhcflurry_binder_call.isin(["SB","WB"]))]
        print(f"  non-canonical tumor binders (SB/WB): {len(tsb)}  "
              f"-- candidate tumor-specific antigens (compare against tissue=='normal' rows).")
    print("  (Sanity check: canonical peptides are MS-confirmed ligands -> should be mostly SB/WB. "
          "If not, re-check the per-sample allele assignment.)")


if __name__ == "__main__":
    main()
