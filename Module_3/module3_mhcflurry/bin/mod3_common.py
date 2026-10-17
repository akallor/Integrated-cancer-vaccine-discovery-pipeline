#!/usr/bin/env python3
"""mod3_common.py -- shared helpers for Module 3 binding predictors (MHCflurry, MHCnuggets).

Encapsulates the sample-matching, tumor/normal tissue axis, non-canonical origin annotation,
peptide cleaning, and the (sample, tissue, peptide) collapse -- so every predictor arm treats
these identically. Import from the per-predictor scripts.
"""
import sys, os, re
import pandas as pd

AA = set("ACDEFGHIKLMNPQRSTVWY")

import bisect as _bisect
import random as _random
_AA_LIST = "ACDEFGHIKLMNPQRSTVWY"


def random_peptides(lengths, per_len, seed=0):
    """A length-matched background of random peptides for percentile-rank calibration.
    `lengths` = iterable of peptide lengths present in the data; `per_len` per length."""
    r = _random.Random(seed)
    return [''.join(r.choice(_AA_LIST) for _ in range(L)) for L in sorted(set(lengths)) for _ in range(per_len)]


def rank_from_background(bg_ic50):
    """Build a percentile-rank function from a background {peptide: ic50} for ONE allele,
    scored with the SAME model MHCnuggets uses for the observed peptides. Returns
    rank_pct(ic50, length) -> % of random peptides that bind at least as strongly (lower
    ic50). Lower rank = stronger. Length-matched (NetMHCpan convention); pools if a length
    is absent. Because the background and the observed peptides go through the same per-allele
    model, ranks are comparable ACROSS alleles even when the raw IC50 scales differ
    (BA vs BA_to_HLAp) -- which is exactly the scale-mixing that breaks a min-IC50 comparison.
    """
    by_len = {}
    for pep, v in bg_ic50.items():
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if fv == fv:
            by_len.setdefault(len(str(pep)), []).append(fv)
    for L in by_len:
        by_len[L].sort()
    pooled = sorted(x for xs in by_len.values() for x in xs)

    def rank_pct(ic50, length):
        arr = by_len.get(length) or pooled
        if not arr:
            return 100.0
        k = _bisect.bisect_left(arr, float(ic50))   # background strictly stronger
        return 100.0 * k / len(arr)
    return rank_pct


# Aliases for haplotype-table sample names whose spelling can't be normalized to the filename
# token automatically (digit transpositions). Value = tokens that appear in filenames.
SAMPLE_ALIASES = {
    "05DP":   ["0D5P", "OD5P"],           # transposition
    "0MM745": ["0MM745", "OMM745", "0MM475", "OMM475"],  # 745<->475 transposition in files
    "OMM745": ["0MM745", "OMM745", "0MM475", "OMM475"],
    "ONVC":   ["0NVC", "ONVC"],
    "0NVC":   ["0NVC", "ONVC"],
    "C3N-02289": ["C3N02289", "C3N-02289"],
    "C3N-02671": ["C3N02671", "C3N-02671"],
}


def norm(s):
    """Uppercase, collapse letter-O to zero, strip non-alphanumerics (removes hyphens/underscores)."""
    return re.sub(r'[^A-Z0-9]', '', str(s).upper().replace('O', '0'))


def il(p):
    """I->L normalization (MS cannot distinguish isoleucine/leucine)."""
    return str(p).upper().replace("I", "L")


def strip_peptide(p):
    """Bare residues only: drop flanks (K.PEP.L), mod annotations [..]/(..), non-residue chars.
    Keeps I and L distinct (predictors need the true sequence)."""
    s = str(p).strip().upper()
    if s.count(".") >= 2:
        parts = s.split(".")
        if len(parts) > 2:
            s = ".".join(parts[1:-1])
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r"\([^)]*\)", "", s)
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    return s


def mhcnuggets_ba_kwargs(predict_fn, verbose=True):
    """Return kwargs that force MHCnuggets to use the binding-affinity (BA) model for ALL alleles,
    if the installed predict() supports it. Prevents silent mixing of BA vs BA_to_HLAp across
    alleles (default falls back to BA only when BA_to_HLAp is missing -> different scales per
    allele). Version-safe: only passes the flag if the signature has it.
    """
    import inspect
    try:
        params = inspect.signature(predict_fn).parameters
    except (TypeError, ValueError):
        return {}
    kw = {"ba_models": True} if "ba_models" in params else {}
    if verbose:
        if kw:
            print("[INFO] forcing MHCnuggets BA model for all alleles (ba_models=True) -> "
                  "consistent scale, no BA_to_HLAp fallback mixing.", file=sys.stderr)
        else:
            print("[WARN] this MHCnuggets version lacks a 'ba_models' param; cannot force one model "
                  "type -- BA vs BA_to_HLAp may be mixed across alleles. Check the run log.",
                  file=sys.stderr)
    return kw


def load_haplotypes(path):
    """Return {sample: [unique allele strings]} from a Sample<TAB>Haplotype table."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    out = {}
    for _, r in df.iterrows():
        sample = str(r["Sample"]).strip()
        alleles = [a.strip() for a in str(r["Haplotype"]).split(";") if a.strip()]
        seen, uniq = set(), []
        for a in alleles:
            if a not in seen:
                seen.add(a); uniq.append(a)
        out[sample] = uniq
    return out


def build_sample_tokens(samples):
    tokens = {}
    for s in samples:
        cands = SAMPLE_ALIASES.get(s, [s])
        tokens[s] = [norm(c) for c in cands]
    return tokens


def assign_sample(run, sample_tokens):
    nrun = norm(run)
    hits = [s for s, toks in sample_tokens.items() if any(t and t in nrun for t in toks)]
    if len(hits) == 1:
        return hits[0], False
    if len(hits) > 1:
        return hits[0], True
    return None, False


def load_origin_map(annotation_path):
    """peptide (I/L-normalized) -> origin subtype. Accepts 'peptide'+'noncanonical_origin'
    (preferred) or 'peptide'+'source' (fallback, e.g. build/peptides.tsv). Returns {} if absent."""
    origin_map = {}
    if annotation_path and os.path.exists(annotation_path):
        ann = pd.read_csv(annotation_path, sep="\t", dtype=str)
        pcol = next((c for c in ann.columns if c.lower() == "peptide"), None)
        ocol = (next((c for c in ann.columns if c.lower() == "noncanonical_origin"), None)
                or next((c for c in ann.columns if c.lower() == "source"), None))
        if pcol and ocol:
            for _, r in ann.iterrows():
                key = il(strip_peptide(r[pcol]))
                if key:
                    origin_map[key] = str(r[ocol]).strip()
            print(f"[INFO] loaded {len(origin_map)} peptide->origin entries from {annotation_path} "
                  f"(peptide='{pcol}', origin='{ocol}')", file=sys.stderr)
        else:
            print(f"[WARN] {annotation_path} needs 'peptide' + 'noncanonical_origin'/'source'; "
                  f"found {list(ann.columns)}. Origins will be 'unknown'.", file=sys.stderr)
    else:
        print(f"[WARN] annotation '{annotation_path}' not found; origins will be 'unknown'.",
              file=sys.stderr)
    return origin_map


def prepare(df, haplo, *, pep_col, run_col, group_col, annotation, tissue_annotation,
            tumor_token, normal_token, min_len, max_len):
    """Full preparation shared by all predictor arms.
    Returns a collapsed dataframe with columns:
        sample, tissue, pep_bare, category, noncanonical_origin
    keeping (sample, tissue, peptide) as the unit so tumor/normal stays split.
    """
    n_in = len(df)
    sample_tokens = build_sample_tokens(list(haplo))

    # --- sample assignment ---
    a = df[run_col].map(lambda r: assign_sample(r, sample_tokens))
    df = df.copy()
    df["sample"] = [x[0] for x in a]
    df["_amb"] = [x[1] for x in a]
    if df["sample"].isna().any():
        runs = sorted(df.loc[df["sample"].isna(), run_col].dropna().unique())[:10]
        print(f"[WARN] {df['sample'].isna().sum()} PSMs matched NO sample. e.g.:", file=sys.stderr)
        for r in runs:
            print(f"        {r}", file=sys.stderr)
    if df["_amb"].any():
        print(f"[WARN] {int(df['_amb'].sum())} PSMs matched MULTIPLE samples (took first).",
              file=sys.stderr)
    df = df[df["sample"].notna()].copy()

    # --- tissue (per run; sample can be both tumor & normal) ---
    tt, nt = tumor_token.lower(), normal_token.lower()
    def tissue_of(run):
        r = str(run).lower(); ht, hn = tt in r, nt in r
        if ht and not hn: return "tumor"
        if hn and not ht: return "normal"
        if ht and hn:     return "ambiguous"
        return "unknown"
    df["tissue"] = df[run_col].map(tissue_of)
    if tissue_annotation and os.path.exists(tissue_annotation) and os.path.getsize(tissue_annotation) > 0:
        try:
            ta = pd.read_csv(tissue_annotation, sep="\t", dtype=str)
        except Exception as e:
            print(f"[WARN] could not read tissue annotation {tissue_annotation}: {e}", file=sys.stderr)
            ta = pd.DataFrame()
        scol = next((c for c in ta.columns if c.lower() == "sample"), None)
        tcol = next((c for c in ta.columns if c.lower() == "tissue"), None)
        if scol and tcol:
            ov = {str(r[scol]).strip(): str(r[tcol]).strip().lower() for _, r in ta.iterrows()}
            mask = df["tissue"].isin(["unknown", "ambiguous"])
            df.loc[mask, "tissue"] = df.loc[mask, "sample"].map(lambda s: ov.get(s, "unknown"))
            print(f"[INFO] applied tissue overrides for {len(ov)} samples from {tissue_annotation}",
                  file=sys.stderr)
    print(f"[INFO] tissue (per PSM): {df['tissue'].value_counts().to_dict()}", file=sys.stderr)

    # --- clean peptides + validity ---
    df["pep_bare"] = df[pep_col].map(strip_peptide)
    valid = df["pep_bare"].map(lambda p: len(p) >= min_len and len(p) <= max_len and set(p) <= AA)
    if (~valid).any():
        print(f"[INFO] dropping {int((~valid).sum())} peptides outside {min_len}-{max_len}aa / "
              f"non-standard residues.", file=sys.stderr)
    df = df[valid].copy()

    # --- category ---
    if group_col in df.columns:
        df["category"] = df[group_col].map(
            lambda g: "non-canonical" if str(g).upper() == "NONCANON" else "canonical")
    else:
        df["category"] = "unknown"

    # --- collapse to (sample, tissue, peptide) ---
    def pick_cat(s):
        return "non-canonical" if (s == "non-canonical").any() else s.iloc[0]
    uniq = (df.groupby(["sample", "tissue", "pep_bare"], as_index=False)
              .agg(category=("category", pick_cat)))

    # --- origin ---
    origin_map = load_origin_map(annotation)
    def origin_of(row):
        if row["category"] == "canonical":
            return "canonical_origin"
        return origin_map.get(il(row["pep_bare"]), "unknown")
    uniq["noncanonical_origin"] = uniq.apply(origin_of, axis=1)

    print(f"[INFO] {n_in} PSMs -> {len(uniq)} unique (sample, tissue, peptide) rows "
          f"({uniq.groupby('tissue').size().to_dict()} by tissue).", file=sys.stderr)
    return uniq
