#!/usr/bin/env python3
"""Path A / step 4 -- group-specific FDR + recovery of known Chong peptides.

Why group-specific FDR (the NewAnce idea): the NONCANON group has a far worse true:false
prior than CANON, so a single global 1% FDR silently lets through many false non-canonical
PSMs. We split PSMs by their FASTA header group tag (grp=CANON / grp=NONCANON) and compute
target-decoy q-values SEPARATELY within each group.

Then we check whether the known-presented Chong peptides (validation_targets.tsv) are
recovered at the chosen q-threshold. That recovery -- not raw PSM counts -- is the success
metric for the Path A validation.

Input: a rescored PSM table (MS2Rescore/mokapot output). Column names are auto-detected but
can be overridden below. Needs, per PSM: a score, a target/decoy label, the peptide, and the
protein/accession string (which carries the |grp= tag).

Usage:  python 04_group_fdr_recovery.py <rescored_psms.tsv> [config.yaml]
"""
import sys, os, re, yaml
import pandas as pd

def pick(df, *cands):
    low = {c.lower(): c for c in df.columns}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    # substring match fallback
    for c in cands:
        for col in df.columns:
            if c.lower() in col.lower():
                return col
    return None

def qvalues(sub, score_col, decoy_col, higher_better=True):
    """Standard target-decoy q-values within one group."""
    s = sub.sort_values(score_col, ascending=not higher_better).copy()
    is_decoy = s[decoy_col].astype(bool).values
    cum_decoy = 0; cum_target = 0; fdr = []
    for d in is_decoy:
        if d: cum_decoy += 1
        else: cum_target += 1
        fdr.append(cum_decoy / max(cum_target, 1))
    # monotonize -> q-value
    q = [0.0]*len(fdr); m = 1.0
    for i in range(len(fdr)-1, -1, -1):
        m = min(m, fdr[i]); q[i] = m
    s["qvalue"] = q
    return s

def strip_peptide(p):
    """Reduce a peptide string to bare residues: drop flanks (K.PEP.L), mods ([+15.99]/(ox)),
    charge suffixes, and whitespace. Used both for length checks and target matching."""
    s = str(p).strip().upper()
    # drop flanking residues like  K.PEPTIDE.L  -> PEPTIDE
    if s.count(".") >= 2:
        parts = s.split(".")
        s = ".".join(parts[1:-1]) if len(parts) > 2 else s
    # remove modification annotations in [..] or (..) and any non-residue chars
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r"\([^)]*\)", "", s)
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    return s

def detect_length_window(df, pep_col, cfg):
    """Determine the [min_len, max_len] the search actually used.
    Priority: config.search.min_len/max_len -> observed target-peptide lengths in the file."""
    s = cfg.get("search", {}) if isinstance(cfg, dict) else {}
    lo, hi = s.get("min_len"), s.get("max_len")
    if lo and hi:
        return int(lo), int(hi), "config"
    lens = df[pep_col].map(lambda p: len(strip_peptide(p)))
    lens = lens[lens > 0]
    if len(lens):
        return int(lens.min()), int(lens.max()), "observed in file"
    return 8, 30, "fallback default"

def main(psm_path, cfg_path="config.yaml"):
    cfg = yaml.safe_load(open(cfg_path)) if os.path.exists(cfg_path) else {}
    qthr = cfg.get("fdr", {}).get("q_threshold", 0.01)

    df = pd.read_csv(psm_path, sep=None, engine="python")
    pep_col   = pick(df, "peptide", "sequence", "pep")
    score_col = pick(df, "mokapot score", "score", "rescore", "hyperscore", "spectral_angle")
    decoy_col = pick(df, "is_decoy", "label", "decoy", "td")
    prot_col  = pick(df, "proteins", "protein", "accession", "proteinids", "protein_list")
    q_col     = pick(df, "mokapot q-value", "q-value", "qvalue")   # use tool q if present
    assert pep_col and prot_col, f"could not find peptide/protein columns in {list(df.columns)}"

    # normalize decoy label to bool (mokapot: label = 1 target / -1 decoy, or is_decoy bool)
    if decoy_col is None:
        # infer from decoy_tag in protein string
        df["_decoy"] = df[prot_col].astype(str).str.contains("rev_|decoy", case=False)
        decoy_col = "_decoy"
    else:
        v = df[decoy_col]
        df["_decoy"] = v.astype(str).str.lower().isin(["true","1","-1","decoy","yes"]) | (v == -1)
        decoy_col = "_decoy"

    # group from |grp= tag in protein string; a PSM mapping to any NONCANON protein is NONCANON
    def grp_of(prot):
        return "NONCANON" if re.search(r'grp=NONCANON', str(prot)) else "CANON"
    df["group"] = df[prot_col].map(grp_of)
    # normalize peptides: strip flanks/mods, then I->L (MS can't distinguish isoleucine/leucine)
    df["peptide_norm"] = df[pep_col].map(strip_peptide).str.replace("I", "L", regex=False)

    # what length window did this search actually use? (drives searchable-vs-total reporting)
    min_len, max_len, src = detect_length_window(df, pep_col, cfg)

    frames = []
    for g, sub in df.groupby("group"):
        if q_col:  # trust the tool's q if it computed one, but still per-group filter
            sub = sub.rename(columns={q_col: "qvalue"})
        else:
            sub = qvalues(sub, score_col, decoy_col)
        keep = sub[(~sub["_decoy"]) & (sub["qvalue"] <= qthr)]
        print(f"[{g}] PSMs={len(sub)}  passing q<= {qthr}: {len(keep)}  "
              f"unique peptides: {keep['peptide_norm'].nunique()}")
        frames.append(sub.assign(group=g))
    allq = pd.concat(frames, ignore_index=True)
    passing = allq[(~allq["_decoy"]) & (allq["qvalue"] <= qthr)]

    # ---- recovery of known Chong peptides ----
    tgt_path = os.path.join(os.path.dirname(__file__), "validation_targets.tsv")
    tgt = pd.read_csv(tgt_path, sep="\t")
    tgt["peptide_bare"] = tgt.peptide.map(strip_peptide)
    tgt["peptide_norm"] = tgt.peptide_bare.str.replace("I", "L", regex=False)
    tgt["length"] = tgt.peptide_bare.str.len()
    # a target is "searchable" only if its length falls within THIS search's window
    tgt["searchable"] = (tgt.length >= min_len) & (tgt.length <= max_len)
    found = set(passing["peptide_norm"])
    # only searchable targets can be 'recovered' by this search; short/long ones are SKIP, not MISS
    tgt["recovered"] = tgt.searchable & tgt.peptide_norm.isin(found)

    print(f"\n=== Search length window: {min_len}-{max_len} aa ({src}) ===")
    print("=== Recovery of known Chong peptides (positive controls) ===")
    for _, r in tgt.iterrows():
        if not r.searchable:
            status = "SKIP  "   # outside the search's length window -> cannot be identified here
            why = f"(len {r.length}, outside {min_len}-{max_len})"
        elif r.recovered:
            status = "FOUND "; why = ""
        else:
            status = "MISS  "; why = ""
        print(f"  {status} {r.peptide:14s} {r['class']:12s} {r['note']} {why}")

    n_total     = len(tgt)
    n_search    = int(tgt.searchable.sum())
    n_skip      = n_total - n_search
    n_recovered = int(tgt.recovered.sum())
    print(f"\nRecovered {n_recovered}/{n_search} SEARCHABLE targets "
          f"(length {min_len}-{max_len}) at group-q<= {qthr}")
    print(f"  {n_skip}/{n_total} targets were outside the search length window and could NOT be "
          f"identified by this search (not a failure -- wrong length class).")
    print(f"  Overall: {n_recovered}/{n_total} of all listed targets.")
    if n_search == 0:
        print("  NOTE: 0 targets were length-eligible for this search. This target list does not "
              "validate this search's peptide-length class (e.g. class-I targets vs an HLA-II search).")

    out = cfg.get("outputs", {}).get("fdr_report", "results/group_fdr_recovery.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    passing.to_csv(out, sep="\t", index=False)
    tgt.to_csv(out.replace(".tsv", ".recovery.tsv"), sep="\t", index=False)
    print(f"\nWrote {out} and recovery table.")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "config.yaml")
