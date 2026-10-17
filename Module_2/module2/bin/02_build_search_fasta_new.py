#!/usr/bin/env python3
"""Build the Sage search FASTA from a Module 1 translation TSV.

Module 1 emits one row per translated ORF peptide with columns:
    id, sequence, class, origin, source, sample, tissue, locus, strand, frame,
    genetic_code, orf_offset_nt, aa_len, variant, gene, readid
(the original 02_build_search_fasta.py read the Chong-paper layout
 peptide/group/source/... — this reader replaces that.)

Design is unchanged from the original:
    target space = { Module 1 peptides (tagged NONCANON) }
                 + { canonical proteome (tagged CANON) as a realistic background }
Each peptide becomes its own short FASTA entry; under a no-enzyme Sage search the peptide
(and its >=min_len sub-peptides) becomes searchable. Headers carry |grp=..|src=.. so
04_group_fdr_recovery.py can split CANON vs NONCANON for group-specific FDR. No decoys here —
Sage generates them. The canonical proteome is what keeps target-decoy competition honest;
omitting it makes FDR optimistic (a warning is printed).

Usage:
    02_build_search_fasta_new.py --peptides peptides.tsv --out build/search.fasta \\
                                 [--proteome canonical.fasta] [--min-len 8]
    # or, for parity with the original: --config config.yaml  (fills any unset paths/values)
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from sequtils import read_fasta, write_fasta, dedup_by_sequence  # noqa: E402

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_id(s):
    """Header uses '|' as the field separator, so strip '|'/'>'/whitespace from the id."""
    return re.sub(r"\s+", "_", str(s)).replace("|", "_").replace(">", "").strip("_")


def _col(df, *names):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n in low:
            return low[n]
    return None


def module1_records(peptides_tsv, min_len):
    import pandas as pd
    df = pd.read_csv(peptides_tsv, sep="\t", dtype=str).fillna("")
    seqcol = _col(df, "sequence", "peptide")
    if seqcol is None:
        sys.exit(f"ERROR: no 'sequence' (or 'peptide') column in {peptides_tsv} "
                 f"— columns are {list(df.columns)}. Is this a Module 1 TSV?")
    idcol  = _col(df, "id", "readid")
    clscol = _col(df, "class")
    orgcol = _col(df, "origin")
    srccol = _col(df, "source")
    n_in = n_out = n_bad = 0
    for i, row in df.iterrows():
        n_in += 1
        pep = str(row[seqcol]).strip().upper()
        if len(pep) < min_len or not set(pep) <= VALID_AA:
            n_bad += 1
            continue
        # Module 1 is the non-canonical translation module: everything is NONCANON unless a
        # row is explicitly marked canonical.
        cls = str(row[clscol]).strip().lower() if clscol else ""
        grp = "CANON" if cls == "canonical" else "NONCANON"
        # provenance label: prefer the rule-based origin (ncorf/variant_snv/...), else biotype
        src = ""
        for c in (orgcol, srccol):
            if c and str(row[c]).strip() and str(row[c]).strip().upper() not in ("NA", "NAN"):
                src = str(row[c]).strip()
                break
        src = (src or "noncanon").replace(" ", "_")[:24]
        base = clean_id(row[idcol]) if idcol else f"m1_{i}"
        yield f"{base}|grp={grp}|src={src}", pep
        n_out += 1
    print(f"  Module 1 peptides: {n_in} rows -> {n_out} FASTA entries "
          f"(dropped {n_bad} shorter than {min_len} aa or with non-standard residues)",
          file=sys.stderr)


def proteome_records(fasta_path):
    for header, seq in read_fasta(fasta_path):
        acc = header.split()[0]
        yield f"{acc}|grp=CANON|src=uniprot", seq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peptides", help="Module 1 translation TSV")
    ap.add_argument("--out", help="output search FASTA")
    ap.add_argument("--proteome", default=None, help="canonical proteome FASTA (CANON background)")
    ap.add_argument("--min-len", type=int, default=8)
    ap.add_argument("--config", default=None, help="optional config.yaml (fills unset values)")
    args = ap.parse_args()

    # optional config.yaml compatibility (same keys as the original)
    if args.config and os.path.exists(args.config) and os.path.getsize(args.config) > 0:
        import yaml
        cfg = yaml.safe_load(open(args.config)) or {}
        b, ins = cfg.get("build", {}), cfg.get("inputs", {})
        args.out = args.out or b.get("out_fasta")
        args.min_len = args.min_len or b.get("min_pep_len", 8)
        args.proteome = args.proteome or ins.get("canonical_proteome_fasta")
        args.peptides = args.peptides or ins.get("peptides_tsv")

    if not args.peptides or not args.out:
        ap.error("need --peptides and --out (or a --config that provides them)")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    def all_records():
        # Module 1 peptides first so a peptide also present in the proteome keeps NONCANON provenance
        yield from module1_records(args.peptides, args.min_len)
        if args.proteome and os.path.exists(args.proteome) and os.path.getsize(args.proteome) > 0:
            yield from proteome_records(args.proteome)
        else:
            print("  [warn] canonical proteome not provided/empty — FDR will be OPTIMISTIC. "
                  "Pass --proteome UP000005640.fasta for a meaningful CANON background.",
                  file=sys.stderr)

    n = 0
    def counted():
        nonlocal n
        for rec in dedup_by_sequence(all_records()):
            n += 1
            yield rec

    write_fasta(counted(), args.out)
    print(f"Wrote {args.out} : {n} entries")
    print("Reminder: point sage_config.json -> database.fasta at this file, and run a NO-ENZYME search.")


if __name__ == "__main__":
    main()
