#!/usr/bin/env python3
"""Path A / step 2 -- build the Sage search FASTA.

Design (see README 'Path A caveat'):
  target space = { published non-canonical peptides (tagged NONCANON) }
               + { full canonical UniProt proteome (tagged CANON) as realistic background }

The canonical proteome is what keeps target-decoy competition meaningful and prevents the
"I only put in true positives" circularity. Do NOT add decoys here -- Sage generates them.

Each published peptide becomes its own short FASTA entry; under a no-enzyme Sage search this
peptide (and its >=min_len sub-peptides) becomes searchable. Header carries the group tag so
04_group_fdr_recovery.py can split CANON vs NONCANON for group-specific FDR.

Usage:  python 02_build_search_fasta.py config.yaml
"""
import sys, os, yaml
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "common"))
from sequtils import read_fasta, write_fasta, dedup_by_sequence, tag


def published_peptide_records(peptides_tsv, min_len):
    df = pd.read_csv(peptides_tsv, sep="\t")
    df = df[df.peptide.str.len() >= min_len]
    for i, row in df.iterrows():
        grp = row["group"]
        src = (str(row.get("source", "")) or ("noncanon" if grp == "NONCANON" else "canon"))
        src = src.replace(" ", "_")[:24]
        # header: <id>|grp=..|src=..   (tag() enforces the grp bucket; we override src label)
        header = f"chong_pub_{i}|grp={grp}|src=pub_{src}"
        yield header, row["peptide"].strip().upper()


def proteome_records(fasta_path):
    for header, seq in read_fasta(fasta_path):
        acc = header.split()[0]
        yield f"{acc}|grp=CANON|src=uniprot", seq


def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    b = cfg["build"]
    min_len = b.get("min_pep_len", 8)
    out_fasta = b["out_fasta"]
    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)

    # order matters for dedup: published peptides first so a peptide that also exists in the
    # proteome keeps its published provenance (rare, but keeps bookkeeping honest)
    def all_records():
        yield from published_peptide_records("build/peptides.tsv", min_len)
        prot = cfg["inputs"].get("canonical_proteome_fasta")
        if prot and os.path.exists(prot):
            yield from proteome_records(prot)
        else:
            print("  [warn] canonical proteome not found -- FDR will be OPTIMISTIC. "
                  "Add UP000005640 for a meaningful background.", file=sys.stderr)

    n = 0
    def counted():
        nonlocal n
        for rec in dedup_by_sequence(all_records()):
            n += 1
            yield rec

    write_fasta(counted(), out_fasta)
    print(f"Wrote {out_fasta} : {n} entries")
    print("Reminder: point sage_config.json -> database.fasta at this file, and run a NO-ENZYME search.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.yaml")
