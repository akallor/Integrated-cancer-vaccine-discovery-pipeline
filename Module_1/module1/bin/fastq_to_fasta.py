#!/usr/bin/env python3
"""fastq_to_fasta.py -- Module 1 input path (Option 2): FASTQ + metadata -> nucleotide FASTA.

Preserves the machine read identifier, strips the quality string, and stamps the biological
origin (+ sample/tissue) into the header so the translator can apply the right rule.

Metadata TSV/CSV: at least two columns
    column 1: fastq read id (or key)      column 2: origin  (SNV/Aberrant_splicing/ERV/ncRNA/...)
optionally 'sample' and 'tissue' columns. If exactly one metadata row is given and its key
matches no read, its values are used as a file-wide default (the common "one origin per file"
case, e.g. a fusion-enriched library).

Origin labels are normalized to the translator's classes (e.g. SNV -> variant_snv,
Aberrant_splicing -> aberrant_splicing, ncRNA/lncRNA -> ncrna, ERV -> erv, ...).
"""
import argparse
import os
import sys

from Bio import SeqIO

ORIGIN_ALIAS = {
    "snv": "variant_snv", "variant": "variant_snv", "snp": "variant_snv",
    "indel": "frameshift_indel", "frameshift": "frameshift_indel",
    "fusion": "fusion",
    "aberrant_splicing": "aberrant_splicing", "splicing": "aberrant_splicing",
    "splice": "aberrant_splicing", "neojunction": "aberrant_splicing",
    "erv": "erv", "herv": "erv",
    "ncrna": "ncrna", "lncrna": "ncrna", "lincrna": "ncrna", "ncorf": "ncorf",
    "riboseq": "ncorf", "viral": "viral", "virus": "viral",
    "microbial": "microbial", "bacterial": "microbial",
    "intronic": "intronic", "intergenic": "intergenic", "out_of_frame": "out_of_frame",
}


def norm_origin(v):
    if v is None:
        return "unknown"
    return ORIGIN_ALIAS.get(str(v).strip().lower(), str(v).strip().lower() or "unknown")


def load_metadata(path):
    """Return (per_read: {id: dict(origin,sample,tissue)}, default: dict|None)."""
    import csv
    per_read, rows = {}, []
    with open(path) as fh:
        sniff = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in sniff else ","
        reader = csv.reader(fh, delimiter=delim)
        data = [r for r in reader if r and not r[0].lstrip().startswith("#")]
    if not data:
        return {}, None
    # detect a header row
    head = [c.strip().lower() for c in data[0]]
    has_header = any(h in ("fastq_header", "header", "read", "read_id", "readid", "id", "origin")
                     for h in head)
    if has_header:
        idx = {name: head.index(name) for name in head}
        key_i = next((idx[k] for k in ("fastq_header", "header", "read", "read_id", "readid", "id")
                      if k in idx), 0)
        org_i = idx.get("origin", 1)
        smp_i = idx.get("sample")
        tis_i = idx.get("tissue")
        body = data[1:]
    else:
        key_i, org_i, smp_i, tis_i = 0, 1, None, None
        body = data
    for r in body:
        if len(r) <= max(key_i, org_i):
            continue
        key = r[key_i].strip().lstrip("@").split()[0] if r[key_i].strip() else ""
        rec = dict(origin=norm_origin(r[org_i]),
                   sample=(r[smp_i].strip() if smp_i is not None and len(r) > smp_i else None),
                   tissue=(r[tis_i].strip() if tis_i is not None and len(r) > tis_i else None))
        rows.append((key, rec))
        if key:
            per_read[key] = rec
    default = rows[0][1] if len(rows) == 1 else None
    return per_read, default


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True, help="output FASTA")
    ap.add_argument("--sample", default=None, help="fallback sample")
    ap.add_argument("--tissue", default=None, help="fallback tissue")
    ap.add_argument("--origin", default=None, help="fallback origin if a read is not in metadata")
    args = ap.parse_args()

    per_read, default = load_metadata(args.metadata)
    fallback_origin = norm_origin(args.origin) if args.origin else (
        default["origin"] if default else "unknown")

    n = 0
    with open(args.out, "w") as out:
        for rec in SeqIO.parse(args.fastq, "fastq"):
            n += 1
            meta = per_read.get(rec.id) or per_read.get(rec.id.split()[0]) or default or {}
            origin = meta.get("origin") or fallback_origin
            sample = meta.get("sample") or args.sample or "NA"
            tissue = meta.get("tissue") or args.tissue or "unknown"
            out.write(f">{rec.id} origin={origin} sample={sample} tissue={tissue} "
                      f"readid={rec.id}\n{str(rec.seq)}\n")
    print(f"[fastq2fasta] {n} reads -> {args.out} (quality stripped; origin stamped)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
