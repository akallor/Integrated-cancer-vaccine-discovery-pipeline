#!/usr/bin/env python3
"""extract_variant_seqs.py -- Module 1 variant path: VCF + reference -> mutant nucleotide FASTA.

For each VCF record, extracts a reference window around the variant, applies the ALT allele,
and emits a nucleotide FASTA record whose header carries the variant, gene, coordinates, and
sample/tissue -- so the translator applies the variant (3-frame) rule and keeps provenance.

WHY A REFERENCE (not a raw FASTQ): a variant is a genomic coordinate + allele; reconstructing
the sequence *around* that coordinate needs the reference base context (or reads already aligned
to it). A raw FASTQ has no coordinates. So this path takes --reference. The read-aware
alternative (align the FASTQ, then translate reads spanning the variant, e.g. via IsoVar) is a
separate arm wired in the Nextflow module and flagged there.

SNV/MNV -> origin=variant_snv ; length-changing indel -> origin=frameshift_indel.
Only the first ALT allele is used per record (split multiallelic upstream if needed).
"""
import argparse
import os
import re
import sys

from Bio import SeqIO


def parse_info_gene(info):
    for key in ("GENE", "SYMBOL", "Gene"):
        m = re.search(rf"(?:^|;){key}=([^;]+)", info)
        if m:
            return m.group(1)
    m = re.search(r"ANN=[^|]*\|[^|]*\|[^|]*\|([^|]+)\|", info)  # SnpEff ANN gene field
    return m.group(1) if m else "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--reference", required=True, help="reference genome FASTA (indexed lazily)")
    ap.add_argument("--out", required=True, help="output nucleotide FASTA")
    ap.add_argument("--flank", type=int, default=30, help="bp of reference context each side (default 30)")
    ap.add_argument("--sample", default=None)
    ap.add_argument("--tissue", default=None)
    args = ap.parse_args()

    ref = SeqIO.index(args.reference, "fasta")
    sample = args.sample or "NA"
    tissue = args.tissue or "unknown"

    n_var = n_out = 0
    with open(args.out, "w") as out:
        for ln in open(args.vcf):
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            n_var += 1
            chrom, pos, vid, refa, alta = f[0], int(f[1]), f[2], f[3].upper(), f[4].split(",")[0].upper()
            info = f[7] if len(f) > 7 else ""
            gene = parse_info_gene(info)
            if chrom not in ref:
                print(f"[WARN] {chrom} not in reference; skipping {chrom}:{pos}", file=sys.stderr)
                continue
            seq = str(ref[chrom].seq)
            vstart = pos - 1                      # 0-based
            if seq[vstart:vstart + len(refa)].upper() != refa:
                print(f"[WARN] REF mismatch at {chrom}:{pos} (vcf={refa}, "
                      f"ref={seq[vstart:vstart+len(refa)]}); applying ALT anyway", file=sys.stderr)
            wstart = max(0, vstart - args.flank)
            wend = min(len(seq), vstart + len(refa) + args.flank)
            left = seq[wstart:vstart]
            right = seq[vstart + len(refa):wend]
            mut = left + alta + right
            origin = "variant_snv" if len(refa) == len(alta) else "frameshift_indel"
            vtype = "SNV" if len(refa) == len(alta) == 1 else ("INDEL" if len(refa) != len(alta) else "MNV")
            vstr = f"{chrom}:{pos}:{refa}:{alta}"
            header = (f"{vid if vid not in ('.', '') else vstr} origin={origin} "
                      f"variant={vstr} gene={gene} coords={chrom}:{wstart+1}-{wend} strand=+ "
                      f"type={vtype} sample={sample} tissue={tissue}")
            out.write(f">{header}\n{mut}\n")
            n_out += 1
    print(f"[variant-extract] {n_var} VCF records -> {n_out} mutant windows (flank={args.flank})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
