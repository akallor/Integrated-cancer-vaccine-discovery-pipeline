#!/usr/bin/env python3
"""translate_orfs.py -- rules-based nucleotide->peptide translation (Module 1 core).

For each record in a nucleotide FASTA:
  1. classify origin from the header (assets/origin_patterns.tsv; explicit origin= wins),
  2. look up the translation rule (assets/translation_rules.tsv): #frames (3 or 6), strand,
     NCBI genetic-code table, min ORF length, whether a start Met is required,
  3. translate the selected frames, split each frame on stop codons into ORF fragments,
  4. emit one protein record per surviving ORF, with full provenance in the header.

Multiple origins in one FASTA are handled per-record (each gets its own rule).

Outputs:
  <out>.fasta : protein sequences, header = key=value provenance tokens
  <out>.tsv   : one row per ORF; sequence column + provenance columns (feeds Module 2's DB)

Frames: 0-2 = forward offsets, 3-5 = reverse-complement offsets. orf_offset_nt is the ORF's
start offset (nt) within the translated strand, so an ORF is traceable back to the record.
"""
import argparse
import os
import re
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m1_common as C

TSV_COLS = ["id", "sequence", "class", "origin", "source", "sample", "tissue",
            "locus", "strand", "frame", "genetic_code", "orf_offset_nt", "aa_len",
            "variant", "gene", "readid"]


def _clean_nt(s):
    return re.sub(r"[^ACGTN]", "", str(s).upper().replace("U", "T"))


def iter_orfs(nt_seq, rule, header_strand):
    """Yield (frame_index, strand_char, orf_offset_nt, peptide) per surviving ORF."""
    s = _clean_nt(nt_seq)
    use_rev = (rule["strand"] == "auto" and header_strand == "-")
    if rule["frames"] == 6:
        plan = [(0, "+"), (1, "+"), (2, "+"), (0, "-"), (1, "-"), (2, "-")]
    elif rule["frames"] == 1:
        plan = [(0, "-" if use_rev else "+")]
    elif use_rev:
        plan = [(0, "-"), (1, "-"), (2, "-")]
    else:
        plan = [(0, "+"), (1, "+"), (2, "+")]

    fwd = Seq(s)
    rev = fwd.reverse_complement()

    # ---- defined-ORF (CDS) mode: single frame, initiator -> M, terminal stop dropped ----
    if rule.get("cds") and rule["frames"] == 1:
        strand = plan[0][1]
        src = fwd if strand == "+" else rev
        sub = str(src)
        sub = sub[: len(sub) // 3 * 3]
        if sub:
            orf = str(Seq(sub).translate(table=rule["genetic_code"])).split("*")[0]
            if orf:
                orf = "M" + orf[1:]        # initiator Met (covers non-AUG starts: TTG/CTG/GTG...)
            if len(orf) >= rule["min_aa_len"]:
                yield 0, strand, 0, orf
        return

    for fi, (off, strand) in enumerate(plan):
        src = fwd if strand == "+" else rev
        sub = str(src)[off:]
        sub = sub[: len(sub) // 3 * 3]
        if not sub:
            continue
        aa = str(Seq(sub).translate(table=rule["genetic_code"]))
        cursor = 0  # aa index within this frame's translation
        for frag in aa.split("*"):
            frag_start_aa = cursor
            cursor += len(frag) + 1  # +1 for the stop that was split out
            pep = frag
            if rule["require_start"]:
                mi = pep.find("M")
                pep = pep[mi:] if mi >= 0 else ""
            if len(pep) >= rule["min_aa_len"]:
                offset_in_pep = (len(frag) - len(pep))  # residues trimmed before the Met
                orf_offset_nt = off + 3 * (frag_start_aa + offset_in_pep)
                # report offset in forward coordinates for '-' frames too, for traceability
                yield fi, strand, orf_offset_nt, pep


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="input nucleotide FASTA")
    ap.add_argument("--out", required=True, help="output prefix (writes <out>.fasta and <out>.tsv)")
    ap.add_argument("--patterns", default=None, help="origin_patterns.tsv (default: alongside script assets)")
    ap.add_argument("--rules", default=None, help="translation_rules.tsv")
    ap.add_argument("--sample", default=None, help="fallback sample name when header lacks one")
    ap.add_argument("--tissue", default=None, help="fallback tissue (tumor/normal) when header lacks one")
    ap.add_argument("--min-aa-len", type=int, default=None,
                    help="override every rule's min ORF length (global floor)")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    patt_path = args.patterns or os.path.join(here, "..", "assets", "origin_patterns.tsv")
    rule_path = args.rules or os.path.join(here, "..", "assets", "translation_rules.tsv")
    patterns = C.load_patterns(patt_path)
    rules = C.load_rules(rule_path)
    if args.min_aa_len is not None:
        for r in rules.values():
            r["min_aa_len"] = args.min_aa_len
    known = set(rules.keys())

    n_in = n_orf = 0
    origin_counts = {}
    fa_out = open(f"{args.out}.fasta", "w")
    rows = []

    for rec in SeqIO.parse(args.fasta, "fasta"):
        n_in += 1
        header = rec.description
        origin = C.classify_header(header, patterns, known)
        origin_counts[origin] = origin_counts.get(origin, 0) + 1
        fields = C.parse_fields(header)
        rule = rules.get(origin, rules["unknown"])
        cls = C.emit_class(origin)
        sample = fields.get("sample") or args.sample or "NA"
        tissue = fields.get("tissue") or args.tissue or "unknown"
        hstrand = fields.get("strand", "+")
        locus = (f"{fields.get('chrom')}:{fields.get('start')}-{fields.get('end')}"
                 if fields.get("chrom") else fields.get("coords", "NA"))

        k = 0
        for fi, strand, orf_nt, pep in iter_orfs(rec.seq, rule, hstrand):
            k += 1
            n_orf += 1
            new_id = f"{fields.get('readid', rec.id)}|orf{k}"
            row = dict(
                id=new_id, sequence=pep, **{"class": cls}, origin=origin,
                source=fields.get("source", "NA"), sample=sample, tissue=tissue,
                locus=locus, strand=strand, frame=fi, genetic_code=rule["genetic_code"],
                orf_offset_nt=orf_nt, aa_len=len(pep),
                variant=fields.get("variant", "NA"), gene=fields.get("gene", "NA"),
                readid=fields.get("readid", rec.id),
            )
            rows.append(row)
            desc = " ".join(f"{c}={row[c]}" for c in TSV_COLS if c not in ("id", "sequence"))
            fa_out.write(f">{new_id} {desc}\n{pep}\n")
    fa_out.close()

    pd.DataFrame(rows, columns=TSV_COLS).to_csv(f"{args.out}.tsv", sep="\t", index=False)
    print(f"[translate] {n_in} input records -> {n_orf} ORF peptides "
          f"(>= min length). Origins: {origin_counts}", file=sys.stderr)
    print(f"[translate] wrote {args.out}.fasta and {args.out}.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
