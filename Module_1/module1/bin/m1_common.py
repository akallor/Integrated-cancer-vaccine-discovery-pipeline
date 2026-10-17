#!/usr/bin/env python3
"""m1_common.py -- shared helpers for Module 1 (nucleotide -> peptide translation).

Origin classification is data-driven:
  * assets/origin_patterns.tsv    header regex  -> origin_class   (extensible, ordered)
  * assets/translation_rules.tsv  origin_class  -> frames/strand/code/filters (extensible)

An explicit `origin=<class>` token in a header always wins over pattern guessing, so an
upstream step that already knows the origin (metadata / alignment) can state it directly.
"""
import os
import re
import sys

import pandas as pd

# origins treated as "variant" (vs "noncanonical") in the emitted class= field
VARIANT_CLASSES = {"variant_snv", "frameshift_indel"}


def load_patterns(path):
    """[(priority, origin_class, compiled_regex)] sorted by priority. '#'/blank lines skipped."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            pr, cls, rx = ln.split("\t", 2)
            rows.append((int(pr), cls.strip(), re.compile(rx, re.IGNORECASE)))
    rows.sort(key=lambda x: x[0])
    return rows


def load_rules(path):
    """origin_class -> dict(frames, strand, genetic_code, min_aa_len, require_start, cds).

    Optional 7th column `cds` (true/false): treat the record as an already-defined ORF -> a single
    forward frame, initiator codon rendered as M (non-AUG starts included, per Ribo-seq/GENCODE
    convention), terminal stop dropped. Use for annotated ORFs (e.g. Ribo-seq ncORFs) to keep the
    search DB lean; leave false for discovery classes that need multi-frame translation."""
    rules = {}
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if p[0].strip().lower() == "origin_class" or len(p) < 6:
                continue
            rules[p[0].strip()] = dict(
                frames=int(p[1]), strand=p[2].strip(), genetic_code=int(p[3]),
                min_aa_len=int(p[4]), require_start=p[5].strip().lower() == "true",
                cds=(len(p) > 6 and p[6].strip().lower() == "true"))
    if "unknown" not in rules:  # always need a fallback
        rules["unknown"] = dict(frames=6, strand="both", genetic_code=1,
                                min_aa_len=8, require_start=False, cds=False)
    return rules


def classify_header(header, patterns, known_classes):
    """Return the origin_class for a header. Explicit origin=<class> wins; else first regex match."""
    m = re.search(r"origin=([A-Za-z0-9_]+)", header, re.IGNORECASE)
    if m and m.group(1).lower() in known_classes:
        return m.group(1).lower()
    for _pr, cls, rx in patterns:
        if cls == "__explicit__":
            continue
        if rx.search(header):
            return cls
    return "unknown"


def parse_fields(header):
    """Best-effort extraction of provenance fields from a heterogeneous header.
    Returns a dict; missing keys are simply absent. Never raises on odd headers."""
    d = {}
    toks = header.split()
    d["readid"] = toks[0] if toks else ""

    # explicit key=value tokens take precedence
    for key in ("sample", "tissue", "gene", "strand", "source", "coords"):
        m = re.search(rf"\b{key}=([^\s|]+)", header, re.IGNORECASE)
        if m:
            d[key] = m.group(1)

    # coordinates: chr:start-end  or  chr:pos(:ref:alt)
    if "coords" not in d:
        m = re.search(r"(chr[\w.]+):(\d+)(?:-(\d+))?", header, re.IGNORECASE)
        if m:
            d["chrom"], d["start"] = m.group(1), m.group(2)
            d["end"] = m.group(3) or m.group(2)
    else:
        m = re.search(r"(chr[\w.]+):(\d+)(?:-(\d+))?", d["coords"], re.IGNORECASE)
        if m:
            d["chrom"], d["start"], d["end"] = m.group(1), m.group(2), (m.group(3) or m.group(2))

    # strand from _strand=+ style if not already set
    if "strand" not in d:
        m = re.search(r"_strand=([+-])", header) or re.search(r"strand[=:]([+-])", header)
        if m:
            d["strand"] = m.group(1)

    # tissue token if not explicit
    if "tissue" not in d:
        if re.search(r"tumou?r", header, re.IGNORECASE):
            d["tissue"] = "tumor"
        elif re.search(r"\bnormal\b", header, re.IGNORECASE):
            d["tissue"] = "normal"

    # sample from a "Source: X" description if not explicit
    if "sample" not in d:
        m = re.search(r"Source:\s*([^\s|]+)", header)
        if m:
            d["sample"] = m.group(1)

    # variant string (chr:pos:ref:alt or rsID or HGVS)
    m = (re.search(r"(chr[\w.]+:\d+:[ACGTN]+:[ACGTN]+)", header, re.IGNORECASE)
         or re.search(r"\b(rs\d+)\b", header, re.IGNORECASE)
         or re.search(r"(_c\.\S+)", header))
    if m:
        d["variant"] = m.group(1)

    # source subtype / biotype if not explicit
    if "source" not in d:
        m = re.search(r"(lncRNA|lincRNA|ncRNA|HERV\w*|ERV\w*|LTR|processed_transcript|antisense|"
                      r"protein_coding)", header, re.IGNORECASE)
        if m:
            d["source"] = m.group(1)
    return d


def emit_class(origin_class):
    """Coarse class= label used in output headers/TSV."""
    if origin_class in VARIANT_CLASSES:
        return "variant"
    if origin_class == "unknown":
        return "unknown"
    return "noncanonical"
