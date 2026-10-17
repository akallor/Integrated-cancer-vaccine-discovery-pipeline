#!/bin/bash -ue
set -euo pipefail
base=$(basename "test5.fasta" | sed 's/^ntseq_//; s/\.[^.]*$//')
S="DEMO"; SARG=""; [ -n "$S" ] && SARG="--sample $S"
T="tumor"; TARG=""; [ -n "$T" ] && TARG="--tissue $T"
M=""; MIN=""; [ -n "$M" ] && MIN="--min-aa-len $M"
python3 "$(command -v translate_orfs.py)"         --fasta "test5.fasta" --patterns "origin_patterns.tsv" --rules "translation_rules.tsv"         $SARG $TARG $MIN --out "peptides_${base}"
NP=$(( $(wc -l < "peptides_${base}.tsv") - 1 ))
echo "[translate] ${base}: ${NP} ORF peptides" >&2
