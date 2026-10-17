# Module_2 test suite (nf-test)

Three tiers, mirroring the pipeline's testing strategy.

    nf-test test                          # everything
    nf-test test tests/GROUP_FDR.nf.test  # one file

Requires: nf-test, Nextflow (>=24.04), and — for the two REAL unit tests — python
with pandas + pyyaml on PATH (the tasks run 02/04 directly).

## Tier 1 — unit (real execution, no external tools)
- BUILD_FASTA.nf.test  — runs 02_build_search_fasta.py: CANON+NONCANON assembly,
  min_pep_len filter, and the dedup-keeps-published-provenance regression.
  (uses common/sequtils.py — reference impl shipped; swap in yours.)
- GROUP_FDR.nf.test    — runs 04_group_fdr_recovery.py: per-group CANON/NONCANON FDR,
  positive-control recovery through flanks + I/L, and the SKIP-vs-MISS length-window
  regression. (uses bin/validation_targets.tsv — reference list; replace with your 36.)

## Tier 2 — process contract (stub run; -stub)
DOWNLOAD_PRIDE / CONVERT_STREAM / SAGE / MS2RESCORE — validate IO contracts, the
per-class tuple shapes, prefix pass-through, and summary emission without pridepy,
ThermoRawFileParser, Sage, ms2rescore, or AWS. MS2RESCORE also pins the post-rescoring
PSM/peptide count (your explainability/robustness number).

## Tier 3 — end-to-end (stub run; -stub)
Module_2.pipeline.nf.test — the whole six-stage DAG with --hla_class both, asserting
both HLA-I and HLA-II branches complete through to group FDR.

Note: fixtures live in tests/fixtures/. Stub tests point every path input at fixtures,
so they are self-contained regardless of what you have placed in assets/.
