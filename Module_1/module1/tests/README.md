# Module 1 tests

Run: `nf-test test` (needs the `nf-test` binary + `nextflow` on PATH; Python deps: biopython, pandas).
All 12 tests execute the real scripts except the two scaffold stubs.

## Tier 1 - real
- **TRANSLATE.nf.test** - rules engine on the 12-record multi-origin fixture: pinned 51 ORF
  peptides; SNV/fusion/viral/microbial/frameshift -> 3 frames, ERV/lncRNA/riboseq -> 6 frames;
  `class` = variant vs noncanonical. Regression: explicit `origin=` overrides header guessing.
- **FASTQ2FASTA_META.nf.test** - origin stamped from metadata, quality stripped, read id kept.
- **VARIANT_EXTRACT.nf.test** - reference window + ALT applied; SNV vs indel class; sample/tissue.
- **VISUALIZE.nf.test** - IGV-style renderer selects a variant + a non-canonical peptide and emits PNGs.

## Tier 2 - scaffold contracts (stub)
- **ALIGN_EXTRACT.nf.test**, **VARIANT_READAWARE.nf.test** - verify each scaffold's IO contract
  (emits an origin/variant-annotated nucleotide FASTA that TRANSLATE can consume). The real bodies
  need an aligner (+ combined DB) / IsoVar and exit 3 until wired.

## Tier 3 - pipeline dispatch + guards
- FASTA -> 51 peptides; FASTQ_META runs; VARIANT propagates `--sample/--tissue` (regression for the
  shell-arg templating fix); guards: `--vcf` without `--reference` fails, `--fastq` alone fails.

Note: `tests/nextflow.config` caps resources to 1 CPU / 3 GB so tests run on small machines.
