Supply your own:
  sage_config_HLAI.json  sage_config_HLAII.json   (no-enzyme; Immuno-HCD rescoring)
  ms2rescore_config.json
  config.yaml            (build.*, inputs.canonical_proteome_fasta: proteome.fasta,
                          search.min_len/max_len, fdr.q_threshold, outputs.fdr_report)
  peptides.tsv           (columns: peptide, group[, source] — 02 hard-codes build/peptides.tsv)
  UP000005640.fasta      (canonical proteome background)
Also: put validation_targets.tsv in bin/ (04 reads it from its own dir).

--- Search-FASTA input (BUILD_FASTA) ------------------------------------------
--peptides_tsv <tsv>     REQUIRED. Peptide table to build the Sage DB from.
--peptides_format        'module1' (default) reads a Module 1 translation TSV
                         (columns: sequence, class, origin, source, ...);
                         'chong' reads the legacy peptide/group/source TSV via
                         02_build_search_fasta.py + a real config.yaml.
--canonical_proteome <f> Optional CANON background (e.g. UniProt UP000005640).
                         Omit -> NONCANON-only DB (FDR will be optimistic).
--min_pep_len N          Min peptide length kept (default 8).
Builder for module1 format: bin/02_build_search_fasta_new.py
