# Module_3 test suite (nf-test)

    nf-test test
    nf-test test tests/Module_3.pipeline.nf.test

## Tier 1 — unit, REAL execution (no MHCflurry/MHCnuggets needed)
- COMBINE_PREDICTIONS — generic TSV concat, single header.
- AGGREGATE — the four summary tables, counts pinned to the fixture.
- PLOT_JACCARD — real CLI on raw predictions -> heatmaps + shared-peptide tables.
- PLOT_STATIC — a static plot (reads a pre-aggregated df_*.tsv); a second case pins the
  broken-shebang `plot_allele_proportions_v2.py` running via python3.
- EXPLAIN — `run_explain.py --backend surrogate` runs the full SHAP/ISM/motif stack.
  **Requires `shap` in the runner env** (`pip install shap`); no MHCflurry/TensorFlow.

## Tier 2 — contract, stub
- PREDICT_MHCFLURRY, PREDICT_MHCNUGGETS_II — each has (a) a stub contract test for output
  schema + dispatched script name, and (b) a REAL `--dry_run` test that runs the actual
  predictor script's sample-matching / tissue / category prep (mod3_common for the MHCnuggets
  arm) up to — but not including — the model call. Needs pandas only; no MHCflurry/MHCnuggets.
  Full model runs still need `mhcflurry`/`mhcnuggets`; add a smoke tier once installed.

## Tier 3 — end-to-end
- Three valid class x predictor routes + the `class II + MHCFlurry` rejection (stub).
- `PLOTS` entry runs the whole plot arm FOR REAL on an existing prediction table
  (real aggregate + all five plot scripts; needs matplotlib/pandas only).

Tests use `hla_class` (alias for --class) to dodge the Groovy `class` keyword in the
params closure; the workflow reads --class and hla_class identically.
