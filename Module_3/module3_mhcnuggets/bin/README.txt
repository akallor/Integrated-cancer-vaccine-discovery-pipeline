bin/ — scripts on PATH (Nextflow adds bin/ automatically). All are invoked as
`python3 $(command -v <script>)` by the pipeline, which (a) is robust to the blank
first line in plot_allele_proportions_v2.py, and (b) keeps sibling imports resolvable.

UPLOADED (real):
  05_mhcflurry_predict.py            HLA-I MHCflurry     (self-contained; needs `mhcflurry`)
  05_b_mhcnuggets_predict.py         HLA-I MHCnuggets    (imports mod3_common; needs `mhcnuggets`)
  05_c_mhcnuggets_predict_hlaii.py   HLA-II MHCnuggets   (imports mod3_common; common-allele PANEL,
                                                          not sample typing; needs `mhcnuggets`)
  run_explain.py + mhcflurry_explain.py   SHAP/ISM/motif explainer (HLA-I MHCflurry only).
                                     --backend surrogate runs the FULL stack without MHCflurry
                                     (needs `shap`; `logomaker` optional -> logo, else heatmap).
  plot_allele_jaccard.py             real CLI: <predictions.tsv> --out-prefix <p>
  plot_binders_counts.py             reads df_binder_raw_counts.tsv   (no CLI, plt.show())
  plot_binders_proportions.py        reads df_binder_raw_counts.tsv   (no CLI, plt.show())
  plot_allele_proportions_v2.py      reads df_allele_pep_counts.tsv   (no CLI; BLANK line 1)
  plot_combined_hla_data_v2.py       reads df_origin_grp.tsv + df_tissue_grp.tsv (no CLI)

INCLUDED:
  mod3_common.py   -- shared prep (load_haplotypes, prepare) imported by 05_b and 05_c.
                      The class-I MHCflurry script has its own copy of that logic.

PROVIDED HERE (inferred, replace with your own if you have one):
  mod3_aggregate.py -- builds the four fixed-name summary tables the static plots read,
                       straight from a combined prediction table. Predictor-agnostic
                       (auto-detects the *_binder_call column), so it also works on
                       MHCnuggets output even though the plot titles/filenames say "mhcflurry".

Predictor CLIs (real):
  05_mhcflurry_predict.py     --peptides <f> --haplotypes <Sample/Haplotype tsv> [--annotation <f>] --out <tsv>
  05_b_mhcnuggets_predict.py  --peptides <f> --haplotypes <f> [--annotation <f>] --out <tsv>
  05_c_...hlaii.py            --peptides <f> --haplotypes <f> [--annotation <f>] [--alleles <panel>]
                              --min-len 9 --max-len 30 --out <tsv>
