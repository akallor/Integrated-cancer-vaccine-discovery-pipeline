Inputs (pass on the CLI):
  --peptides    glob or file (1..N) of identified peptides (Module 2 group-FDR tables).
                Fans out to one prediction task per file.
  --haplotypes  Sample<TAB>Haplotype table (per-sample genotype; for class II used only
                for the sample list).
  --annotation  build/peptides.tsv (peptide -> group/noncanonical_origin). Optional; the
                predictors fall back if omitted.
  --panel       class II only: common HLA-II allele panel (file, one allele/line, or comma
                list). Optional; the class-II script defaults to a top-10 DRB1 panel.

NO_ANNOTATION / NO_PANEL are empty sentinels used to represent "optional input not given"
(distinct names avoid an input-filename collision in the class-II process). Do not remove.
