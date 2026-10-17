#!/bin/bash -ue
set -euo pipefail
if [ "python" = "gviz" ]; then
    Rscript "$(command -v plot_peptide_gviz.R)" --tsv *.tsv --outdir . --n 2 --seed 12
else
    python3 "$(command -v plot_peptide_igv.py)" --tsv *.tsv --outdir . --n 2 --seed 12
fi
