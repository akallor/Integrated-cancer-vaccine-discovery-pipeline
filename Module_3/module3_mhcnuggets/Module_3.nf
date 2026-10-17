#!/usr/bin/env nextflow
// ===============================================================================
//  Module_3.nf  --  Binding affinity / presentation (Module 3)
// -------------------------------------------------------------------------------
//  Dispatch on --class {I|II} and --predictor {MHCFlurry|MHCnuggets}:
//    --class I  --predictor MHCFlurry   -> 05_mhcflurry_predict.py       (+ SHAP explain + plots)
//    --class I  --predictor MHCnuggets  -> 05_b_mhcnuggets_predict.py
//    --class II --predictor MHCnuggets  -> 05_c_mhcnuggets_predict_hlaii.py (common-allele panel)
//    --class II --predictor MHCFlurry   -> REJECTED (no class-II MHCflurry model)
//
//  Real script interfaces (from the uploaded code):
//    predictors : --peptides <f> --haplotypes <Sample/Haplotype tsv> [--annotation build/peptides.tsv] --out <tsv>
//                 (class-II adds --alleles <panel> --min-len --max-len)
//    explain    : run_explain.py --input <pred.tsv> --outdir <dir> --backend {mhcflurry|surrogate}
//    plots      : plot_allele_jaccard.py <pred.tsv> --out-prefix <p>       (real CLI, raw predictions)
//                 plot_binders_counts.py / plot_binders_proportions.py     (read df_binder_raw_counts.tsv)
//                 plot_allele_proportions_v2.py                            (reads df_allele_pep_counts.tsv)
//                 plot_combined_hla_data_v2.py                             (reads df_origin_grp.tsv + df_tissue_grp.tsv)
//
//  DEPENDENCIES YOU MUST SUPPLY IN bin/ (not uploaded):
//    mod3_common.py   -- imported by BOTH MHCnuggets predictors (05_b, 05_c)
//  PROVIDED (inferred, replace if you have your own):
//    mod3_aggregate.py -- builds the 4 fixed-name summary tables the static plots read
//
//  NOTES baked in from the real code:
//    * plot_allele_proportions_v2.py has a blank first line (shebang on line 2), so all
//      python scripts are invoked via `python3 $(command -v <script>)` -- robust to that,
//      and keeps sibling imports (mod3_common / mhcflurry_explain) resolvable.
//    * static plot scripts call plt.show(); processes set MPLBACKEND=Agg so that's a no-op.
//    * run_explain.py --backend surrogate runs the full SHAP/ISM/motif stack WITHOUT
//      MHCflurry/TensorFlow -- use it for CI (needs `shap`).
//
//  RUN:
//    nextflow run Module_3.nf --class I --predictor MHCFlurry \
//      --peptides 'results_module2/06_group_fdr/**/group_fdr_recovery*.tsv' \
//      --haplotypes assets/PXD013649_hlaI_haplotypes.tsv --annotation assets/peptides.tsv
//
//    # re-run only the diagnostic plots on an existing prediction table:
//    nextflow run Module_3.nf -entry PLOTS --predictions results_module3/predictions/pred_x.tsv
// ===============================================================================

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------- params
params.peptides         = null        // glob or file (1..N) of identified peptides
params.haplotypes       = null        // Sample<TAB>Haplotype table (per-sample genotype)
params.annotation       = null        // build/peptides.tsv (group/origin); optional
params.panel            = null        // HLA-II common-allele panel (class II only; optional)
params.tissue_annotation = null
params.predictions      = null        // for -entry PLOTS: existing prediction table(s)

params.hla_class        = null        // alias for --class
params.predictor        = null        // MHCFlurry | MHCnuggets
params.run_plots        = null        // default: on for both predictors; set false to skip diagnostic plots
params.explain_backend  = 'mhcflurry' // or 'surrogate' (no MHCflurry needed)
params.explain_seed     = 7
params.min_len_ii       = 9
params.max_len_ii       = 30
params.outdir           = 'results_module3'
params.dry_run          = false       // predictors: verify sample/tissue/peptide assignment, no model

// script names (params so they can be corrected without touching processes)
params.mhcflurry_predict_script     = '05_mhcflurry_predict.py'
params.mhcnuggets_predict_i_script  = '05_b_mhcnuggets_predict.py'
params.mhcnuggets_predict_ii_script = '05_c_mhcnuggets_predict_hlaii.py'
params.explain_script               = 'run_explain.py'
params.aggregate_script             = 'mod3_aggregate.py'
params.plot_jaccard_script          = 'plot_allele_jaccard.py'
params.plot_combined_script         = 'plot_combined_hla_data_v2.py'
params.plot_binder_counts_script    = 'plot_binders_counts.py'
params.plot_binder_prop_script      = 'plot_binders_proportions.py'
params.plot_allele_prop_script      = 'plot_allele_proportions_v2.py'

// ---------------------------------------------------------------------------- progress bar
def renderBar(int n, int total) {
    int w = 30
    int f   = ((w * n) / total + 0.5) as int
    int pct = ((100 * n) / total + 0.5) as int
    return "${pct}%|" + ('█' * f) + ('░' * (w - f)) + "| ${n}/${total}"
}
def stageBar(int n, int total, String label, String detail) {
    println ''
    println "  ${label}"
    println "  ${renderBar(n, total)}"
    if (detail?.trim()) println "  \u21B3 ${detail}"
}

// ============================================================================ processes

process PREDICT_MHCFLURRY {
    tag "${peptides.simpleName}"
    publishDir "${params.outdir}/predictions", mode: 'copy', pattern: 'pred_*.tsv'
    publishDir "${params.outdir}/predictions", mode: 'copy', pattern: '*.summary'
    input:
    path peptides
    path haplotypes
    path annotation
    output:
    path 'pred_*.tsv', emit: pred
    path '*.summary',  emit: sm
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{peptides}" | sed 's/\\.[^.]*$//')
    echo "[predict] MHCflurry HLA-I on ${base} ..."
    ANN=""; [ "!{annotation}" != "NO_ANNOTATION" ] && ANN="--annotation !{annotation}"
    SCRIPT="$(command -v !{params.mhcflurry_predict_script})"
    if [ "!{params.dry_run}" = "true" ]; then
        python3 "$SCRIPT" --peptides "!{peptides}" --haplotypes "!{haplotypes}" $ANN --dry-run 2>&1 | tee "pred_${base}.tsv"
        echo "DRY_RUN MHCFlurry_HLA-I [script=!{params.mhcflurry_predict_script}] (${base})" | tee "pred_${base}.summary"
    else
        python3 "$SCRIPT" --peptides "!{peptides}" --haplotypes "!{haplotypes}" $ANN --out "pred_${base}.tsv"
        NP=$(( $(wc -l < "pred_${base}.tsv") - 1 ))
        echo "MHCFlurry_HLA-I [script=!{params.mhcflurry_predict_script}]: ${NP} rows (${base})" | tee "pred_${base}.summary"
    fi
    '''
    stub:
    """
    base=\$(basename "${peptides}" | sed 's/\\.[^.]*\$//')
    printf 'peptide\\tsample\\ttissue\\tcategory\\tnoncanonical_origin\\tbinding_allele\\tmhcflurry_affinity\\tmhcflurry_affinity_percentile\\tmhcflurry_binder_call\\n' > pred_\${base}.tsv
    printf 'KYKDRTNILF\\tS1\\ttumor\\tnon-canonical\\tlncRNA\\tHLA-A*02:01\\t45.2\\t0.30\\tSB\\n' >> pred_\${base}.tsv
    printf 'SLLQHLIGL\\tS1\\ttumor\\tcanonical\\t\\tHLA-A*02:01\\t780\\t1.50\\tWB\\n'            >> pred_\${base}.tsv
    printf 'GILGFVFTL\\tS1\\tnormal\\tcanonical\\t\\tHLA-B*15:01\\t22000\\t45.0\\tNB\\n'        >> pred_\${base}.tsv
    echo "MHCFlurry_HLA-I [script=${params.mhcflurry_predict_script}]: 3 rows (\${base})" > pred_\${base}.summary
    """
}

process PREDICT_MHCNUGGETS_I {
    tag "${peptides.simpleName}"
    publishDir "${params.outdir}/predictions", mode: 'copy', pattern: 'pred_*.tsv'
    publishDir "${params.outdir}/predictions", mode: 'copy', pattern: '*.summary'
    input:
    path peptides
    path haplotypes
    path annotation
    output:
    path 'pred_*.tsv', emit: pred
    path '*.summary',  emit: sm
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{peptides}" | sed 's/\\.[^.]*$//')
    echo "[predict] MHCnuggets HLA-I on ${base} ..."
    ANN=""; [ "!{annotation}" != "NO_ANNOTATION" ] && ANN="--annotation !{annotation}"
    SCRIPT="$(command -v !{params.mhcnuggets_predict_i_script})"
    if [ "!{params.dry_run}" = "true" ]; then
        python3 "$SCRIPT" --peptides "!{peptides}" --haplotypes "!{haplotypes}" $ANN --dry-run 2>&1 | tee "pred_${base}.tsv"
        echo "DRY_RUN MHCnuggets_HLA-I [script=!{params.mhcnuggets_predict_i_script}] (${base})" | tee "pred_${base}.summary"
    else
        python3 "$SCRIPT" --peptides "!{peptides}" --haplotypes "!{haplotypes}" $ANN --out "pred_${base}.tsv"
        NP=$(( $(wc -l < "pred_${base}.tsv") - 1 ))
        echo "MHCnuggets_HLA-I [script=!{params.mhcnuggets_predict_i_script}]: ${NP} rows (${base})" | tee "pred_${base}.summary"
    fi
    '''
    stub:
    """
    base=\$(basename "${peptides}" | sed 's/\\.[^.]*\$//')
    printf 'peptide\\tsample\\ttissue\\tcategory\\tnoncanonical_origin\\tbinding_allele\\tmhcnuggets_ic50\\tmhcnuggets_binder_call\\n' > pred_\${base}.tsv
    printf 'KYKDRTNILF\\tS1\\ttumor\\tnon-canonical\\tlncRNA\\tHLA-A02:01\\t42.0\\tSB\\n' >> pred_\${base}.tsv
    printf 'SLLQHLIGL\\tS1\\ttumor\\tcanonical\\t\\tHLA-A02:01\\t320.0\\tWB\\n'           >> pred_\${base}.tsv
    echo "MHCnuggets_HLA-I [script=${params.mhcnuggets_predict_i_script}]: 2 rows (\${base})" > pred_\${base}.summary
    """
}

process PREDICT_MHCNUGGETS_II {
    tag "${peptides.simpleName}"
    publishDir "${params.outdir}/predictions", mode: 'copy', pattern: 'pred_*.tsv'
    publishDir "${params.outdir}/predictions", mode: 'copy', pattern: '*.summary'
    input:
    path peptides
    path haplotypes
    path annotation
    path panel
    output:
    path 'pred_*.tsv', emit: pred
    path '*.summary',  emit: sm
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{peptides}" | sed 's/\\.[^.]*$//')
    echo "[predict] MHCnuggets HLA-II on ${base} (common-allele panel) ..."
    ANN="";   [ "!{annotation}" != "NO_ANNOTATION" ] && ANN="--annotation !{annotation}"
    PANEL=""; [ "!{panel}"      != "NO_PANEL" ] && PANEL="--alleles !{panel}"
    SCRIPT="$(command -v !{params.mhcnuggets_predict_ii_script})"
    LEN="--min-len !{params.min_len_ii} --max-len !{params.max_len_ii}"
    if [ "!{params.dry_run}" = "true" ]; then
        python3 "$SCRIPT" --peptides "!{peptides}" --haplotypes "!{haplotypes}" $ANN $PANEL $LEN --dry-run 2>&1 | tee "pred_${base}.tsv"
        echo "DRY_RUN MHCnuggets_HLA-II [script=!{params.mhcnuggets_predict_ii_script}] (${base})" | tee "pred_${base}.summary"
    else
        python3 "$SCRIPT" --peptides "!{peptides}" --haplotypes "!{haplotypes}" $ANN $PANEL $LEN --out "pred_${base}.tsv"
        NP=$(( $(wc -l < "pred_${base}.tsv") - 1 ))
        echo "MHCnuggets_HLA-II [script=!{params.mhcnuggets_predict_ii_script}]: ${NP} rows (${base})" | tee "pred_${base}.summary"
    fi
    '''
    stub:
    """
    base=\$(basename "${peptides}" | sed 's/\\.[^.]*\$//')
    printf 'peptide\\tsample\\ttissue\\tcategory\\tnoncanonical_origin\\tbinding_allele\\tbinding_allele_source\\tmhcnuggets_ic50\\tmhcnuggets_binder_call\\n' > pred_\${base}.tsv
    printf 'NLVPMVATVQGQ\\tS1\\ttumor\\tcanonical\\t\\tHLA-DRB1*15:01\\tcommon_panel_NOT_sample_specific\\t120.0\\tWB\\n' >> pred_\${base}.tsv
    echo "MHCnuggets_HLA-II [script=${params.mhcnuggets_predict_ii_script}]: 1 rows (\${base})" > pred_\${base}.summary
    """
}

process EXPLAIN {
    tag "${predictions.simpleName}"
    publishDir "${params.outdir}/explain", mode: 'copy'
    input:
    path predictions
    output:
    path 'explain_*', emit: expl
    path '*.summary', emit: sm
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{predictions}" | sed 's/^pred_//; s/\\.[^.]*$//')
    echo "[explain] SHAP (backend=!{params.explain_backend}) for ${base} ..."
    python3 "$(command -v !{params.explain_script})" \
        --input "!{predictions}" --outdir "explain_${base}" \
        --backend "!{params.explain_backend}" --seed !{params.explain_seed}
    NF=$(find "explain_${base}" -type f | wc -l)
    echo "SHAP explanation: ${base} (${NF} artefacts: positional SHAP + ISM weight matrices + motifs)" \
        | tee "explain_${base}.summary"
    '''
    stub:
    """
    base=\$(basename "${predictions}" | sed 's/^pred_//; s/\\.[^.]*\$//')
    mkdir -p explain_\${base}
    : > explain_\${base}/positional_shap.tsv
    : > explain_\${base}/positional_shap_grid.png
    : > explain_\${base}/ism_weights_grid.png
    echo "SHAP explanation: \${base} (3 artefacts: positional SHAP + ISM weight matrices + motifs)" > explain_\${base}.summary
    """
}

process COMBINE_PREDICTIONS {
    tag "combine"
    publishDir "${params.outdir}", mode: 'copy', pattern: 'combined_predictions.tsv'
    input:
    path 'pred_*.tsv'
    output:
    path 'combined_predictions.tsv', emit: combined
    shell:
    '''
    set -euo pipefail
    awk 'FNR==1 && NR!=1 {next} {print}' pred_*.tsv > combined_predictions.tsv
    '''
    stub:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' pred_*.tsv > combined_predictions.tsv
    """
}

process AGGREGATE {
    tag "aggregate"
    publishDir "${params.outdir}/summary_tables", mode: 'copy', pattern: 'df_*.tsv'
    input:
    path combined
    output:
    path 'df_*.tsv', emit: tables
    shell:
    '''
    set -euo pipefail
    python3 "$(command -v !{params.aggregate_script})" --input "!{combined}" --outdir .
    '''
    stub:
    """
    printf 'sample\\tSB\\tWB\\tNB\\nS1\\t1\\t1\\t1\\n'                 > df_binder_raw_counts.tsv
    printf 'sample\\tpeptide\\tHLA-A*02:01\\nS1\\t3\\t3\\n'           > df_allele_pep_counts.tsv
    printf 'sample\\tcanonical\\tnon_canonical\\nS1\\t2\\t1\\n'       > df_origin_grp.tsv
    printf 'tissue\\tcanonical\\tnon_canonical\\ntumor\\t2\\t1\\n'    > df_tissue_grp.tsv
    """
}

process PLOT_JACCARD {
    tag "jaccard"
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:
    path combined
    output:
    path 'jaccard_*', emit: fig
    path 'jaccard.summary', emit: sm
    shell:
    '''
    set -euo pipefail
    export MPLBACKEND=Agg
    python3 "$(command -v !{params.plot_jaccard_script})" "!{combined}" --out-prefix jaccard_
    NIMG=$(ls jaccard_*.png 2>/dev/null | wc -l)
    echo "Plot jaccard: ${NIMG} heatmap(s) + shared-peptide tables" | tee jaccard.summary
    '''
    stub:
    """
    : > jaccard_allele_jaccard_heatmap.png
    : > jaccard_shared_peptides_by_allele.tsv
    echo "Plot jaccard: 1 heatmap(s) + shared-peptide tables" > jaccard.summary
    """
}

process PLOT_STATIC {
    tag "${plot_name}"
    publishDir "${params.outdir}/plots", mode: 'copy'
    input:
    tuple val(plot_name), val(script)
    path agg
    output:
    path "fig_${plot_name}",       emit: fig
    path "${plot_name}.summary",   emit: sm
    shell:
    '''
    set -euo pipefail
    export MPLBACKEND=Agg
    echo "[plot] !{plot_name} ..."
    # the df_*.tsv summary tables are staged in CWD under their fixed names; run the script here
    python3 "$(command -v !{script})" || { echo "[warn] !{plot_name} produced no output"; }
    mkdir -p "fig_!{plot_name}"
    shopt -s nullglob
    PRED=$(echo "!{params.predictor}" | tr 'A-Z' 'a-z')
    for f in *.png *.pdf; do mv "$f" "fig_!{plot_name}/$(echo "$f" | sed "s/mhcflurry/${PRED}/")"; done
    NF=$(find "fig_!{plot_name}" -type f | wc -l)
    echo "Plot !{plot_name}: ${NF} figure(s)" | tee "!{plot_name}.summary"
    '''
    stub:
    """
    mkdir -p fig_${plot_name}
    : > fig_${plot_name}/${plot_name}.png
    echo "Plot ${plot_name}: 1 figure (stub)" > ${plot_name}.summary
    """
}

// ============================================================================ plot sub-workflow (reusable)
workflow PLOT_ARM {
    take:
    predictions_ch
    main:
    comb = COMBINE_PREDICTIONS(predictions_ch.collect())
    agg  = AGGREGATE(comb.combined)
    jac  = PLOT_JACCARD(comb.combined)
    static_ch = Channel.fromList([
        ['combined_hla',       params.plot_combined_script],
        ['binder_counts',      params.plot_binder_counts_script],
        ['binder_proportions', params.plot_binder_prop_script],
        ['allele_proportions', params.plot_allele_prop_script],
    ])
    st = PLOT_STATIC(static_ch, agg.tables)
    emit:
    plots_sm = jac.sm.mix(st.sm)
}

// entry to re-run only the diagnostic plots on an existing prediction table
workflow PLOTS {
    pin = Channel.fromPath(params.predictions ?: params.peptides, checkIfExists: true)
    PLOT_ARM(pin)
}

// ============================================================================ main
workflow {

    def hlaClass  = (params.containsKey('class') ? params['class'] : params.hla_class)?.toString()
    def predictor = params.predictor?.toString()

    if (!params.peptides)   { error "Missing --peptides (glob or file of identified peptides)" }
    if (!params.haplotypes) { error "Missing --haplotypes (Sample<TAB>Haplotype table)" }
    if (!hlaClass)  { error "Missing --class (I | II)" }
    if (!predictor) { error "Missing --predictor (MHCFlurry | MHCnuggets)" }
    if (!(hlaClass in ['I', 'II'])) { error "--class must be I or II (got '${hlaClass}')" }
    def predLc = predictor.toLowerCase()
    if (!(predLc in ['mhcflurry', 'mhcnuggets'])) { error "--predictor must be MHCFlurry or MHCnuggets" }
    if (hlaClass == 'II' && predLc == 'mhcflurry') {
        error "MHCflurry has no HLA class-II model -- use --predictor MHCnuggets for --class II"
    }

    def flurry   = (predLc == 'mhcflurry')
    def runPlots = (params.run_plots != null) ? (params.run_plots as boolean) : true   // plots for both predictors

    println """
    ===========================================================================
     Module_3  |  Binding affinity / presentation
     class     : ${hlaClass}
     predictor : ${predictor}
     explain   : ${flurry ? "yes (SHAP, backend=${params.explain_backend})" : 'no'}
     plots     : ${runPlots ? 'yes' : 'no'}
    ===========================================================================
    """.stripIndent()

    int total      = 1 + (flurry ? 1 : 0) + (runPlots ? 1 : 0)
    int idxPredict = 1
    int idxExplain = flurry ? 2 : 0
    int idxPlots   = total

    peptides_ch = Channel.fromPath(params.peptides, checkIfExists: true)
    hap = file(params.haplotypes, checkIfExists: true)
    ann = params.annotation ? file(params.annotation, checkIfExists: true) : file("${projectDir}/assets/NO_ANNOTATION")
    pan = params.panel      ? file(params.panel, checkIfExists: true)      : file("${projectDir}/assets/NO_PANEL")

    if (flurry) {
        pr = PREDICT_MHCFLURRY(peptides_ch, hap, ann)
    } else if (hlaClass == 'I') {
        pr = PREDICT_MHCNUGGETS_I(peptides_ch, hap, ann)
    } else {
        pr = PREDICT_MHCNUGGETS_II(peptides_ch, hap, ann, pan)
    }
    predictions = pr.pred

    pr.sm.collect().subscribe { l ->
        stageBar(idxPredict, total, "Predicting binding affinity (HLA-${hlaClass}, ${predictor})", "${l.size()} file(s)")
    }

    if (flurry && !params.dry_run) {
        ex = EXPLAIN(predictions)
        ex.sm.collect().subscribe { l -> stageBar(idxExplain, total, 'Explaining predictions (SHAP)', "${l.size()} explained") }
    }

    if (runPlots && !params.dry_run) {
        arm = PLOT_ARM(predictions)
        arm.plots_sm.collect().subscribe { l -> stageBar(idxPlots, total, 'Diagnostic plots', "${l.size()} plot jobs") }
    }
    workflow.onComplete = {
        if (workflow.success) println "  ${renderBar(total, total)}"
        println ''
        println (workflow.success
            ? "  Module_3 complete  ->  ${params.outdir}/  (status: OK)"
            : "  Module_3 FAILED  (see .nextflow.log)")
    }
}
