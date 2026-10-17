#!/usr/bin/env nextflow
/*
===============================================================================
 Module_2.nf  —  MS/MS validation (Module 2 of the neoantigen pipeline)
-------------------------------------------------------------------------------
 Chains six existing scripts into one runnable DSL2 pipeline:

   [1] download_pride_data_with_pridepy_arg.sh   PRIDE  -> local -> S3 (.raw)
   [2] convert_and_stream_raw_files_arg.sh       .raw   -> .mzML (streamed via S3)
   [3] 02_build_search_fasta.py                  build the Sage search FASTA
   [4] run_sage.sh  (parameterised)              Sage database search
   [5] run_ms2rescore.sh (parameterised)         MS2Rescore rescoring
   [6] 04_group_fdr_recovery.py                  group-specific FDR + recovery

 Data plane = S3 (consistent with the streaming design of your scripts).
 Nextflow orchestrates; heavy data lives in S3 and is pulled only when needed.

-------------------------------------------------------------------------------
 EXPECTED LAYOUT (place your assets here before running):

   Module_2.nf
   nextflow.config
   bin/
     download_pride_data_with_pridepy_arg.sh   (chmod +x)
     convert_and_stream_raw_files_arg.sh       (chmod +x)
     02_build_search_fasta.py                  (chmod +x)
     04_group_fdr_recovery.py                  (chmod +x)
     validation_targets.tsv        <-- 04 reads this from its OWN dir (bin/)
   common/
     sequtils.py                   <-- imported by 02 as ../common/sequtils
   assets/
     sage_config_HLAI.json  sage_config_HLAII.json
     ms2rescore_config.json
     config.yaml                   <-- for 02 (build.*) and 04 (fdr.*, outputs.*)
     peptides.tsv                  <-- 02 hard-codes build/peptides.tsv
     UP000005640.fasta             <-- canonical proteome background

 config.yaml MUST contain (minimal):
   build:
     out_fasta: build/search.fasta          # relative -> emitted by BUILD_FASTA
     min_pep_len: 8
   inputs:
     canonical_proteome_fasta: proteome.fasta   # staged name (see BUILD_FASTA)
   search:
     min_len: 8
     max_len: 15                             # 8-15 for HLA-I; widen for HLA-II
   fdr:
     q_threshold: 0.01
   outputs:
     fdr_report: results/group_fdr_recovery.tsv

-------------------------------------------------------------------------------
 RUN (single file or many — count is discovered from PRIDE/S3, not hard-coded):

   nextflow run Module_2.nf \
     --pride    PXD013649 \
     --s3_raw   s3://my-bucket/PXD013649/raw \
     --s3_mzml  s3://my-bucket/PXD013649/mzml \
     --hla_class I            # 'I' | 'II' | 'both'

 Data sources (choose ONE):
   --pride  PXD013649                                  PRIDE accession/URL
   --repo   https://repository.jpostdb.org/entry/JPST003393
                                                       ANY PX-family repo: PRIDE (PXD/RPXD),
                                                       jPOST (JPST), MassIVE (MSV), iProX (IPX),
                                                       PDC (PDC) accession OR entry URL. Resolved
                                                       via pridepy (jPOST/MassIVE/iProX native),
                                                       so jPOST entries whose PRIDE mirror is empty
                                                       still work. Empty repo -> the run fails with
                                                       "No files exist in the provided repository link".
   --local  /data/inhouse/raw   (or)  s3://bucket/inhouse
                                                       in-house data: a local directory is uploaded to
                                                       --s3_raw; an s3:// prefix is used in place (no
                                                       download, no --s3_raw needed).
   --from_mzml --s3_mzml s3://b/m                      skip download+convert; reuse existing mzML

 Vendor / extension (arm c): --vendor thermo|bruker|sciex picks .raw/.d/.wiff;
   override explicitly with --ext .raw|.d|.wiff (EXT wins over VENDOR).

 Validate wiring without executing tools:
   nextflow run Module_2.nf -preview \
     --repo JPST003393 --s3_raw s3://b/r --s3_mzml s3://b/m --peptides_tsv p.tsv
===============================================================================
*/

nextflow.enable.dsl = 2

// ----------------------------------------------------------------------------
// Parameters
// ----------------------------------------------------------------------------
params.pride    = null                       // PRIDE accession or full URL (pridepy path)
params.repo     = null                       // ANY repository: PXD/JPST/MSV accession or entry URL (ProteomeXchange-resolved)
params.local    = null                       // in-house data: local directory OR s3:// prefix (no public download)
params.ext      = null                       // optional raw extension override (.raw/.d/.wiff); else derived from --vendor
params.s3_raw   = null                        // s3://bucket/.../raw
params.s3_mzml  = null                        // s3://bucket/.../mzml
params.pride_files = null                     // comma-sep filenames -> fetch only these (single-file demo)
params.from_mzml   = false                    // skip download+convert; start from an existing --s3_mzml
params.vendor      = 'thermo'                  // 'thermo' (.raw) | 'bruker' (.d) | 'sciex' (.wiff)
params.hla_class = 'I'                         // 'I' | 'II' | 'both'

params.sage_config_I  = "${projectDir}/assets/sage_config_HLAI.json"
params.sage_config_II = "${projectDir}/assets/sage_config_HLAII.json"
params.batch_size_I   = 4
params.batch_size_II  = 2
params.mzml_I  = null                          // optional per-class mzML prefix (default s3_mzml)
params.mzml_II = null

params.build_config       = null                       // legacy Chong config.yaml (only for --peptides_format chong)
params.peptides_tsv       = null                       // Module 1 translation TSV (REQUIRED for BUILD_FASTA)
params.canonical_proteome = null                       // optional CANON background (UniProt); omit -> optimistic FDR
params.peptides_format    = 'module1'                  // 'module1' (Module 1 TSV) | 'chong' (peptide/group/source TSV)
params.min_pep_len        = 8
params.ms2rescore_config  = "${projectDir}/assets/ms2rescore_config.json"
params.fdr_config         = "${projectDir}/assets/config.yaml"

params.min_ms2         = 1                      // exported for convert QC
params.rescore_threads = 4
params.rescore_use_deeplc = true   // DeepLC RT features; auto-falls back to MS²PIP-only if calibration fails
params.id_diag_script = '07_id_diagnostics.py'
params.id_fdr         = 0.01       // FDR for the pre/post-rescoring ID diagnostics
params.run_id_plots   = true       // Module 2 identification diagnostics (set false to skip)
params.spectrum_script = '08_spectrum_plots.py'
params.spectrum_seed   = 0
params.run_spectrum_plots = true   // annotated MS2 spectrum + XIC for 1 canonical & 1 non-canonical peptide
params.outdir          = "results_module2"

// ----------------------------------------------------------------------------
// tqdm-style stage progress bar (advances live on the console as stages finish)
// ----------------------------------------------------------------------------
def renderBar(int n) {
    int w = 30
    int total = 6
    int f   = ((w * n) / total + 0.5) as int
    int pct = ((100 * n) / total + 0.5) as int
    return "${pct}%|" + ('█' * f) + ('░' * (w - f)) + "| ${n}/${total}"
}

def stageBar(int n, String label, String detail) {
    println ''
    println "  ${label}"
    println "  ${renderBar(n)}"
    if (detail?.trim()) println "  \u21B3 ${detail}"
}

// ----------------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------------

process DOWNLOAD_PRIDE {
    tag { "${pride}" }
    publishDir "${params.outdir}/01_download", mode: 'copy', pattern: '*.summary'

    input:
    val pride
    val s3_raw
    val pride_files

    output:
    val s3_raw,               emit: s3_raw
    path 'download.summary',  emit: sm

    shell:
    '''
    set -euo pipefail
    echo "[1/6] Downloading files from PRIDE (!{pride})!{pride_files ? ' [files: '+pride_files+']' : ''} ..."

    export VENDOR="!{params.vendor}"
    export EXT="!{params.ext ?: ''}"
    if [ "$VENDOR" = "bruker" ]; then PAT='\\.d(\\.zip)?$'; else PAT='\\.raw$'; fi
    download_pride_data_with_pridepy_arg.sh "!{pride}" "!{s3_raw}" "!{pride_files ?: ''}"

    # ---- machine-readable summary from S3 (vendor-aware pattern) ----
    S3="!{s3_raw}"; S3="${S3%/}"
    aws s3 ls "${S3}/" --recursive --summarize \
      | awk -v pat="$PAT" '$0 ~ pat {c++} /Total Size:/{b=$3} \
             END{gb=b/1073741824; \
                 printf "Files downloaded from PRIDE: Total number %d and total size %.2f GB\\n", c, gb}' \
      | tee download.summary
    '''
    stub:
    """
    echo "Files downloaded from PRIDE: Total number 1 and total size 1.50 GB" > download.summary
    """

}

// ---- generic repository download (ProteomeXchange-resolved: PRIDE / JPOST / MassIVE / iProX) ----
process DOWNLOAD_REPO {
    tag { "${repo}" }
    publishDir "${params.outdir}/01_download", mode: 'copy', pattern: '*.summary'
    input:
    val repo
    val s3_raw
    val files
    output:
    val s3_raw,              emit: s3_raw
    path 'download.summary', emit: sm
    shell:
    '''
    set -euo pipefail
    echo "[1/6] Downloading from repository (!{repo}) ..."
    export VENDOR="!{params.vendor}"
    export EXT="!{params.ext ?: ''}"
    download_repo.sh "!{repo}" "!{s3_raw}" "!{files ?: ''}" | tee download.summary
    '''
    stub:
    """
    echo "Files downloaded from repository: Total number 1 and total size 0.90 GB" > download.summary
    """
}

// ---- in-house data: ingest from a local directory into the S3 raw prefix (no public download) ----
process LOCAL_INGEST {
    tag { "${local}" }
    publishDir "${params.outdir}/01_download", mode: 'copy', pattern: '*.summary'
    input:
    val local
    val s3_raw
    output:
    val s3_raw,             emit: s3_raw
    path 'ingest.summary',  emit: sm
    shell:
    '''
    set -euo pipefail
    echo "[1/6] Ingesting in-house data from !{local} -> !{s3_raw} ..."
    EXT="!{params.ext ?: ''}"
    if [ -z "$EXT" ]; then case "!{params.vendor}" in bruker) EXT=".d";; sciex) EXT=".wiff";; *) EXT=".raw";; esac; fi
    [ -d "!{local}" ] || { echo "ERROR: local path !{local} not found (use an s3:// prefix or a real directory)" >&2; exit 1; }
    N=$(find "!{local}" -maxdepth 3 -name "*${EXT}*" 2>/dev/null | wc -l)
    if [ "$N" -eq 0 ]; then echo "No files with extension ${EXT} in !{local}" >&2; exit 1; fi
    aws s3 sync "!{local}" "!{s3_raw}" --exclude "*" --include "*${EXT}" --include "*${EXT}/*" --include "*${EXT}.zip"
    echo "Ingested ${N} local file(s) (${EXT}) -> !{s3_raw}" | tee ingest.summary
    '''
    stub:
    """
    echo "Ingested 1 file -> !{s3_raw}" > ingest.summary
    """
}

process CONVERT_STREAM {
    tag { "${s3_mzml}" }
    publishDir "${params.outdir}/02_convert", mode: 'copy', pattern: '*.summary'

    input:
    val s3_raw
    val s3_mzml

    output:
    val s3_mzml,             emit: s3_mzml
    path 'convert.summary',  emit: sm

    shell:
    '''
    set -euo pipefail
    echo "[2/6] Converting RAW -> mzML (streamed via S3) ..."

    export MIN_MS2="!{params.min_ms2}"
    export VENDOR="!{params.vendor}"
    export EXT="!{params.ext ?: ''}"
    convert_and_stream_raw_files_arg.sh "!{s3_raw}" "!{s3_mzml}"

    S3="!{s3_mzml}"; S3="${S3%/}"
    aws s3 ls "${S3}/" --recursive --summarize \
      | awk '/\\.mzML$/{c++} /Total Size:/{b=$3} \
             END{gb=b/1073741824; \
                 printf "RAW->mzML conversion complete: %d mzML files, %.2f GB\\n", c, gb}' \
      | tee convert.summary
    '''
    stub:
    """
    echo "RAW->mzML conversion complete: 3 mzML files, 4.20 GB" > convert.summary
    """

}

process BUILD_FASTA {
    tag "search-space"
    publishDir "${params.outdir}/03_fasta", mode: 'copy', pattern: 'build/search.fasta'
    publishDir "${params.outdir}/03_fasta", mode: 'copy', pattern: '*.summary'

    input:
    val  ready                                   // gate: run after CONVERT_STREAM
    path build_cfg, stageAs: 'config.yaml'
    path peptides
    path proteome,  stageAs: 'proteome.fasta'

    output:
    path 'build/search.fasta',   emit: fasta
    path 'buildfasta.summary',   emit: sm

    shell:
    '''
    set -euo pipefail
    echo "[3/6] Preparing search space (building Sage FASTA) ..."
    mkdir -p build
    PROT=""; [ -s proteome.fasta ] && PROT="--proteome proteome.fasta"
    if [ "!{params.peptides_format}" = "module1" ]; then
        02_build_search_fasta_new.py --peptides "!{peptides}" $PROT \
            --out build/search.fasta --min-len !{params.min_pep_len}
    else
        cp "!{peptides}" build/peptides.tsv
        # legacy Chong layout (peptide/group/source): use a provided --build_config, else synthesise one
        if [ -s config.yaml ]; then
            CFG=config.yaml
        elif [ -s proteome.fasta ]; then
            printf 'build:\n  min_pep_len: %s\n  out_fasta: build/search.fasta\ninputs:\n  canonical_proteome_fasta: proteome.fasta\n' "!{params.min_pep_len}" > config.gen.yaml
            CFG=config.gen.yaml
        else
            printf 'build:\n  min_pep_len: %s\n  out_fasta: build/search.fasta\ninputs: {}\n' "!{params.min_pep_len}" > config.gen.yaml
            CFG=config.gen.yaml
        fi
        02_build_search_fasta.py "$CFG"
    fi
    N=$(grep -c '^>' build/search.fasta || true)
    echo "Search FASTA ready: ${N} sequences (CANON + NONCANON)" | tee buildfasta.summary
    '''
    stub:
    """
    mkdir -p build
    printf '>chong_pub_0|grp=NONCANON|src=pub_test\\nKYKDRTNILF\\n>sp|P0|grp=CANON|src=uniprot\\nPEPTIDEK\\n' > build/search.fasta
    echo "Search FASTA ready: 2 sequences (CANON + NONCANON)" > buildfasta.summary
    """

}

process SAGE {
    tag { "HLA-${cls}" }
    publishDir path: { "${params.outdir}/04_sage/HLA-${cls}" }, mode: 'copy', pattern: 'sage_*/*.tsv'
    publishDir "${params.outdir}/04_sage", mode: 'copy', pattern: '*.summary'

    input:
    tuple val(cls), path(sage_cfg), val(batch), val(mzml_prefix)
    path fasta

    output:
    tuple val(cls), path("sage_${cls}/results.sage.tsv"), val(mzml_prefix), emit: res
    path "sage_${cls}.summary", emit: sm

    shell:
    '''
    set -euo pipefail
    echo "[4/6] Sage search (HLA-!{cls}) ..."

    # ---- INFERRED STEP: pull mzML from S3 to local (your convert step deletes local mzML) ----
    mkdir -p mzml_!{cls} sage_!{cls}
    aws s3 sync "!{mzml_prefix}" mzml_!{cls}/ --exclude "*" --include "*.mzML"
    NMZ=$(find mzml_!{cls} -name '*.mzML' | wc -l)
    [ "$NMZ" -eq 0 ] && { echo "ERROR: no mzML under !{mzml_prefix}"; exit 1; }

    # ---- point the Sage config at OUR fasta + output dir (leaves the rest of your config intact) ----
    python3 - <<'PY'
import json
c = json.load(open("!{sage_cfg}"))
c.setdefault("database", {})["fasta"] = "!{fasta}"
c["output_directory"] = "sage_!{cls}"
json.dump(c, open("sage_run.json", "w"), indent=2)
PY

    # ---- run Sage exactly as run_sage.sh does (no-enzyme search lives in your config) ----
    sage sage_run.json mzml_!{cls}/*.mzML --batch-size !{batch}

    NP=$(( $(wc -l < sage_!{cls}/results.sage.tsv) - 1 ))
    echo "Sage HLA-!{cls}: ${NP} PSMs from ${NMZ} mzML files" | tee sage_!{cls}.summary
    '''
    stub:
    """
    mkdir -p sage_${cls}
    printf 'peptide\\tproteins\\trank\\tsage_discriminant_score\\tlabel\\n'  > sage_${cls}/results.sage.tsv
    printf 'KYKDRTNILF\\tchong_pub_0|grp=NONCANON\\t1\\t2.5\\t1\\n'         >> sage_${cls}/results.sage.tsv
    printf 'PEPTIDEK\\tsp|P0|grp=CANON\\t1\\t1.8\\t1\\n'                    >> sage_${cls}/results.sage.tsv
    echo "Sage HLA-${cls}: 2 PSMs from 3 mzML files" > sage_${cls}.summary
    """

}

process MS2RESCORE {
    tag { "HLA-${cls}" }
    publishDir path: { "${params.outdir}/05_rescore/HLA-${cls}" }, mode: 'copy', pattern: 'rescore_*/**'
    publishDir "${params.outdir}/05_rescore", mode: 'copy', pattern: '*.summary'

    input:
    tuple val(cls), path(sage_tsv), val(mzml_prefix)
    path ms2_cfg

    output:
    tuple val(cls), path("rescore_${cls}/rescored_psms.tsv"), emit: res
    path "rescore_${cls}.summary", emit: sm

    shell:
    '''
    set -euo pipefail
    echo "[5/6] MS2Rescore rescoring (HLA-!{cls}) ..."

    mkdir -p rescore_!{cls} mzml_!{cls}
    aws s3 sync "!{mzml_prefix}" mzml_!{cls}/ --exclude "*" --include "*.mzML"

    # ---- rank-1 filter (generalised from run_ms2rescore.sh: locate 'rank' column dynamically) ----
    RANKCOL=$(head -1 "!{sage_tsv}" | tr '\\t' '\\n' | grep -n '^rank$' | cut -d: -f1 || true)
    if [ -z "${RANKCOL:-}" ]; then
        echo "   no 'rank' column found -> rescoring all rows"
        cp "!{sage_tsv}" rank1.sage.tsv
    else
        awk -F'\\t' -v r="$RANKCOL" 'NR==1 || $r==1' "!{sage_tsv}" > rank1.sage.tsv
    fi

    # ---- run MS2Rescore; DeepLC RT calibration needs confident target PSMs, so fall back to
    #      MS2PIP-only if it fails (e.g. a run with too few IDs -> "no target PSMs for calibration") ----
    strip_deeplc () {   # $1 in-config -> $2 out-config (removes the deeplc feature generator)
        python3 - "$1" "$2" <<'PYEOF'
import json, sys
c = json.load(open(sys.argv[1]))
((c.get("ms2rescore") or {}).get("feature_generators") or {}).pop("deeplc", None)
json.dump(c, open(sys.argv[2], "w"), indent=2)
PYEOF
    }
    CFG="!{ms2_cfg}"
    if [ "!{params.rescore_use_deeplc}" = "false" ]; then
        strip_deeplc "!{ms2_cfg}" cfg_use.json && CFG=cfg_use.json && echo "   DeepLC disabled (--rescore_use_deeplc false)"
    fi
    if ! ms2rescore -p rank1.sage.tsv -s mzml_!{cls}/ -o rescore_!{cls}/ -c "$CFG" -n !{params.rescore_threads}; then
        echo "[WARN] MS2Rescore failed with the primary config — commonly DeepLC retention-time calibration" >&2
        echo "       when a run has too few confident target PSMs (make sure your search DB includes a" >&2
        echo "       canonical proteome: --canonical_proteome). Retrying WITHOUT DeepLC ..." >&2
        rm -rf rescore_!{cls}; mkdir -p rescore_!{cls}
        strip_deeplc "!{ms2_cfg}" cfg_nodeeplc.json
        ms2rescore -p rank1.sage.tsv -s mzml_!{cls}/ -o rescore_!{cls}/ -c cfg_nodeeplc.json -n !{params.rescore_threads}
    fi

    # ---- resolve the rescored PSM table (name varies by MS2Rescore/mokapot version+config) ----
    PSM=""
    for pat in 'rescore_!{cls}'/*.psms.tsv 'rescore_!{cls}'/*mokapot*psms* 'rescore_!{cls}'/*.psms.txt ; do
        for f in $pat ; do [ -e "$f" ] && { PSM="$f"; break; } ; done
        [ -n "$PSM" ] && break
    done
    if [ -z "$PSM" ]; then
        echo "ERROR: no rescored PSM table in rescore_!{cls}/ — contents:"; ls -la rescore_!{cls}/; exit 1
    fi
    cp "$PSM" rescore_!{cls}/rescored_psms.tsv

    NP=$(( $(wc -l < rescore_!{cls}/rescored_psms.tsv) - 1 ))
    echo "MS2Rescore HLA-!{cls}: ${NP} rescored PSMs (from ${PSM##*/})" | tee rescore_!{cls}.summary
    '''
    stub:
    """
    mkdir -p rescore_${cls}
    printf 'peptide\\tproteins\\tmokapot score\\tmokapot q-value\\tis_decoy\\n'   > rescore_${cls}/rescored_psms.tsv
    printf 'KYKDRTNILF\\tchong_pub_0|grp=NONCANON\\t3.10\\t0.001\\tFalse\\n'     >> rescore_${cls}/rescored_psms.tsv
    printf 'PEPTIDEK\\tsp|P0|grp=CANON\\t2.40\\t0.004\\tFalse\\n'                >> rescore_${cls}/rescored_psms.tsv
    echo "MS2Rescore HLA-${cls}: 2 rescored PSMs (from rescored_psms.tsv)" > rescore_${cls}.summary
    """

}

process GROUP_FDR {
    tag { "HLA-${cls}" }
    publishDir path: { "${params.outdir}/06_group_fdr/HLA-${cls}" }, mode: 'copy', pattern: 'results/*.tsv'
    publishDir "${params.outdir}/06_group_fdr", mode: 'copy', pattern: '*.summary'

    input:
    tuple val(cls), path(rescored)
    path fdr_cfg, stageAs: 'config.yaml'

    output:
    tuple val(cls), path('results/*.tsv'), emit: tables
    path "fdr_${cls}.summary", emit: sm

    shell:
    '''
    set -euo pipefail
    echo "[6/6] Group-specific FDR + recovery (HLA-!{cls}) ..."
    # NOTE: 04_group_fdr_recovery.py reads validation_targets.tsv from ITS OWN dir (bin/).
    mkdir -p results
    04_group_fdr_recovery.py "!{rescored}" config.yaml

    OUT=$(ls -1 results/*.tsv | head -1)
    NT=$(ls results/*.tsv | wc -l)
    NROWS=$(( $(wc -l < "$OUT") - 1 ))
    echo "Group FDR HLA-!{cls}: ${NT} output tables, top table ${NROWS} rows" | tee fdr_!{cls}.summary
    '''
    stub:
    """
    mkdir -p results
    printf 'peptide_norm\\tgroup\\tqvalue\\n'    > results/group_fdr_recovery.tsv
    printf 'KYKDRTNLLF\\tNONCANON\\t0.001\\n'   >> results/group_fdr_recovery.tsv
    echo "Group FDR HLA-${cls}: 1 output tables, top table 1 rows" > fdr_${cls}.summary
    """

}

// ----------------------------------------------------------------------------
// Per-class spec builder  (fan out at SAGE onward for 'both')
// ----------------------------------------------------------------------------
def classSpec(String c) {
    def cfg   = (c == 'I') ? file(params.sage_config_I)  : file(params.sage_config_II)
    def batch = (c == 'I') ? params.batch_size_I         : params.batch_size_II
    def mz    = (c == 'I') ? (params.mzml_I  ?: params.s3_mzml)
                           : (params.mzml_II ?: params.s3_mzml)
    return tuple(c, cfg, batch, mz)
}

// ----------------------------------------------------------------------------
// Workflow
// ----------------------------------------------------------------------------
// ---- spectrum + chromatogram for one random canonical & one non-canonical identified peptide ----
process SPECTRUM_PLOTS {
    tag { "HLA-${cls}" }
    publishDir path: { "${params.outdir}/07_id_diagnostics/HLA-${cls}/spectra" }, mode: 'copy'
    input:
    tuple val(cls), path(sage_tsv), val(mzml_prefix)
    output:
    path 'spectrum_*', emit: fig, optional: true
    shell:
    '''
    set -euo pipefail
    echo "[QC] Spectrum + chromatogram plots (HLA-!{cls}) ..."
    mkdir -p mzml_!{cls}
    aws s3 sync "!{mzml_prefix}" mzml_!{cls}/ --exclude "*" --include "*.mzML"
    python3 "$(command -v !{params.spectrum_script})" \
        --sage "!{sage_tsv}" --mzml-dir mzml_!{cls} --hla-class "!{cls}" \
        --seed !{params.spectrum_seed} --outdir .
    '''
    stub:
    """
    : > spectrum_HLA-${cls}_noncanonical_DEMO_scan1_spectrum.png
    : > spectrum_HLA-${cls}_canonical_DEMO_scan2_spectrum.png
    """
}

// ---- identification diagnostics: Sage (pre) vs MS2Rescore (post), per HLA class ----
process ID_DIAGNOSTICS {
    tag { "HLA-${cls}" }
    publishDir path: { "${params.outdir}/07_id_diagnostics/HLA-${cls}" }, mode: 'copy'
    input:
    tuple val(cls), path(sage_tsv), path(rescored)
    output:
    path 'id_diagnostics_*.png', emit: fig, optional: true
    path 'id_diagnostics_*.tsv', emit: tab
    shell:
    '''
    set -euo pipefail
    echo "[QC] Identification diagnostics (HLA-!{cls}) ..."
    python3 "$(command -v !{params.id_diag_script})" \
        --sage "!{sage_tsv}" --rescored "!{rescored}" \
        --hla-class "!{cls}" --fdr !{params.id_fdr} --outdir .
    '''
    stub:
    """
    : > id_diagnostics_HLA-${cls}.png
    printf 'stage\\thla_class\\tpsms\\tpeptides\\tcanonical\\tnoncanonical\\n' > id_diagnostics_HLA-${cls}.tsv
    printf 'post\\t${cls}\\t3850\\t3850\\t2731\\t1119\\n' >> id_diagnostics_HLA-${cls}.tsv
    """
}

workflow {

    // ---- validate the essentials up front ----
    if (params.from_mzml) {
        if (!params.s3_mzml) { error "--from_mzml needs --s3_mzml (an existing prefix of .mzML files)" }
    } else {
        def srcs = [params.pride, params.repo, params.local].findAll { it }
        if (srcs.size() == 0) { error "Provide a data source: --pride <acc>, --repo <acc|URL>, --local <dir|s3://>, or --from_mzml" }
        if (srcs.size() > 1)  { error "Provide only ONE of --pride / --repo / --local" }
        if (!params.s3_mzml)  { error "Missing --s3_mzml (s3:// prefix for .mzML files)" }
        def localIsS3 = params.local && params.local.toString().startsWith('s3://')
        if (!params.s3_raw && !localIsS3) { error "Missing --s3_raw (s3:// prefix to stage raw files)" }
    }
    if (!(params.hla_class in ['I', 'II', 'both'])) {
        error "--hla_class must be one of: I | II | both (got '${params.hla_class}')"
    }
    if (!(params.vendor in ['thermo', 'bruker', 'sciex'])) {
        error "--vendor must be 'thermo' (.raw), 'bruker' (.d) or 'sciex' (.wiff) (got '${params.vendor}')"
    }
    if (!params.peptides_tsv) {
        error "Missing --peptides_tsv (Module 1 translation TSV; for the legacy layout use --peptides_format chong)"
    }
    if (!(params.peptides_format in ['module1', 'chong'])) {
        error "--peptides_format must be 'module1' or 'chong' (got '${params.peptides_format}')"
    }

    println """
    ===========================================================================
     Module_2  |  MS/MS validation
     mode      : ${params.from_mzml ? 'from existing mzML (skip download+convert)' : (params.local ? 'local/S3 ingest -> mzML' : 'repository -> raw -> mzML')}
     vendor    : ${params.vendor}${params.ext ? ' (ext ' + params.ext + ')' : ''}${params.pride_files ? '   files: ' + params.pride_files : ''}
     source    : ${params.pride ?: params.repo ?: params.local ?: '(from mzML)'}
     S3 raw    : ${params.s3_raw ?: '(n/a)'}
     S3 mzML   : ${params.s3_mzml}
     HLA class : ${params.hla_class}
    ===========================================================================
    """.stripIndent()

    def classes = (params.hla_class == 'both') ? ['I', 'II'] : [params.hla_class]
    class_ch = Channel.fromList(classes.collect { classSpec(it) })

    // ---- resolve the mzML prefix: either build it (download+convert) or reuse an existing one ----
    if (params.from_mzml) {
        mzml = Channel.value(params.s3_mzml)
    } else {
        def localIsS3 = params.local && params.local.toString().startsWith('s3://')
        def raw_ch
        if (params.local && localIsS3) {
            // in-house data already staged in S3: convert straight from it (no download, no copy)
            raw_ch = Channel.value(params.local)
            stageBar(1, 'Using in-house S3 raw prefix (no download)', "source: ${params.local}")
        } else if (params.local) {
            li = LOCAL_INGEST(params.local, params.s3_raw)
            raw_ch = li.s3_raw
            li.sm.map { it.text.trim() }.subscribe { d -> stageBar(1, 'Ingesting in-house data', d) }
        } else if (params.repo) {
            dl = DOWNLOAD_REPO(params.repo, params.s3_raw, params.pride_files ?: '')
            raw_ch = dl.s3_raw
            dl.sm.map { it.text.trim() }.subscribe { d -> stageBar(1, "Downloading from repository (${params.repo})", d) }
        } else {
            dl = DOWNLOAD_PRIDE(params.pride, params.s3_raw, params.pride_files ?: '')
            raw_ch = dl.s3_raw
            dl.sm.map { it.text.trim() }.subscribe { d -> stageBar(1, 'Downloading files from PRIDE', d) }
        }
        cv = CONVERT_STREAM(raw_ch, params.s3_mzml)
        mzml = cv.s3_mzml
        cv.sm.map { it.text.trim() }.subscribe { d -> stageBar(2, "Converting ${params.vendor} raw to mzML", d) }
    }

    pep  = file(params.peptides_tsv, checkIfExists: true)
    prot = params.canonical_proteome ? file(params.canonical_proteome, checkIfExists: true)
                                     : file("${projectDir}/assets/NO_PROTEOME")
    cfg  = params.build_config       ? file(params.build_config, checkIfExists: true)
                                     : file("${projectDir}/assets/NO_CONFIG")
    bf = BUILD_FASTA(mzml, cfg, pep, prot)

    sg = SAGE(class_ch, bf.fasta)
    rs = MS2RESCORE(sg.res, file(params.ms2rescore_config))
    fd = GROUP_FDR(rs.res, file(params.fdr_config))

    // identification diagnostics: join Sage (pre) with MS2Rescore (post) by HLA class
    if (params.run_id_plots) {
        ID_DIAGNOSTICS(sg.res.map { c, sage, mz -> tuple(c, sage) }.join(rs.res))
    }
    if (params.run_spectrum_plots) {
        SPECTRUM_PLOTS(sg.res)
    }

    // ---- tqdm-style stage bars (fire on completion, so they print in order) ----
    bf.sm.map { it.text.trim() }.subscribe { d -> stageBar(3, 'Preparing search space',        d) }
    sg.sm.map { it.text.trim() }.subscribe { d -> stageBar(4, 'Sage database search',          d) }
    rs.sm.map { it.text.trim() }.subscribe { d -> stageBar(5, 'MS2Rescore rescoring',          d) }
    fd.sm.map { it.text.trim() }.subscribe { d -> stageBar(6, 'Group FDR + recovery',          d) }
    workflow.onComplete = {
        if (workflow.success) println "  ${renderBar(6)}"
        println ''
        println (workflow.success
            ? "  Module_2 complete  ->  ${params.outdir}/  (status: OK)"
            : "  Module_2 FAILED  (see .nextflow.log)")
    }
}
