#!/usr/bin/env nextflow
// ===============================================================================
//  Module_1.nf  --  Nucleotide -> peptide translation (Module 1)
// -------------------------------------------------------------------------------
//  Produces the origin-annotated protein FASTA (+ TSV) that seeds Module 2's search DB.
//  Rules-based, extensible translation: origin (from the header) -> #frames / strand /
//  genetic code / filters, via two editable tables in assets/.
//
//  INPUT MODES (resolved from which args are given; each fans out over 1..N files):
//    --fasta <glob>                          already a nucleotide FASTA -> translate by header rules
//    --fastq <glob> --metadata <tsv>         FASTQ + metadata -> FASTA (origin stamped) -> translate
//    --fastq <glob> --reference <db>         FASTQ only -> align to a combined origin-annotated DB,
//                                            derive origin from the mapped reference -> translate  [SCAFFOLD]
//    --vcf <glob>   --reference <genome>     VCF -> reference-window + ALT -> translate (variant rule)
//    --vcf <glob> --fastq <glob> --reference --readaware
//                                            read-aware variant arm (align + IsoVar-style)          [SCAFFOLD]
//
//  Why the variant path needs --reference: a variant is a coordinate + allele; reconstructing the
//  sequence around it needs reference context (or reads already aligned to it). A raw FASTQ has no
//  coordinates, so coordinate extraction requires a reference (or the read-aware arm, which aligns).
//
//  Translation rules (assets/translation_rules.tsv), per your spec:
//    variant_snv, frameshift_indel, fusion, aberrant_splicing, viral -> 3 frames
//    microbial -> 3 frames (genetic code 11);  erv, ncrna, ncorf, intronic, intergenic,
//    out_of_frame, unknown -> 6 frames.  Add/modify rows to add or change rules.
//  Origin is inferred from the header (assets/origin_patterns.tsv); an explicit origin=<class>
//  token always wins. Multiple origins in one FASTA are handled per-record.
//
//  OUTPUTS (per input file):
//    peptides_<base>.fasta   protein sequences; header = key=value provenance
//    peptides_<base>.tsv     one row per ORF: sequence + provenance columns (feeds Module 2)
//    visualization/          IGV-style panels for 2 non-canonical/SNV peptides (PNG + PDF) +
//                            viz_selected.tsv   (--visualize false to skip; --viz_backend gviz for R/Gviz)
// ===============================================================================

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------- params
params.fasta      = null
params.fastq      = null
params.vcf        = null
params.metadata   = null
params.reference  = null
params.sample     = null
params.tissue     = null
params.min_aa_len = null
params.flank      = 30
params.readaware  = false          // vcf+fastq: use the read-aware (IsoVar-style) scaffold
params.patterns   = null           // override assets/origin_patterns.tsv
params.rules      = null           // override assets/translation_rules.tsv
params.outdir     = 'results_module1'

params.visualize  = true            // IGV-style view of 2 non-canonical/SNV peptides (set false to skip)
params.viz_backend = 'python'       // 'python' (matplotlib, default) or 'gviz' (R/Bioconductor)
params.viz_n      = 2
params.viz_seed   = 12              // -1 for nondeterministic selection

params.fastq2fasta_script = 'fastq_to_fasta.py'
params.variant_script     = 'extract_variant_seqs.py'
params.translate_script   = 'translate_orfs.py'
params.viz_script         = 'plot_peptide_igv.py'
params.gviz_script        = 'plot_peptide_gviz.R'
params.aligner            = 'STAR'  // or bwa-mem / minimap2 (scaffold)

// ============================================================================ processes

process FASTQ2FASTA_META {
    tag "${fastq.simpleName}"
    publishDir "${params.outdir}/nt_fasta", mode: 'copy'
    input:
    path fastq
    path metadata
    output:
    path 'ntseq_*.fasta', emit: fasta
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{fastq}" | sed 's/\\.[^.]*$//')
    S="!{params.sample ?: ''}"; SARG=""; [ -n "$S" ] && SARG="--sample $S"
    T="!{params.tissue ?: ''}"; TARG=""; [ -n "$T" ] && TARG="--tissue $T"
    python3 "$(command -v !{params.fastq2fasta_script})" \
        --fastq "!{fastq}" --metadata "!{metadata}" $SARG $TARG \
        --out "ntseq_${base}.fasta"
    '''
    stub:
    """
    base=\$(basename "${fastq}" | sed 's/\\.[^.]*\$//')
    printf '>r1 origin=fusion sample=NA tissue=unknown readid=r1\\nATGGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA\\n' > ntseq_\${base}.fasta
    """
}

process VARIANT_EXTRACT {
    tag "${vcf.simpleName}"
    publishDir "${params.outdir}/nt_fasta", mode: 'copy'
    input:
    path vcf
    path reference
    output:
    path 'ntseq_*.fasta', emit: fasta
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{vcf}" | sed 's/\\.[^.]*$//')
    S="!{params.sample ?: ''}"; SARG=""; [ -n "$S" ] && SARG="--sample $S"
    T="!{params.tissue ?: ''}"; TARG=""; [ -n "$T" ] && TARG="--tissue $T"
    python3 "$(command -v !{params.variant_script})" \
        --vcf "!{vcf}" --reference "!{reference}" --flank !{params.flank} \
        $SARG $TARG --out "ntseq_${base}.fasta"
    '''
    stub:
    """
    base=\$(basename "${vcf}" | sed 's/\\.[^.]*\$//')
    SEQ="ATG\$(python3 -c "print('GCT'*18)")TAA"
    printf '>chrT:40:T:A origin=variant_snv variant=chrT:40:T:A gene=TESTG coords=chrT:10-70 strand=+ type=SNV sample=NA tissue=unknown\\n%s\\n' "\$SEQ" > ntseq_\${base}.fasta
    """
}

// ---- SCAFFOLD: FASTQ-only -> align to a combined origin-annotated DB -> origin-tagged FASTA ----
process ALIGN_EXTRACT {
    tag "${fastq.simpleName}"
    publishDir "${params.outdir}/nt_fasta", mode: 'copy'
    input:
    path fastq
    path reference
    output:
    path 'ntseq_*.fasta', emit: fasta
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{fastq}" | sed 's/\\.[^.]*$//')
    # >>> SCAFFOLD: needs an aligner + samtools + a combined DB whose headers carry origin. <<<
    # The reference DB must contain human + ncORF + fusion + ERV + viral/microbial contigs, each
    # header origin-annotated (e.g. `... origin=erv` or a name the classifier recognizes). Then:
    #   1) align reads          e.g. STAR --genomeDir <idx> --readFilesIn !{fastq} ...
    #                            or  bwa-mem2 mem <idx> !{fastq} | samtools sort -o aln.bam
    #   2) keep primary mapped  samtools view -bF 0x904 aln.bam > mapped.bam
    #   3) for each mapped read, emit a FASTA record: read sequence, header = the mapped
    #      reference's origin-annotated header (so TRANSLATE applies the right rule).
    # Replace the placeholder below with that alignment + header-join step.
    echo "[align-extract] ALIGNER=!{params.aligner} ref=!{reference} reads=!{fastq}" >&2
    : > "ntseq_${base}.fasta"
    if [ ! -s "ntseq_${base}.fasta" ]; then
        echo "ERROR: ALIGN_EXTRACT is a scaffold - wire your aligner + DB header-join here." >&2
        exit 3
    fi
    '''
    stub:
    """
    base=\$(basename "${fastq}" | sed 's/\\.[^.]*\$//')
    printf '>read1 origin=erv sample=NA tissue=unknown readid=read1\\nATG%sTAA\\n' "\$(python3 -c "print('GCT'*12)")" > ntseq_\${base}.fasta
    """
}

// ---- SCAFFOLD: read-aware variant arm (VCF + FASTQ + reference), IsoVar-style ----
process VARIANT_READAWARE {
    tag "${vcf.simpleName}"
    publishDir "${params.outdir}/nt_fasta", mode: 'copy'
    input:
    path vcf
    path fastq
    path reference
    output:
    path 'ntseq_*.fasta', emit: fasta
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{vcf}" | sed 's/\\.[^.]*$//')
    # >>> SCAFFOLD: align !{fastq} to !{reference}, then translate reads spanning each variant in
    #     !{vcf} (e.g. openvax/isovar ProteinSequenceCreator), emitting variant-spanning contigs
    #     with origin=variant_snv / frameshift_indel headers. Wire isovar (needs aligned BAM). <<<
    echo "[variant-readaware] vcf=!{vcf} reads=!{fastq} ref=!{reference}" >&2
    echo "ERROR: VARIANT_READAWARE is a scaffold - wire alignment + IsoVar here (or use the" >&2
    echo "       reference-based path: --vcf --reference without --readaware)." >&2
    exit 3
    '''
    stub:
    """
    base=\$(basename "${vcf}" | sed 's/\\.[^.]*\$//')
    printf '>chrT:40:T:A origin=variant_snv variant=chrT:40:T:A gene=TESTG coords=chrT:10-70 strand=+ sample=NA tissue=unknown\\nATG%sTAA\\n' "\$(python3 -c "print('GCT'*18)")" > ntseq_\${base}.fasta
    """
}

process TRANSLATE {
    tag "${ntfasta.simpleName}"
    publishDir "${params.outdir}/peptides", mode: 'copy'
    input:
    path ntfasta
    path patterns
    path rules
    output:
    path 'peptides_*.fasta', emit: protein
    path 'peptides_*.tsv',   emit: tsv
    shell:
    '''
    set -euo pipefail
    base=$(basename "!{ntfasta}" | sed 's/^ntseq_//; s/\\.[^.]*$//')
    S="!{params.sample ?: ''}"; SARG=""; [ -n "$S" ] && SARG="--sample $S"
    T="!{params.tissue ?: ''}"; TARG=""; [ -n "$T" ] && TARG="--tissue $T"
    M="!{params.min_aa_len ?: ''}"; MIN=""; [ -n "$M" ] && MIN="--min-aa-len $M"
    python3 "$(command -v !{params.translate_script})" \
        --fasta "!{ntfasta}" --patterns "!{patterns}" --rules "!{rules}" \
        $SARG $TARG $MIN --out "peptides_${base}"
    NP=$(( $(wc -l < "peptides_${base}.tsv") - 1 ))
    echo "[translate] ${base}: ${NP} ORF peptides" >&2
    '''
    stub:
    """
    base=\$(basename "${ntfasta}" | sed 's/^ntseq_//; s/\\.[^.]*\$//')
    printf 'id\\tsequence\\tclass\\torigin\\tsource\\tsample\\ttissue\\tlocus\\tstrand\\tframe\\tgenetic_code\\torf_offset_nt\\taa_len\\tvariant\\tgene\\treadid\\n' > peptides_\${base}.tsv
    printf 'r1|orf1\\tMAAAAAAAAAA\\tnoncanonical\\tfusion\\tNA\\tNA\\tunknown\\tNA\\t+\\t0\\t1\\t0\\t11\\tNA\\tNA\\tr1\\n' >> peptides_\${base}.tsv
    printf '>r1|orf1 class=noncanonical origin=fusion\\nMAAAAAAAAAA\\n' > peptides_\${base}.fasta
    """
}

process VISUALIZE {
    tag "igv"
    publishDir "${params.outdir}/visualization", mode: 'copy'
    input:
    path tsvs
    output:
    path 'igv_*.png',        emit: png, optional: true
    path 'igv_peptides.pdf', emit: pdf, optional: true
    path 'viz_selected.tsv', emit: sel
    shell:
    '''
    set -euo pipefail
    if [ "!{params.viz_backend}" = "gviz" ]; then
        Rscript "$(command -v !{params.gviz_script})" --tsv *.tsv --outdir . --n !{params.viz_n} --seed !{params.viz_seed}
    else
        python3 "$(command -v !{params.viz_script})" --tsv *.tsv --outdir . --n !{params.viz_n} --seed !{params.viz_seed}
    fi
    '''
    stub:
    """
    printf 'id\\torigin\\tclass\\n' > viz_selected.tsv
    printf 'stub|orf1\\tvariant_snv\\tvariant\\n' >> viz_selected.tsv
    : > igv_stub.png
    : > igv_peptides.pdf
    """
}

// ============================================================================ main
workflow {

    // ---- resolve input mode ----
    def haveFasta = params.fasta as boolean
    def haveFastq = params.fastq as boolean
    def haveVcf   = params.vcf as boolean
    def haveMeta  = params.metadata as boolean
    def haveRef   = params.reference as boolean

    def mode
    if (haveFasta)                          { mode = 'FASTA' }
    else if (haveVcf)                       { mode = (haveFastq && params.readaware) ? 'VARIANT_READAWARE' : 'VARIANT' }
    else if (haveFastq && haveMeta)         { mode = 'FASTQ_META' }
    else if (haveFastq && haveRef)          { mode = 'FASTQ_ALIGN' }
    else if (haveFastq)                     { error "With --fastq you must give either --metadata (Option 2) or --reference (Option 3)" }
    else                                    { error "No input. Provide --fasta, or --fastq (+--metadata|--reference), or --vcf (+--reference)" }

    if (mode in ['VARIANT','VARIANT_READAWARE','FASTQ_ALIGN'] && !haveRef) {
        error "Mode ${mode} needs --reference (coordinate/alignment context)."
    }

    def patt  = file(params.patterns ?: "${projectDir}/assets/origin_patterns.tsv", checkIfExists: true)
    def ruleF = file(params.rules    ?: "${projectDir}/assets/translation_rules.tsv", checkIfExists: true)

    println """
    ===========================================================================
     Module_1  |  Nucleotide -> peptide translation
     mode      : ${mode}
     sample    : ${params.sample ?: '(from headers)'}   tissue: ${params.tissue ?: '(from headers)'}
    ===========================================================================
    """.stripIndent()

    // ---- produce the nucleotide FASTA channel per mode (fan out over 1..N inputs) ----
    if (mode == 'FASTA') {
        nt = Channel.fromPath(params.fasta, checkIfExists: true)
    }
    else if (mode == 'FASTQ_META') {
        fq  = Channel.fromPath(params.fastq, checkIfExists: true)
        md  = file(params.metadata, checkIfExists: true)
        nt  = FASTQ2FASTA_META(fq, md).fasta
    }
    else if (mode == 'FASTQ_ALIGN') {
        fq  = Channel.fromPath(params.fastq, checkIfExists: true)
        ref = file(params.reference, checkIfExists: true)
        nt  = ALIGN_EXTRACT(fq, ref).fasta
    }
    else if (mode == 'VARIANT') {
        vc  = Channel.fromPath(params.vcf, checkIfExists: true)
        ref = file(params.reference, checkIfExists: true)
        nt  = VARIANT_EXTRACT(vc, ref).fasta
    }
    else { // VARIANT_READAWARE
        vc  = Channel.fromPath(params.vcf, checkIfExists: true)
        fq  = file(params.fastq, checkIfExists: true)
        ref = file(params.reference, checkIfExists: true)
        nt  = VARIANT_READAWARE(vc, fq, ref).fasta
    }

    // ---- translate (rules applied per record) ----
    out = TRANSLATE(nt, patt, ruleF)

    // ---- IGV-style visualization of 2 non-canonical/SNV peptides ----
    if (params.visualize) {
        VISUALIZE(out.tsv.collect())
    }

    out.tsv.collect().subscribe { l ->
        println "\n  Translation complete: ${l.size()} peptide table(s) -> ${params.outdir}/peptides/"
    }

    workflow.onComplete = {
        println ''
        println (workflow.success
            ? "  Module_1 complete  ->  ${params.outdir}/  (status: OK)"
            : "  Module_1 FAILED  (see .nextflow.log)")
    }
}
