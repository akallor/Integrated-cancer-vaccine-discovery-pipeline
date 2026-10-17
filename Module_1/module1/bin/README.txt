Module 1 - bin/ scripts (nucleotide -> peptide translation)
===========================================================

m1_common.py            Shared, imported by translate_orfs.py (and reusable): loads the two rule
                        tables, classifies a header -> origin_class (explicit `origin=` wins),
                        and parses provenance fields (coords, strand, gene, sample, tissue,
                        variant, source) from heterogeneous headers.

translate_orfs.py       CORE translator. FASTA in -> protein FASTA + TSV out.
  --fasta F             input nucleotide FASTA
  --out PREFIX          writes PREFIX.fasta and PREFIX.tsv
  --patterns / --rules  override assets/origin_patterns.tsv / assets/translation_rules.tsv
  --sample / --tissue   fallback values when a header lacks them
  --min-aa-len N        global override of every rule's min ORF length
  Behaviour: per record, classify origin -> look up rule (#frames 3|6, strand, NCBI genetic
  code, min length, require-start) -> translate the chosen frames -> split each frame on stop
  codons -> emit one protein record per surviving ORF. Frames 0-2 = forward, 3-5 = reverse.

fastq_to_fasta.py       INPUT Option 2 (FASTQ + metadata -> FASTA).
  --fastq / --metadata / --out [--sample --tissue --origin]
  Metadata TSV/CSV: col1 = read id (or key), col2 = origin (SNV/ERV/ncRNA/Fusion/...); optional
  sample, tissue columns. A single metadata row acts as a file-wide default. Origin labels are
  normalized to translator classes. Preserves the read id, strips quality, stamps origin=.

extract_variant_seqs.py INPUT variant path (VCF + reference genome -> mutant nucleotide FASTA).
  --vcf / --reference / --out [--flank 30 --sample --tissue]
  Per VCF record: pull a reference window around the variant, apply the first ALT allele, emit a
  nucleotide record headed with variant/gene/coords. SNV/MNV -> variant_snv; length-changing ->
  frameshift_indel. (Split multiallelic sites upstream.) Reference is indexed lazily.

plot_peptide_igv.py     IGV-style visualization (default backend, matplotlib).
  --tsv F [F ...] --outdir D [--n 2 --seed 12]
  Picks n peptides (prefers one SNV/variant + one non-canonical) with a chr locus and draws each
  on genomic coordinates: coordinate axis, source locus, ORF block (strand arrow + residue ruler),
  and a variant lollipop for SNV/indel. Writes igv_<id>.png per peptide, igv_peptides.pdf, and
  viz_selected.tsv. Peptide genomic span = locus + orf_offset_nt + strand + aa_len (indel windows
  approximate). Writes only a note (exit 0) if nothing has coordinates.

plot_peptide_gviz.R     Same view via Bioconductor Gviz (alternate backend; --viz_backend gviz).
  Needs R + Gviz; NOT run in this build env (no CRAN/Bioconductor there). Install with
  BiocManager::install("Gviz").

Dependencies: Python 3, biopython, pandas, matplotlib (viz).  R + Gviz only for the gviz backend.  All scripts are invoked from Module_1.nf via PATH
(`command -v <script>`); keep them executable (chmod +x).

Scaffolded arms (in Module_1.nf, not here): ALIGN_EXTRACT (FASTQ-only: STAR/BWA-MEM + samtools
against a combined origin-annotated DB) and VARIANT_READAWARE (VCF+FASTQ read-aware, IsoVar-style).
Both carry the intended tool invocations and exit 3 until wired.
