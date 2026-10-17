Module 1 - assets/ rule tables (both are plain TSV, edit freely)
================================================================

origin_patterns.tsv   priority <TAB> origin_class <TAB> regex
  Lower priority = tested first; the first regex that matches the FASTA header wins. Regexes are
  matched case-insensitively against the whole header. Add a row to recognize a new header style,
  or raise a row's priority (lower number) to make it win over a more general pattern.
  An explicit `origin=<class>` token in a header ALWAYS overrides these patterns.

translation_rules.tsv origin_class <TAB> frames(1|3|6) <TAB> strand(forward|both|auto)
                      <TAB> genetic_code <TAB> min_aa_len <TAB> require_start <TAB> cds
  frames=1 = one forward frame; frames=3 = one strand (forward, or header strand when auto);
  frames=6 = both strands. genetic_code is an NCBI table (1=standard, 11=bacterial).
  cds=true: treat the record as an ALREADY-DEFINED ORF -> single frame, initiator codon rendered
  as M (non-AUG starts TTG/CTG/GTG... included, per Ribo-seq/GENCODE convention), terminal stop
  dropped. This is the lean choice for the MS search DB when the sequence IS the ORF (e.g. ncorf /
  Ribo-seq ncORFs). Leave cds=false for discovery classes that need multi-frame translation.
  'unknown' is the fallback (6-frame) and is always present.

  Frame policy (why it matters for the MS/MS search DB): every extra frame that isn't the real ORF
  inflates the target DB and the FDR penalty, which can bury true IDs. So emit the MINIMUM frames
  justified by annotation: defined ORFs (ncorf) -> cds=true (1 frame); genuinely unannotated
  transcript/region (ncrna, intronic, intergenic, erv) -> 6-frame discovery.

To add a brand-new origin end to end: (1) add a pattern row mapping its header to a class, and
(2) add a rule row for that class. No code change needed.
