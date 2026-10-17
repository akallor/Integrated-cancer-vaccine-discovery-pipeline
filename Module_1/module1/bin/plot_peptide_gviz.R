#!/usr/bin/env Rscript
# plot_peptide_gviz.R -- IGV-style peptide visualization using Bioconductor Gviz.
#
# Alternate backend to plot_peptide_igv.py (identical inputs/idea). Requires R + Gviz
# (Bioconductor). It is NOT exercised in the pipeline's build/test environment because
# CRAN/Bioconductor are not reachable there; enable it at runtime with
#   --viz_backend gviz   (Module_1.nf)  in an environment where Gviz is installed:
#   if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#   BiocManager::install("Gviz")
#
# Usage: plot_peptide_gviz.R --tsv peptides_A.tsv [peptides_B.tsv ...] --outdir DIR [--n 2] [--seed 12]
#
# Selects up to n peptides with class in {variant, noncanonical} that carry a chr locus
# (preferring one variant + one non-canonical) and draws each on its genomic coordinates:
# a GenomeAxisTrack, the source locus, the ORF/peptide block, and -- for SNV/indel -- a
# HighlightTrack at the variant position. Writes igv_<id>.png per peptide and igv_peptides.pdf.

suppressMessages(library(Gviz))

## ---- minimal arg parsing ----
argv <- commandArgs(trailingOnly = TRUE)
getv <- function(flag, default = NULL, multi = FALSE) {
  i <- which(argv == flag)
  if (length(i) == 0) return(default)
  vals <- c(); j <- i + 1
  while (j <= length(argv) && !startsWith(argv[j], "--")) { vals <- c(vals, argv[j]); j <- j + 1 }
  if (multi) vals else vals[1]
}
tsvs   <- getv("--tsv", multi = TRUE)
outdir <- getv("--outdir", ".")
n      <- as.integer(getv("--n", "2"))
seed   <- as.integer(getv("--seed", "12"))
if (is.null(tsvs)) stop("--tsv is required")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (seed >= 0) set.seed(seed)

## ---- load + select ----
df <- do.call(rbind, lapply(tsvs, function(f) read.delim(f, stringsAsFactors = FALSE)))
has_locus <- grepl("^chr[0-9A-Za-z.]+:[0-9]+-[0-9]+", df$locus)
df <- df[df$class %in% c("variant", "noncanonical") & has_locus, ]
df <- df[!duplicated(df$id), ]
if (nrow(df) == 0) {
  writeLines("# no non-canonical/variant peptide with genomic coordinates to visualize",
             file.path(outdir, "viz_selected.tsv"))
  quit(save = "no", status = 0)
}
df <- df[sample(nrow(df)), ]
vv <- df[df$class == "variant", ]; nc <- df[df$class == "noncanonical", ]
sel <- if (n >= 2 && nrow(vv) && nrow(nc)) rbind(vv[1, ], nc[1, ]) else df[seq_len(min(n, nrow(df))), ]
write.table(sel, file.path(outdir, "viz_selected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parse_locus <- function(s) {
  m <- regmatches(s, regexec("^(chr[0-9A-Za-z.]+):([0-9]+)-([0-9]+)", s))[[1]]
  list(chrom = m[2], start = as.integer(m[3]), end = as.integer(m[4]))
}
variant_pos <- function(s) {
  m <- regmatches(s, regexec("^(chr[0-9A-Za-z.]+):([0-9]+):([ACGTNacgtn]+):([ACGTNacgtn]+)", s))[[1]]
  if (length(m) == 0) return(NULL)
  list(chrom = m[2], pos = as.integer(m[3]), ref = m[4], alt = m[5])
}

## ---- draw ----
draw_one <- function(row) {
  loc <- parse_locus(row$locus); nt <- 3 * row$aa_len
  if (identical(row$strand, "-")) { gend <- loc$end - row$orf_offset_nt; gstart <- gend - nt + 1
  } else { gstart <- loc$start + row$orf_offset_nt; gend <- gstart + nt - 1 }
  col <- if (row$class == "variant") "firebrick" else "steelblue"

  gax   <- GenomeAxisTrack()
  locusT <- AnnotationTrack(start = loc$start, width = max(1, loc$end - loc$start),
                            chromosome = loc$chrom, name = "locus", fill = "grey80",
                            col = "grey50")
  orfT  <- AnnotationTrack(start = gstart, width = max(1, gend - gstart), chromosome = loc$chrom,
                           strand = ifelse(row$strand == "-", "-", "+"),
                           id = substr(row$sequence, 1, 22), name = paste0("ORF (", row$origin, ")"),
                           fill = col, col = "black", featureAnnotation = "id",
                           fontcolor.feature = col)
  tracks <- list(gax, locusT, orfT)
  lo <- min(loc$start, gstart); hi <- max(loc$end, gend)
  v <- variant_pos(row$variant)
  if (!is.null(v)) { lo <- min(lo, v$pos); hi <- max(hi, v$pos) }
  pad <- max(5, as.integer(0.12 * (hi - lo + 1)))
  ttl <- paste0(row$id, " | ", row$origin, " | sample=", row$sample, " tissue=", row$tissue,
                " | strand ", row$strand, " frame ", row$frame, " ", row$aa_len, "aa")

  plot_it <- function() {
    if (!is.null(v)) {
      ht <- HighlightTrack(trackList = tracks, start = v$pos, width = 1, chromosome = loc$chrom,
                           col = "darkorange", fill = "#ffe0b2")
      plotTracks(ht, from = lo - pad, to = hi + pad, main = ttl, cex.main = 0.8,
                 background.title = "white", col.title = "black", col.axis = "black")
    } else {
      plotTracks(tracks, from = lo - pad, to = hi + pad, main = ttl, cex.main = 0.8,
                 background.title = "white", col.title = "black", col.axis = "black")
    }
  }
  safe <- gsub("[^A-Za-z0-9]+", "_", row$id)
  png(file.path(outdir, paste0("igv_", safe, ".png")), width = 1350, height = 420, res = 150)
  plot_it(); dev.off()
  plot_it()  # also into the open PDF device
}

pdf(file.path(outdir, "igv_peptides.pdf"), width = 9, height = 2.8)
invisible(apply(sel, 1, function(r) draw_one(as.list(r))))
dev.off()
cat(sprintf("[viz-gviz] rendered %d peptide(s)\n", nrow(sel)))
