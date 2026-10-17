#!/usr/bin/env python3
"""plot_peptide_igv.py -- IGV-style genomic visualization of translated peptides (Module 1).

Picks a small number of peptides (default 2) from the translation TSV(s) -- preferring one
SNV/variant plus one other non-canonical origin -- and draws each on its genomic coordinates:
a coordinate axis, the source locus, the ORF/peptide block (with strand arrow and a residue
ruler), and, for SNV/indel records, a lollipop at the variant position.

Genomic span of an ORF is derived from the record locus + orf_offset_nt + strand + aa_len:
  strand + : gstart = locus_start + orf_offset_nt ;      gend = gstart + 3*aa_len - 1
  strand - : gend   = locus_end   - orf_offset_nt ;      gstart = gend - 3*aa_len + 1
(For indel-derived windows the mutant length differs from reference, so coords are approximate;
the panel notes this.)

Outputs (in --outdir): one igv_<id>.png per peptide, a multi-page igv_peptides.pdf, and
viz_selected.tsv listing what was chosen. If no eligible peptide has coordinates, writes a
note and exits 0 (never fails the pipeline).
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrow, Rectangle
import pandas as pd

ELIGIBLE = {"variant", "noncanonical"}


def parse_locus(locus):
    m = re.match(r"(chr[\w.]+):(\d+)-(\d+)", str(locus))
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def variant_pos(variant):
    m = re.match(r"(chr[\w.]+):(\d+):([ACGTN]+):([ACGTN]+)", str(variant), re.IGNORECASE)
    if m:
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return None


def orf_span(row):
    loc = parse_locus(row["locus"])
    if not loc:
        return None
    chrom, lstart, lend = loc
    nt = 3 * int(row["aa_len"])
    off = int(row["orf_offset_nt"])
    if str(row["strand"]) == "-":
        gend = lend - off
        gstart = gend - nt + 1
    else:
        gstart = lstart + off
        gend = gstart + nt - 1
    return chrom, lstart, lend, gstart, gend


def pick(df, n, seed):
    df = df[df["class"].isin(ELIGIBLE)].copy()
    df = df[df["locus"].map(lambda x: parse_locus(x) is not None)]
    df = df.drop_duplicates(subset=["id"])
    if df.empty:
        return df
    rng = df.sample(frac=1.0, random_state=(None if seed < 0 else seed))  # shuffle
    var = rng[rng["class"] == "variant"].copy()
    non = rng[rng["class"] == "noncanonical"]
    # prefer a variant whose position falls inside its ORF block (nicer IGV view)
    def in_span(row):
        v = variant_pos(row.get("variant", "NA"))
        sp = orf_span(row)
        if not v or not sp:
            return 0
        _c, _ls, _le, gs, ge = sp
        return 1 if min(gs, ge) <= v[1] <= max(gs, ge) else 0
    if len(var):
        var["_in"] = var.apply(in_span, axis=1)
        var = var.sort_values("_in", ascending=False, kind="stable")
    if n >= 2 and len(var) and len(non):        # prefer one SNV/variant + one non-canonical
        chosen = [var.iloc[0], non.iloc[0]]
    else:
        chosen = [rng.iloc[i] for i in range(min(n, len(rng)))]
    return pd.DataFrame(chosen)


def draw(ax, row):
    chrom, lstart, lend, gstart, gend = orf_span(row)
    var = variant_pos(row.get("variant", "NA"))
    xs = [lstart, lend, gstart, gend] + ([var[1]] if var else [])
    lo, hi = min(xs), max(xs)
    pad = max(5, int(0.12 * (hi - lo + 1)))
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(0, 10)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)

    # locus / record extent
    ax.add_patch(Rectangle((lstart, 6.2), max(1, lend - lstart), 0.5,
                           facecolor="#d9d9d9", edgecolor="#9e9e9e", lw=0.6))
    ax.text(lstart, 6.95, f"locus {chrom}:{lstart:,}-{lend:,}"
            + (f"  gene={row['gene']}" if str(row.get("gene", "NA")) != "NA" else ""),
            fontsize=7, color="#555555")

    # ORF / peptide block
    col = "#c0392b" if row["class"] == "variant" else "#2c7fb8"
    ax.add_patch(Rectangle((gstart, 3.6), max(1, gend - gstart), 1.0,
                           facecolor=col, edgecolor="black", lw=0.8, alpha=0.9))
    # strand arrow
    amid = (gstart + gend) / 2
    dx = (gend - gstart) * 0.35 * (1 if str(row["strand"]) != "-" else -1)
    ax.add_patch(FancyArrow(amid - dx / 2, 4.1, dx, 0, width=0.02,
                            head_width=0.5, head_length=max(1, abs(dx) * 0.25),
                            length_includes_head=True, color="white"))
    pep = str(row["sequence"])
    peplab = pep if len(pep) <= 22 else pep[:19] + "..."
    ax.text(amid, 3.05, peplab, ha="center", va="top", fontsize=7.5,
            family="monospace", color=col)

    # residue ruler under the ORF (ticks every 5 aa)
    aa = int(row["aa_len"])
    for i in range(0, aa + 1, 5):
        x = gstart + 3 * i if str(row["strand"]) != "-" else gend - 3 * i
        ax.plot([x, x], [3.6, 3.4], color="#888", lw=0.5)
        ax.text(x, 2.55, str(i), ha="center", va="bottom", fontsize=5.5, color="#888")

    # SNV / indel lollipop
    if var:
        vchrom, vpos, vref, valt = var
        ax.plot([vpos, vpos], [4.6, 8.2], color="#e67e22", lw=1.0, zorder=5)
        ax.plot([vpos], [8.2], marker="o", ms=6, color="#e67e22", zorder=6)
        ax.text(vpos, 8.5, f"{vref}>{valt} @ {vpos:,}", ha="center", fontsize=6.5,
                color="#d35400")

    # coordinate axis
    ax.spines["bottom"].set_position(("data", 1.4))
    ax.tick_params(axis="x", labelsize=6)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{int(v):,}"))
    ax.set_xlabel(f"{chrom} position (bp)", fontsize=7)

    approx = " (approx: indel)" if row["origin"] == "frameshift_indel" else ""
    ax.set_title(f"{row['id']}   |   {row['origin']}   |   "
                 f"sample={row.get('sample','NA')} tissue={row.get('tissue','NA')}   |   "
                 f"strand {row['strand']}  frame {row['frame']}  {aa} aa{approx}",
                 fontsize=8, loc="left")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", nargs="+", required=True, help="one or more peptides TSVs from TRANSLATE")
    ap.add_argument("--outdir", default=".", help="output directory")
    ap.add_argument("--n", type=int, default=2, help="how many peptides to visualize (default 2)")
    ap.add_argument("--seed", type=int, default=12, help="random seed (use -1 for nondeterministic)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.concat([pd.read_csv(t, sep="\t") for t in args.tsv], ignore_index=True)
    for c in ("gene", "sample", "tissue", "source", "variant"):
        if c in df.columns:
            df[c] = df[c].fillna("NA")
    sel = pick(df, args.n, args.seed)

    if sel.empty:
        note = os.path.join(args.outdir, "viz_selected.tsv")
        with open(note, "w") as fh:
            fh.write("# no non-canonical/variant peptide with genomic coordinates to visualize\n")
        print("[viz] nothing eligible (need class=variant/noncanonical with a chr locus); wrote note.",
              file=sys.stderr)
        return

    sel.to_csv(os.path.join(args.outdir, "viz_selected.tsv"), sep="\t", index=False)
    pdf = PdfPages(os.path.join(args.outdir, "igv_peptides.pdf"))
    for _, row in sel.iterrows():
        fig, ax = plt.subplots(figsize=(9, 2.8))
        draw(ax, row)
        fig.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9]+", "_", str(row["id"]))
        fig.savefig(os.path.join(args.outdir, f"igv_{safe}.png"), dpi=150)
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()
    print(f"[viz] rendered {len(sel)} peptide(s): "
          + ", ".join(f"{r['id']}({r['origin']})" for _, r in sel.iterrows()), file=sys.stderr)


if __name__ == "__main__":
    main()
