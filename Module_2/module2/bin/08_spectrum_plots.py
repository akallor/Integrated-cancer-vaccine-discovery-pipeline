#!/usr/bin/env python3
"""08_spectrum_plots.py -- Module 2 diagnostic: annotated MS2 spectrum + XIC chromatogram
for one randomly chosen NON-CANONICAL and one CANONICAL identified peptide (per HLA class).

Picks from the identified PSMs in the Sage table (target, rank-1, spectrum_q <= --fdr), resolves
each peptide's scan in the matching mzML, and writes an annotated b/y fragment spectrum plus an
extracted-ion chromatogram. If a category has no peptide it still emits the other one. Robust to
missing files / scans / unparseable modified peptides: it shuffles candidates within a category
and moves on until one plots.

    08_spectrum_plots.py --sage results.sage.tsv --mzml-dir mzml/ --hla-class I --outdir .

Requires pyopenms + matplotlib (pyopenms imported lazily so the candidate-selection logic can be
unit-tested without it).
"""
import argparse
import os
import re
import sys
import pandas as pd

CANON, NONCANON = "canonical", "non-canonical"


def _find(df, *names, contains=None):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    if contains:
        for c in df.columns:
            if contains in c.lower():
                return c
    return None


def category(proteins):
    p = str(proteins).upper()
    if "GRP=CANON" in p:
        return CANON
    if "GRP=NONCANON" in p:
        return NONCANON
    return "unknown"


def clean_mzml_name(fn):
    fn = str(fn)
    return fn.split(".mzML")[0] + ".mzML" if ".mzML" in fn else fn


def scan_of(scannr):
    m = re.search(r"scan=(\d+)", str(scannr))
    if m:
        return int(m.group(1))
    try:
        return int(float(str(scannr)))
    except (TypeError, ValueError):
        return None


# ============================ pyOpenMS-backed plotting (lazy import) ============================
def get_theoretical_fragments(peptide_seq, max_charge):
    import pyopenms as oms
    peptide = oms.AASequence.fromString(peptide_seq)
    tsg = oms.TheoreticalSpectrumGenerator()
    p = tsg.getParameters()
    p.setValue("add_b_ions", "true"); p.setValue("add_y_ions", "true"); p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    theo = oms.MSSpectrum()
    tsg.getSpectrum(theo, peptide, 1, max_charge)
    strings = theo.getStringDataArrays()[0]
    out = []
    for i, peak in enumerate(theo):
        lab = strings[i].decode("utf-8") if isinstance(strings[i], bytes) else strings[i]
        out.append({"mz": peak.getMZ(), "label": lab})
    return out


def plot_spectrum(mzml_path, scan_number, peptide_seq, charge, out_base, frag_tol=0.02, formats=("png",)):
    import pyopenms as oms
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    exp = oms.OnDiscMSExperiment(); exp.openFile(mzml_path)
    spectrum = None
    for i in range(exp.getNrSpectra()):
        if f"scan={scan_number}" in exp.getSpectrum(i).getNativeID():
            spectrum = exp.getSpectrum(i); break
    if spectrum is None:
        print(f"   scan {scan_number} not in {os.path.basename(mzml_path)}", file=sys.stderr)
        return None
    exp_mz, exp_int = spectrum.get_peaks()
    if len(exp_mz) == 0:
        return None
    theo = get_theoretical_fragments(peptide_seq, charge)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.vlines(exp_mz, 0, exp_int, colors="gray", linewidth=0.8, alpha=0.5)
    n_ann = 0
    for frag in theo:
        diffs = [abs(x - frag["mz"]) for x in exp_mz]
        if not diffs:
            continue
        j = diffs.index(min(diffs))
        if diffs[j] <= frag_tol:
            col = "blue" if str(frag["label"]).startswith("b") else "red"
            ax.vlines(exp_mz[j], 0, exp_int[j], colors=col, linewidth=1.2)
            ax.text(exp_mz[j], exp_int[j] + max(exp_int) * 0.01, frag["label"], color=col,
                    fontsize=8, rotation=90, va="bottom", ha="center")
            n_ann += 1
    ax.set_title(f"{peptide_seq} (+{charge})  |  {os.path.basename(mzml_path)}  |  scan {scan_number}  |  {n_ann} annotated")
    ax.set_xlabel("m/z"); ax.set_ylabel("Intensity"); ax.grid(True, ls="--", alpha=0.3)
    ax.set_ylim(0, max(exp_int) * 1.15)
    ax.legend(handles=[Line2D([0], [0], color="blue", lw=1.5, label="b-ions"),
                       Line2D([0], [0], color="red", lw=1.5, label="y-ions"),
                       Line2D([0], [0], color="gray", lw=0.8, alpha=0.5, label="unassigned")], loc="upper right")
    plt.tight_layout()
    for ext in formats:
        plt.savefig(f"{out_base}_spectrum.{ext}", dpi=300 if ext == "png" else None, format=ext)
    plt.close()
    return spectrum


def plot_chromatogram(mzml_path, target_rt, peptide_seq, charge, out_base,
                      rt_window=180, mz_ppm=15, formats=("png",)):
    import pyopenms as oms
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    prec_mz = oms.AASequence.fromString(peptide_seq).getMZ(charge)
    exp = oms.OnDiscMSExperiment(); exp.openFile(mzml_path)
    rts, ints = [], []
    lo, hi = target_rt - rt_window, target_rt + rt_window
    for i in range(exp.getNrSpectra()):
        si = exp.getSpectrum(i)
        if si.getMSLevel() != 1:
            continue
        rt = si.getRT()
        if lo <= rt <= hi:
            mzs, iy = si.get_peaks()
            best = 0.0
            for mz, it in zip(mzs, iy):
                if abs(mz - prec_mz) / prec_mz * 1e6 <= mz_ppm and it > best:
                    best = it
            rts.append(rt / 60.0); ints.append(best)
    if not rts:
        print(f"   no XIC data for {peptide_seq}", file=sys.stderr)
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(rts, ints, color="purple", lw=1.5, label=f"XIC {prec_mz:.4f} m/z")
    ax.fill_between(rts, ints, color="purple", alpha=0.1)
    ax.axvline(target_rt / 60.0, color="darkorange", ls="--", lw=1.2, label=f"MS2 RT {target_rt/60.0:.2f} min")
    ax.set_title(f"XIC  |  {peptide_seq} (+{charge})")
    ax.set_xlabel("Retention time (min)"); ax.set_ylabel("Intensity"); ax.grid(True, ls="--", alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()
    for ext in formats:
        plt.savefig(f"{out_base}_chromatogram.{ext}", dpi=300 if ext == "png" else None, format=ext)
    plt.close()


# ============================ candidate selection ============================
def pick_and_plot(df, cat, mzml_dir, outdir, cls, seed, opts, plot_fn=plot_spectrum, chrom_fn=plot_chromatogram):
    cand = df[df["_cat"] == cat].drop_duplicates("_pep")
    if cand.empty:
        print(f"[spectra] no {cat} peptides identified — skipping", file=sys.stderr)
        return None
    cand = cand.sample(frac=1.0, random_state=(None if seed < 0 else seed))   # shuffle
    for _, row in cand.iterrows():
        pep = str(row["_pep"]); scan = scan_of(row["_scan"])
        if scan is None:
            continue
        try:
            charge = int(float(row["_charge"]))
        except (TypeError, ValueError):
            continue
        mp = os.path.join(mzml_dir, clean_mzml_name(row["_file"]))
        if not os.path.exists(mp):
            continue
        base = os.path.join(outdir, f"spectrum_HLA-{cls}_{cat.replace('-', '')}_"
                                    f"{re.sub(r'[^A-Za-z0-9]', '_', pep)}_scan{scan}")
        try:
            spec = plot_fn(mp, scan, pep, charge, base, frag_tol=opts["frag_tol"], formats=opts["formats"])
            if spec is None:
                continue
            chrom_fn(mp, spec.getRT(), pep, charge, base,
                     rt_window=opts["rt_window"], mz_ppm=opts["mz_ppm"], formats=opts["formats"])
            print(f"[spectra] {cat}: {pep}  (scan {scan}, {os.path.basename(mp)})")
            return pep
        except Exception as e:
            print(f"[spectra] skip {pep}: {e}", file=sys.stderr)
            continue
    print(f"[spectra] could not plot any {cat} peptide (no resolvable scan/file)", file=sys.stderr)
    return None


def load_identified(path, fdr):
    d = pd.read_csv(path, sep="\t", dtype=str)
    pep = _find(d, "peptide"); prot = _find(d, "proteins", "protein", contains="protein")
    fn = _find(d, "filename", "file", "run"); scan = _find(d, "scannr", "scan", "spectrum_id")
    ch = _find(d, "charge"); q = _find(d, "spectrum_q", "psm_q") or _find(d, contains="_q")
    lab = _find(d, "label")
    for req, nm in [(pep, "peptide"), (prot, "proteins"), (fn, "filename"), (scan, "scannr/scan"), (ch, "charge")]:
        if req is None:
            sys.exit(f"{path}: missing '{nm}' column (cols={list(d.columns)})")
    d = d.rename(columns={pep: "_pep", fn: "_file", scan: "_scan", ch: "_charge"})
    d["_cat"] = d[prot].map(category)
    if lab:
        lv = pd.to_numeric(d[lab], errors="coerce")
        d = d[(lv > 0) | lv.isna()]
    if q:
        d = d[pd.to_numeric(d[q], errors="coerce").fillna(1.0) <= fdr]
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sage", required=True, help="Sage results.sage.tsv (or rank1.sage.tsv)")
    ap.add_argument("--mzml-dir", required=True, help="directory of mzML files referenced by the Sage 'filename' column")
    ap.add_argument("--hla-class", default="I")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=0, help="shuffle seed (-1 = nondeterministic)")
    ap.add_argument("--frag-tol", type=float, default=0.02, help="fragment match tolerance (Da)")
    ap.add_argument("--rt-window", type=float, default=180, help="XIC half-window (s)")
    ap.add_argument("--mz-ppm", type=float, default=15, help="XIC precursor tolerance (ppm)")
    ap.add_argument("--formats", default="png", help="comma-separated output formats (png,svg,pdf)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    opts = dict(frag_tol=a.frag_tol, rt_window=a.rt_window, mz_ppm=a.mz_ppm,
                formats=tuple(x.strip() for x in a.formats.split(",") if x.strip()))

    d = load_identified(a.sage, a.fdr)
    nc = pick_and_plot(d, NONCANON, a.mzml_dir, a.outdir, a.hla_class, a.seed, opts)
    cn = pick_and_plot(d, CANON, a.mzml_dir, a.outdir, a.hla_class, a.seed, opts)
    if nc is None and cn is None:
        print("[spectra] nothing plotted (no identified peptides with a resolvable scan/mzML).", file=sys.stderr)
    else:
        print(f"[spectra] done: non-canonical={nc or 'none'}  canonical={cn or 'none'}  -> {a.outdir}")


if __name__ == "__main__":
    main()
