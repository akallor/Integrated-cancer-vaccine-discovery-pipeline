#!/usr/bin/env python
"""
Peptide-sharing analysis across a grouping variable (HLA allele or tissue).

For any grouping column this produces:
  1. A heatmap of Jaccard SIMILARITY between each group's set of unique peptides
        |A n B| / |A u B|   (darker = more sharing, colour-blind-safe map)
     with raw shared-peptide counts annotated and the colour-bar on the right.
  2. A ranked table: each shared peptide and the list of groups it appears in,
     sorted by how many groups share it (most-shared first).

Run as a script to generate both the allele and tissue outputs at once.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']


# ----------------------------------------------------------------------------
# Core set logic
# ----------------------------------------------------------------------------
def build_peptide_sets(df, peptide_col, group_col):
    """group -> set(unique peptides)."""
    return {
        grp: set(sub[peptide_col].unique())
        for grp, sub in df.groupby(group_col)
    }


def compute_matrices(pep_sets, groups):
    """Return (jaccard_similarity, shared_count) square matrices."""
    n = len(groups)
    jac = np.zeros((n, n), dtype=float)
    shared = np.zeros((n, n), dtype=int)

    for i, a in enumerate(groups):
        A = pep_sets[a]
        for j, b in enumerate(groups):
            B = pep_sets[b]
            inter = len(A & B)
            union = len(A | B)
            shared[i, j] = inter
            jac[i, j] = inter / union if union else 0.0

    return jac, shared


# ----------------------------------------------------------------------------
# Ranked shared-peptide table
# ----------------------------------------------------------------------------
def rank_shared_peptides(df, peptide_col, group_col, shared_only=True):
    """
    One row per peptide:
        <peptide_col> | n_<group> | <group>_list
    ranked by the number of distinct groups the peptide appears in (desc),
    then by peptide name for a stable order.
    """
    grouped = (
        df.groupby(peptide_col)[group_col]
          .apply(lambda s: sorted(set(s)))
          .reset_index(name=f'{group_col}_list')
    )
    grouped[f'n_{group_col}'] = grouped[f'{group_col}_list'].apply(len)

    if shared_only:
        grouped = grouped[grouped[f'n_{group_col}'] > 1]

    grouped = grouped.sort_values(
        [f'n_{group_col}', peptide_col], ascending=[False, True]
    ).reset_index(drop=True)

    # Present the list as a readable, single-cell string
    grouped[f'{group_col}_list'] = grouped[f'{group_col}_list'].apply('; '.join)

    return grouped[[peptide_col, f'n_{group_col}', f'{group_col}_list']]


# ----------------------------------------------------------------------------
# Heatmap
# ----------------------------------------------------------------------------
def plot_jaccard_heatmap(
    df,
    peptide_col='peptide',
    group_col='binding_allele',
    entity_label='allele',
    order='size',            # 'size' (by peptide count desc) or 'name'
    cmap='cividis_r',        # colour-blind-safe; reversed so dark = high
    annotate=True,
    save_path='jaccard_heatmap.png',
    dpi=300,
    title=None,
):
    pep_sets = build_peptide_sets(df, peptide_col, group_col)

    if order == 'size':
        groups = sorted(pep_sets, key=lambda a: (-len(pep_sets[a]), a))
    else:
        groups = sorted(pep_sets)

    n = len(groups)
    jac, shared = compute_matrices(pep_sets, groups)

    # Mask diagonal so the colour scale reflects only between-group overlap
    jac_off = jac.copy()
    np.fill_diagonal(jac_off, np.nan)

    vmax = np.nanmax(jac_off) if n > 1 else np.nan
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad('#d9d9d9')  # diagonal colour

    cell = 0.7 if n <= 4 else 0.55
    fig_w = max(6.0, cell * n + 3.5)
    fig_h = max(5.0, cell * n + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(jac_off, cmap=cmap_obj, vmin=0.0, vmax=vmax, aspect='equal')

    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(groups, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(groups, fontsize=10)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.2)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0)

    if annotate:
        ann_fs = max(7, min(12, int(160 / n)))
        for i in range(n):
            for j in range(n):
                if i == j:
                    txt, color = f"{shared[i, j]:,}", 'black'
                else:
                    if shared[i, j] == 0:
                        continue
                    txt = f"{shared[i, j]:,}"
                    color = 'white' if (jac[i, j] / vmax) > 0.55 else 'black'
                ax.text(j, i, txt, ha='center', va='center',
                        fontsize=ann_fs, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Jaccard similarity  (|A \u2229 B| / |A \u222a B|)', fontsize=10)

    if title is None:
        title = f'Peptide Sharing Between {entity_label.capitalize()}s (Jaccard similarity)'
    ax.set_title(title, fontsize=13, weight='bold', pad=14)
    ax.text(
        1.0, -0.16,
        f'Colour: Jaccard similarity (darker = more sharing; diagonal masked, grey)\n'
        f'Numbers: peptides in common (diagonal = total unique peptides for that {entity_label})',
        transform=ax.transAxes, ha='right', va='top', fontsize=8, color='#444444'
    )

    plt.tight_layout()
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    return groups, jac, shared


# ----------------------------------------------------------------------------
# Convenience: run both allele and tissue analyses
# ----------------------------------------------------------------------------
def run_all(
    df,
    peptide_col='peptide',
    allele_col='binding_allele',
    tissue_col='tissue',
    out_prefix='',
):
    # --- Allele ---
    plot_jaccard_heatmap(
        df, peptide_col, allele_col, entity_label='allele',
        order='size', save_path=f'{out_prefix}allele_jaccard_heatmap.png',
    )
    allele_tbl = rank_shared_peptides(df, peptide_col, allele_col)
    allele_tbl.to_csv(f'{out_prefix}shared_peptides_by_allele.tsv', sep='\t', index=False)

    # --- Tissue ---
    plot_jaccard_heatmap(
        df, peptide_col, tissue_col, entity_label='tissue',
        order='name', save_path=f'{out_prefix}tissue_jaccard_heatmap.png',
        title='Peptide Sharing Between Tissue Types (Jaccard similarity)',
    )
    tissue_tbl = rank_shared_peptides(df, peptide_col, tissue_col)
    tissue_tbl.to_csv(f'{out_prefix}shared_peptides_by_tissue.tsv', sep='\t', index=False)

    return allele_tbl, tissue_tbl


def _load(path):
    sep = '\t' if path.lower().endswith(('.tsv', '.txt')) else ','
    return pd.read_csv(path, sep=sep)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', help='CSV/TSV with peptide, allele and tissue columns')
    p.add_argument('--peptide-col', default='peptide')
    p.add_argument('--allele-col', default='binding_allele')
    p.add_argument('--tissue-col', default='tissue')
    p.add_argument('--out-prefix', default='')
    args = p.parse_args()

    df = _load(args.input)
    allele_tbl, tissue_tbl = run_all(
        df,
        peptide_col=args.peptide_col,
        allele_col=args.allele_col,
        tissue_col=args.tissue_col,
        out_prefix=args.out_prefix,
    )
    print('Wrote:')
    print(f'  {args.out_prefix}allele_jaccard_heatmap.png')
    print(f'  {args.out_prefix}shared_peptides_by_allele.tsv  ({len(allele_tbl)} peptides shared by >1 allele)')
    print(f'  {args.out_prefix}tissue_jaccard_heatmap.png')
    print(f'  {args.out_prefix}shared_peptides_by_tissue.tsv  ({len(tissue_tbl)} peptides shared by >1 tissue)')
