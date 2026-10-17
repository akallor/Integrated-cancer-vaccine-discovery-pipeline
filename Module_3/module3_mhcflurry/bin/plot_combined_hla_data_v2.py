#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==============================================================================
# GLOBAL SETUP
# ==============================================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

colors_binary = ["#0077BB", "#33BBEE"]


def clean_and_load_dataframe(file_path, index_col_name, expected_cols):
    """
    Return a WIDE dataframe:
        index   = sample / tissue
        columns = expected_cols  (e.g. ['canonical', 'non_canonical'])

    Handles both layouts automatically:
      * WIDE : one row per sample with 'canonical'/'non_canonical' columns.
      * LONG : one row per (sample, origin) pair plus a count column -> pivoted.
    """
    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} not found. Skipping these figures.")
        return None

    df = pd.read_csv(file_path, sep='\t')

    # Drop a leading unnamed column ONLY if it's a throwaway 0..N-1 row index
    first = df.columns[0]
    if str(first).startswith('Unnamed'):
        col0 = pd.to_numeric(df[first], errors='coerce')
        if col0.notna().all() and list(col0.astype(int)) == list(range(len(df))):
            df = df.iloc[:, 1:]

    origin_set = set(expected_cols)

    # ---- Case 1: already wide ----------------------------------------------
    if origin_set.issubset(df.columns):
        idx_candidates = [c for c in df.columns if c not in origin_set]
        idx_col = idx_candidates[0] if idx_candidates else index_col_name
        wide = df.set_index(idx_col)[list(expected_cols)]

    # ---- Case 2: long / tidy -> pivot to wide ------------------------------
    else:
        # origin column: the one whose values ARE the expected origin labels
        origin_col = next(
            (c for c in df.columns
             if df[c].notna().any()
             and set(df[c].dropna().astype(str)).issubset(origin_set)),
            None
        )
        if origin_col is None:
            raise ValueError(
                f"{file_path}: found neither {expected_cols} as columns nor an "
                f"origin column containing those labels. Columns: {list(df.columns)}"
            )

        # count column: first numeric column that isn't the origin column
        count_col = next(
            c for c in df.columns
            if c != origin_col and pd.to_numeric(df[c], errors='coerce').notna().any()
        )

        # index column: whatever is left over
        idx_col = next(c for c in df.columns if c not in (origin_col, count_col))

        wide = df.pivot_table(
            index=idx_col,
            columns=origin_col,
            values=count_col,
            aggfunc='sum'
        ).reindex(columns=list(expected_cols))   # force canonical, non_canonical order

    wide.index.name = index_col_name
    return wide.apply(pd.to_numeric, errors='coerce')


def log_transform(df):
    """
    Apply log10(count + 1).

    NaN values remain NaN.
    Zero values become 0.
    """
    return np.log10(df + 1)


def make_count_labels(df):
    """
    Create labels from the ORIGINAL count dataframe.

    NaN values are given an empty label instead of attempting
    to convert NaN to an integer.
    """
    labels = []

    for value in df:
        if pd.isna(value):
            labels.append('')
        else:
            labels.append(f"{int(value):,}")

    return labels


# ==============================================================================
# DATASET 1: Sample Origin Data
# ==============================================================================
print("Processing Sample Origin Data...")

df_origin = clean_and_load_dataframe(
    file_path='df_origin_grp.tsv',
    index_col_name='sample',
    expected_cols=['canonical', 'non_canonical']
)

if df_origin is not None:

    # --------------------------------------------------------------------------
    # Figure 1: Log10-Transformed Counts
    # --------------------------------------------------------------------------
    print(" -> Generating Figure 1: Sample Origin Log-Transformed Counts...")

    df_origin_log = log_transform(df_origin)

    fig1, ax1 = plt.subplots(figsize=(12, 7))

    df_origin_log.plot(
        kind='bar',
        width=0.8,
        color=colors_binary,
        ax=ax1
    )

    # Original count labels, NaN-safe
    for i, container in enumerate(ax1.containers):

        labels = make_count_labels(df_origin.iloc[:, i])

        ax1.bar_label(
            container,
            labels=labels,
            fontsize=8,
            padding=3,
            rotation=0
        )

    ax1.set_title(
        'Log10-Transformed Counts of Canonical and Non-Canonical Peptides per Sample',
        fontsize=14,
        weight='bold'
    )

    ax1.set_xlabel('Sample', fontsize=12)

    ax1.set_ylabel(
        'Log10(Number of Unique Peptides + 1)',
        fontsize=12
    )

    ax1.tick_params(axis='x', rotation=0)

    ax1.legend(
        title='Origin',
        labels=['Canonical', 'Non-Canonical'],
        loc='upper right'
    )

    ax1.text(
        0.99,
        0.02,
        'Bars: log10(count + 1)\nLabels: original counts',
        transform=ax1.transAxes,
        ha='right',
        va='bottom',
        fontsize=9
    )

    plt.tight_layout()

    plt.savefig(
        'fig1_sample_origin_log10.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        'fig1_sample_origin_log10.pdf',
        dpi=300,
        bbox_inches='tight')


    plt.close()


    # --------------------------------------------------------------------------
    # Figure 2: Percentage-Normalized Stacked Barplot
    # --------------------------------------------------------------------------
    print(" -> Generating Figure 2: Sample Origin Percentages...")

    df_origin_pct = (
        df_origin.div(df_origin.sum(axis=1), axis=0) * 100
    )

    fig2, ax2 = plt.subplots(figsize=(11, 7))

    df_origin_pct.plot(
        kind='bar',
        stacked=True,
        width=0.7,
        color=colors_binary,
        ax=ax2
    )

    for i, column in enumerate(df_origin_pct.columns):

        for j, val in enumerate(df_origin_pct[column]):

            if pd.notna(val) and val > 1.5:

                rect = ax2.containers[i].patches[j]

                x = rect.get_x() + rect.get_width() / 2
                y = rect.get_y() + rect.get_height() / 2

                ax2.text(
                    x,
                    y,
                    f"{val:.1f}%",
                    ha='center',
                    va='center',
                    fontsize=8,
                    color='white',
                    weight='bold'
                )

    ax2.set_title(
        'Percentage Normalized Distribution of Peptide Origins per Sample',
        fontsize=14,
        weight='bold'
    )

    ax2.set_xlabel('Sample', fontsize=12)

    ax2.set_ylabel(
        'Percentage of Unique Peptides (%)',
        fontsize=12
    )

    ax2.set_ylim(0, 100)

    ax2.tick_params(axis='x', rotation=0)

    ax2.legend(
        title='Origin',
        labels=['Canonical', 'Non-Canonical'],
        bbox_to_anchor=(1.02, 1),
        loc='upper left'
    )

    plt.tight_layout()

    plt.savefig(
        'fig2_sample_origin_normalized.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        'fig2_sample_origin_normalized.pdf',
        dpi=300,
        bbox_inches='tight'
    )


    plt.close()


# ==============================================================================
# DATASET 2: Tissue Type Data
# ==============================================================================
print("\nProcessing Tissue Type Data...")

df_tissue = clean_and_load_dataframe(
    file_path='df_tissue_grp.tsv',
    index_col_name='tissue',
    expected_cols=['canonical', 'non_canonical']
)

if df_tissue is not None:

    # --------------------------------------------------------------------------
    # Figure 3: Log10-Transformed Counts
    # --------------------------------------------------------------------------
    print(" -> Generating Figure 3: Tissue Log-Transformed Counts...")

    df_tissue_log = log_transform(df_tissue)

    fig3, ax3 = plt.subplots(figsize=(8, 6))

    df_tissue_log.plot(
        kind='bar',
        width=0.6,
        color=colors_binary,
        ax=ax3
    )

    # Original count labels, NaN-safe
    for i, container in enumerate(ax3.containers):

        labels = make_count_labels(df_tissue.iloc[:, i])

        ax3.bar_label(
            container,
            labels=labels,
            fontsize=9,
            padding=3
        )

    ax3.set_title(
        'Log10-Transformed Counts of Canonical and Non-Canonical Peptides per Tissue Status',
        fontsize=14,
        weight='bold'
    )

    ax3.set_xlabel(
        'Tissue of Origin',
        fontsize=12
    )

    ax3.set_ylabel(
        'Log10(Number of Unique Peptides + 1)',
        fontsize=12
    )

    ax3.tick_params(axis='x', rotation=0)

    ax3.legend(
        title='Origin',
        labels=['Canonical', 'Non-Canonical'],
        loc='upper right'
    )

    ax3.text(
        0.99,
        0.02,
        'Bars: log10(count + 1)\nLabels: original counts',
        transform=ax3.transAxes,
        ha='right',
        va='bottom',
        fontsize=9
    )

    plt.tight_layout()

    plt.savefig(
        'fig3_tissue_log10.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        'fig3_tissue_log10.pdf',
        dpi=300,
        bbox_inches='tight'
    )

	

    plt.close()


    # --------------------------------------------------------------------------
    # Figure 4: Percentage-Normalized Stacked Barplot
    # --------------------------------------------------------------------------
    print(" -> Generating Figure 4: Tissue Percentages...")

    df_tissue_pct = (
        df_tissue.div(df_tissue.sum(axis=1), axis=0) * 100
    )

    fig4, ax4 = plt.subplots(figsize=(8, 6))

    df_tissue_pct.plot(
        kind='bar',
        stacked=True,
        width=0.5,
        color=colors_binary,
        ax=ax4
    )

    for i, column in enumerate(df_tissue_pct.columns):

        for j, val in enumerate(df_tissue_pct[column]):

            if pd.notna(val) and val > 0.5:

                rect = ax4.containers[i].patches[j]

                x = rect.get_x() + rect.get_width() / 2
                y = rect.get_y() + rect.get_height() / 2

                ax4.text(
                    x,
                    y,
                    f"{val:.1f}%",
                    ha='center',
                    va='center',
                    fontsize=9,
                    color='white',
                    weight='bold'
                )

    ax4.set_title(
        'Percentage Normalized Distribution of Peptide Origins per Tissue Status',
        fontsize=14,
        weight='bold'
    )

    ax4.set_xlabel(
        'Tissue of Origin',
        fontsize=12
    )

    ax4.set_ylabel(
        'Percentage of Unique Peptides (%)',
        fontsize=12
    )

    ax4.set_ylim(0, 100)

    ax4.tick_params(axis='x', rotation=0)

    ax4.legend(
        title='Origin',
        labels=['Canonical', 'Non-Canonical'],
        bbox_to_anchor=(1.02, 1),
        loc='upper left'
    )

    plt.tight_layout()

    plt.savefig(
        'fig4_tissue_normalized.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        'fig4_tissue_normalized.pdf',
        dpi=300,
        bbox_inches='tight'
    )


    plt.close()


print("\nAll plots generated successfully!")
print("Generated files:")
print("  - fig1_sample_origin_log10.png")
print("  - fig2_sample_origin_normalized.png")
print("  - fig3_tissue_log10.png")
print("  - fig4_tissue_normalized.png")
