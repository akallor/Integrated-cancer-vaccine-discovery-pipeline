#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

# 1. Set font globally to Arial
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# 2. Read the TSV data
# Replace 'df_binder_raw_counts.tsv' with your actual file path if needed
df = pd.read_csv('df_binder_raw_counts.tsv', sep='\t')

# 3. Set 'sample' as index
df = df.set_index('sample')

# 4. Explicitly sort columns in descending order of counts (SB > WB > NB)
# This ensures that within each sample's group, the bars go from largest to smallest.
ordered_categories = ['SB', 'WB', 'NB']
df_sorted = df[ordered_categories]

# 5. Define Tol Vibrant color-blind friendly palette for the 3 categories
# SB = Dark Blue, WB = Cyan, NB = Orange
colors = ["#0077BB", "#33BBEE", "#009988"]

# 6. Plot adjacent/grouped barplot
ax = df_sorted.plot(kind='bar', figsize=(12, 7), width=0.8, color=colors)

# 7. Add text labels on top of the bars for exact counts
for container in ax.containers:
    ax.bar_label(container, fmt='%d', fontsize=8, padding=3, rotation=0)

# 8. Format layout with Arial typography
plt.title('Number of Unique Binders per Binding Class per Sample', fontsize=14, weight='bold')
plt.xlabel('Sample', fontsize=12)
plt.ylabel('Number of Unique Peptides', fontsize=12)
plt.xticks(rotation=0, ha='right', fontsize=10)
plt.yticks(fontsize=10)

# Expand y-axis slightly so labels don't get cut off at the top
plt.ylim(0, df_sorted.max().max() * 1.15)

# Style and place the legend
plt.legend(title='Binding Class', labels=['Strong Binder (SB)', 'Weak Binder (WB)', 'Non Binder (NB)'], 
           loc='upper left', frameon=True)

plt.tight_layout()
plt.show()
plt.savefig("mhcflurry_predicted_binders_raw_counts.png",dpi = 300)
plt.savefig("mhcflurry_predicted_binders_raw_counts.pdf")
