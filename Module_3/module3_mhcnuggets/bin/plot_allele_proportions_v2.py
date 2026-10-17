
#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle

# 1. Set font globally to Arial
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# 2. Define the Paul Tol Vibrant color-blind friendly palette
tol_vibrant_colors = [
    "#EE7733", # Orange
    "#0077BB", # Blue
    "#33BBEE", # Cyan
    "#EE3377", # Magenta
    "#CC3311", # Red
    "#009988", # Teal
    "#BBBBBB"  # Grey
]
color_cycler = cycle(tol_vibrant_colors)

# 3. Setup data and normalize to 100%

df_allele_pep = pd.read_csv("df_allele_pep_counts.tsv",sep = "\t")

df_plot = df_allele_pep.set_index('sample')
total_peptides = df_plot['peptide']
allele_columns = df_plot.drop(columns=['peptide'])
df_norm = allele_columns.div(total_peptides, axis=0) * 100

# Remove completely empty columns across all samples
df_norm = df_norm.loc[:, (df_norm != 0).any(axis=0)] 

# Assign a solid color from the Tol cycle to each remaining allele column
assigned_colors = [next(color_cycler) for _ in range(len(df_norm.columns))]

# 4. Plot the chart (without legend)
ax = df_norm.plot(kind='bar', stacked=True, figsize=(11, 8), width=0.7, color=assigned_colors, legend=False)

# 5. Add small allele text labels inside each stacked block
# Pandas plots columns bottom-to-top, so we track heights per sample (container)
for i, column in enumerate(df_norm.columns):
    for j, val in enumerate(df_norm[column]):
        if val > 3.5:  # Only write text if the block is large enough to read ( > 3.5% height)
            # Find the bottom baseline of the current block
            # Matplotlib stores rectangles grouped by column sequence
            rect = ax.containers[i].patches[j]
            x = rect.get_x() + rect.get_width() / 2
            y = rect.get_y() + rect.get_height() / 2
            
            # Draw the allele label inside the segment
            ax.text(x, y, column, ha='center', va='center', 
                    fontsize=7, color='white', weight='bold', rotation=90)

# 6. Format layout with Arial typography
plt.title('Percentage Normalized Proportional Distribution of Peptides per Allele', fontsize=14, weight='bold', fontname='Arial')
plt.xlabel('Sample', fontsize=12, fontname='Arial')
plt.ylabel('Percentage of Peptides (%)', fontsize=12, fontname='Arial')
plt.xticks(rotation=45, ha='right', fontsize=10, fontname='Arial')
plt.yticks(fontsize=10, fontname='Arial')
plt.ylim(0, 100)

plt.tight_layout()
plt.show()
plt.savefig('mhcflurry_predicted_binders_proportion_v2.png',dpi = 300)
