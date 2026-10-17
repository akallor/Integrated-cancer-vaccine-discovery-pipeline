#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

# 1. Set font globally to Arial
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# 2. Read the raw binder counts table
# Replace 'df_binder_raw_counts.tsv' with your actual file path if needed
df = pd.read_csv('df_binder_raw_counts.tsv', sep='\t')
df = df.set_index('sample')

# 3. Order categories in descending order of size (SB, WB, NB)
ordered_categories = ['SB', 'WB', 'NB']
df_ordered = df[ordered_categories]

# 4. Normalize the rows to percentage sums (0 to 100%)
# .sum(axis=1) dynamically calculates the total peptides per sample directly from this table
df_percentage = df_ordered.div(df_ordered.sum(axis=1), axis=0) * 100

# 5. Define the requested color palette
# SB = Dark Blue, WB = Cyan, NB = Teal (Updated to #009988)
colors = ["#0077BB", "#33BBEE", "#009988"]

# 6. Plot the percentage normalized stacked bar chart
ax = df_percentage.plot(kind='bar', stacked=True, figsize=(11, 7), width=0.7, color=colors)

# 7. Add percentage labels inside each stacked block
for i, column in enumerate(df_percentage.columns):
    for j, val in enumerate(df_percentage[column]):
        if val > 3.0:  # Only display labels if the block is larger than 3%
            rect = ax.containers[i].patches[j]
            x = rect.get_x() + rect.get_width() / 2
            y = rect.get_y() + rect.get_height() / 2
            
            # Print the percentage value centered inside the block
            ax.text(x, y, f"{val:.1f}%", ha='center', va='center', 
                    fontsize=9, color='white', weight='bold')

# 8. Format the final layout using Arial typography
plt.title('Percentage Normalized Distribution of Binding Classes per Sample', fontsize=14, weight='bold')
plt.xlabel('Sample', fontsize=12)
plt.ylabel('Percentage of Unique Peptides (%)', fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)
plt.ylim(0, 100)

# Position the legend on the upper right outside the plotting grid
plt.legend(title='Binding Class', labels=['Strong Binder (SB)', 'Weak Binder (WB)', 'No Binder (NB)'], 
           bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.show()
plt.savefig("mhcflurry_predicted_binders_proportion.png",dpi = 300)
plt.savefig("mhcflurry_predicted_binders_proportion.pdf")
