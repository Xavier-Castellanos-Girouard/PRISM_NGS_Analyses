import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import numpy as np

# Usage: python generate_report.py <csv_path> <terminus>
csv_path = sys.argv[1]
terminus = sys.argv[2]
output_dir = os.path.dirname(csv_path)

print(f"Loading {csv_path}...")
df = pd.read_csv(csv_path)

# =========================================================
# 1. EXPORT UNIQUE GENES (Using ALL Reads)
# =========================================================
# We use ALL on-target reads (including R2-only) to maximize gene discovery.
df_all_target = df[df['onTarget'] == True]

unique_genes_detected = df_all_target[df_all_target['gene'].notna() & (df_all_target['gene'] != "None")]['gene'].unique()
unique_genes_detected = sorted(unique_genes_detected)
num_detected = len(unique_genes_detected)

gene_list_path = os.path.join(output_dir, f"{terminus}_detected_genes.txt")
with open(gene_list_path, "w") as f:
    for gene in unique_genes_detected:
        f.write(f"{gene}\n")

print(f"Gene list saved. Total Unique Genes Detected: {num_detected}")


# =========================================================
# 2. FILTER FOR GRAPHS (Precise Reads Only)
# =========================================================
df_precise = df[
    (df['onTarget'] == True) & 
    (df['category'].isin(['forward_only', 'proper_pairs']))
].copy()

total_reads_for_graph = len(df_precise)
# Calculate unique genes specifically for the Pie Chart data
genes_in_pie = df_precise['gene'].nunique()

print(f"Total Reads for Profiling: {total_reads_for_graph}")
print(f"Genes represented in Profile: {genes_in_pie}")


# =========================================================
# 3. DATA CLASSIFICATION
# =========================================================
LARGE_DEL_LIMIT = 60

def classify_read(row):
    loss = row['bp_loss']
    #if loss > LARGE_DEL_LIMIT:
    #    return "Large deletion"
    if loss == 0:
        return "In-frame as designed"
    if row['is_in_frame']:
        return "In-frame with indel"
    return "Out-of-frame"

df_precise['Category'] = df_precise.apply(classify_read, axis=1)
category_counts = df_precise['Category'].value_counts()


# =========================================================
# 4. GENERATE PIE CHART
# =========================================================
plt.figure(figsize=(10, 8))

blue_colors = {
    "In-frame as designed": "#08306b",
    "In-frame with indel":  "#2171b5",
    "Out-of-frame":         "#6baed6"#,
    #"Large deletion":       "#c6dbef"
}

sorted_index = category_counts.index
plot_colors = [blue_colors.get(l, "#999999") for l in sorted_index]

def smart_autopct(pct):
    return ('%1.1f%%' % pct) if pct > 5 else ''

total_counts = category_counts.sum()
legend_labels = [f'{idx} ({ (val/total_counts*100):.1f}%)' for idx, val in category_counts.items()]

plt.pie(
    category_counts, 
    labels=None, 
    autopct=smart_autopct, 
    startangle=90, 
    colors=plot_colors,
    textprops={'fontsize': 14, 'weight': 'bold', 'color': 'white'},
    wedgeprops={"edgecolor":"k",'linewidth': 0.5, 'antialiased': True}
)

plt.title(f"Edit Profile: Reads ({terminus})", fontsize=18)

# UPDATED: Added gene count to Legend Title
plt.legend(
    legend_labels, 
    title=f"Reads (n={total_reads_for_graph:,} | {genes_in_pie:,} genes)",
    title_fontsize=13,
    loc="center left", 
    bbox_to_anchor=(1, 0.5), 
    fontsize=12
)

pie_output = os.path.join(output_dir, f"{terminus}_Profile_PieChart.png")
plt.tight_layout()
plt.savefig(pie_output, dpi=300, bbox_inches='tight')
print(f"Pie Chart saved to: {pie_output}")


# =========================================================
# 5. GENERATE HISTOGRAM (All In-Frame)
# =========================================================
df_inframe = df_precise[
    (df_precise['is_in_frame'] == True)
].copy()

total_inframe = len(df_inframe)
genes_in_hist = df_inframe['gene'].nunique()

# Invert sign (Negative = Deletion)
df_inframe['aa_change'] = -1 * df_inframe['aa_removed']

AA_VIEW_LIMIT = 15
viz_df = df_inframe[df_inframe['aa_change'].between(-AA_VIEW_LIMIT, AA_VIEW_LIMIT)].copy()

plt.figure(figsize=(12, 7))

sns.histplot(
    data=viz_df,
    x='aa_change',
    binwidth=1,
    color='#2171b5', 
    edgecolor="black",
    linewidth=0.5,
    discrete=True
)

plt.title(f"In-Frame Indel Distribution ({terminus})\n(n = {total_inframe:,} reads | {genes_in_hist:,} genes represented)", fontsize=16)
plt.xlabel("Amino Acid Change\n(- = Removed, + = Added)", fontsize=12)
plt.ylabel("Read Count", fontsize=12)

plt.xlim(-AA_VIEW_LIMIT, AA_VIEW_LIMIT)
plt.xticks(range(-AA_VIEW_LIMIT, AA_VIEW_LIMIT + 1, 1))

plt.axvline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

hist_output = os.path.join(output_dir, f"{terminus}_AA_Change_Distribution.png")
plt.tight_layout()
plt.savefig(hist_output, dpi=300)

print(f"Histogram saved to: {hist_output}")
print("Done.")
