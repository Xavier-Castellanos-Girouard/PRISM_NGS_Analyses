import pandas as pd
import sys
import os

# Usage: python combine_gene_lists.py <Nterm_list.txt> <Cterm_list.txt> <Output.csv>

if len(sys.argv) < 4:
    print("Usage: python combine_gene_lists.py <Nterm_list.txt> <Cterm_list.txt> <Output.csv>")
    sys.exit(1)

nterm_file = sys.argv[1]
cterm_file = sys.argv[2]
output_file = sys.argv[3]

def load_genes(filepath):
    try:
        with open(filepath, 'r') as f:
            # Read non-empty lines and strip whitespace
            return sorted(list(set([line.strip() for line in f if line.strip()])))
    except FileNotFoundError:
        print(f"Error: File not found {filepath}")
        return []

print(f"Reading N-term list: {nterm_file}")
n_genes = load_genes(nterm_file)

print(f"Reading C-term list: {cterm_file}")
c_genes = load_genes(cterm_file)

n_set = set(n_genes)
c_set = set(c_genes)

# 1. Union: Genes in N or C (The comprehensive list)
union_genes = sorted(list(n_set | c_set))

# 2. Intersection: Genes in N and C (Detected in both)
intersection_genes = sorted(list(n_set & c_set))

# 3. Prepare data for DataFrame (Pad short lists with empty strings)
# We find the longest list to ensure the DataFrame is large enough
max_len = max(len(n_genes), len(c_genes), len(union_genes), len(intersection_genes))

def pad_list(lst, length):
    return lst + [""] * (length - len(lst))

data = {
    "N_terminal_Only": pad_list(n_genes, max_len),
    "C_terminal_Only": pad_list(c_genes, max_len),
    "Union (Any Detected)": pad_list(union_genes, max_len),
    "Intersection (Both Detected)": pad_list(intersection_genes, max_len)
}

df = pd.DataFrame(data)

# Export
df.to_csv(output_file, index=False)

print("-" * 30)
print(f"Comparison Complete.")
print(f"Saved to: {output_file}")
print("-" * 30)
print(f"Stats:")
print(f"  N-term detected: {len(n_set)}")
print(f"  C-term detected: {len(c_set)}")
print(f"  Total Unique Genes (Union): {len(union_genes)}")
print(f"  Shared Genes (Intersection): {len(intersection_genes)}")
print("-" * 30)
