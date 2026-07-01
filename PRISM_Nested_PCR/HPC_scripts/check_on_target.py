# Xavier Castellanos-Girouard
# 
# Annotate Bowtie2 reads with sgRNA targets
#
# Date Last Modified: November 26th 2025

import pandas as pd
import numpy as np
import sys
import os
from tqdm import tqdm

# 1. SETUP AND LOAD
# ---------------------------------------------------------
read_file = sys.argv[1]
sgrna_file = sys.argv[2]
gene_table_file = sys.argv[3]
TERMINUS = sys.argv[4]        # "Nterm" or "Cterm"

print(f"--- Processing {TERMINUS} Samples ---")

# Load Reads
print(f"Loading reads from {read_file}...")
reads_df = pd.read_csv(read_file)

# Drop any reads with NaN coordinates (safety check)
reads_df = reads_df.dropna(subset=['chromosome', 'start', 'end'])
# Ensure coordinates are integers
reads_df['start'] = reads_df['start'].astype(int)
reads_df['end'] = reads_df['end'].astype(int)

# Load sgRNAs
print(f"Loading sgRNAs from {sgrna_file}...")
sgrna_df = pd.read_csv(sgrna_file, sep="\t")

# Load Gene Table
print(f"Loading Gene Annotation from {gene_table_file}...")
gene_df = pd.read_csv(gene_table_file, sep="\t")

# 2. PREPARE REFERENCE DATA
# ---------------------------------------------------------

# A. Clean Gene Table
gene_df_clean = gene_df[['name', 'chrom', 'strand', 'cdsStart', 'cdsEnd']].copy()
gene_df_clean.rename(columns={'name': 'accession', 'strand': 'gene_strand'}, inplace=True)

# B. Prepare sgRNA Library
sgrna_df[['chromosome', 'position']] = sgrna_df['genomic_location'].str.split(':', expand=True)
sgrna_df['position'] = sgrna_df['position'].astype(int)

# Calculate sgRNA Cut Sites (Cas9: 3bp upstream of PAM)
sgrna_df['cut_site'] = np.where(
    sgrna_df['strand'] == '+', 
    sgrna_df['position'] + 17, 
    sgrna_df['position'] - 17
)

# C. Merge Gene Info onto sgRNA Library
sgrna_merged = pd.merge(sgrna_df, gene_df_clean, on='accession', how='left')

# Drop sgRNAs with missing chromosome info
sgrna_merged = sgrna_merged.dropna(subset=['chromosome', 'cut_site'])
sgrna_merged['cut_site'] = sgrna_merged['cut_site'].astype(int)


# 3. PROCESSING LOOP (Split-Apply-Combine)
# ---------------------------------------------------------
# We process each chromosome separately. This guarantees 'merge_asof' 
# sees perfectly sorted data every time.

annotated_chunks = []
unique_chroms = sorted(list(set(reads_df['chromosome'].unique()) & set(sgrna_merged['chromosome'].unique())))

print(f"Processing {len(unique_chroms)} chromosomes...")

for chrom in tqdm(unique_chroms):
    
    # A. Subset Data for this Chromosome
    reads_chrom = reads_df[reads_df['chromosome'] == chrom].copy()
    sgrna_chrom = sgrna_merged[sgrna_merged['chromosome'] == chrom].copy()
    
    if reads_chrom.empty or sgrna_chrom.empty:
        continue

    # B. Sort STRICTLY (Required for merge_asof)
    # Since we are inside a chromosome loop, we don't need to sort by chrom
    sgrna_chrom = sgrna_chrom.sort_values(by='cut_site')
    
    # --- PROCESS PRECISE READS (R1) ---
    precise_mask = reads_chrom['category'].isin(['proper_pairs', 'forward_only'])
    precise_reads = reads_chrom[precise_mask].copy()
    
    if not precise_reads.empty:
        precise_reads = precise_reads.sort_values(by='start')
        
        # Merge (No 'by' argument needed here)
        merged_precise = pd.merge_asof(
            precise_reads,
            sgrna_chrom[['cut_site', 'gene', 'accession', 'cdsStart', 'cdsEnd', 'gene_strand']],
            left_on='start',
            right_on='cut_site',
            tolerance=50,
            direction='nearest'
        )
        
        merged_precise['onTarget'] = merged_precise['gene'].notna()
        
        # Calculations
        merged_precise['bp_loss'] = np.nan
        merged_precise['aa_removed'] = np.nan
        merged_precise['is_in_frame'] = False
        merged_precise['dist_to_cut_site'] = merged_precise['start'] - merged_precise['cut_site']

        valid_idx = merged_precise['onTarget']
        
        if TERMINUS == "Nterm":
            # Nterm Logic
            cond_plus = (merged_precise['gene_strand'] == '+') & valid_idx
            merged_precise.loc[cond_plus, 'bp_loss'] = merged_precise.loc[cond_plus, 'start'] - merged_precise.loc[cond_plus, 'cdsStart']
            
            cond_minus = (merged_precise['gene_strand'] == '-') & valid_idx
            merged_precise.loc[cond_minus, 'bp_loss'] = merged_precise.loc[cond_minus, 'cdsEnd'] - merged_precise.loc[cond_minus, 'start']
            
        elif TERMINUS == "Cterm":
            # Cterm Logic
            cond_plus = (merged_precise['gene_strand'] == '+') & valid_idx
            merged_precise.loc[cond_plus, 'bp_loss'] = merged_precise.loc[cond_plus, 'cdsEnd'] - merged_precise.loc[cond_plus, 'start']
            
            cond_minus = (merged_precise['gene_strand'] == '-') & valid_idx
            merged_precise.loc[cond_minus, 'bp_loss'] = merged_precise.loc[cond_minus, 'start'] - merged_precise.loc[cond_minus, 'cdsStart']

        merged_precise['aa_removed'] = (merged_precise['bp_loss'] / 3).round(2)
        merged_precise['is_in_frame'] = (merged_precise['bp_loss'].fillna(0).astype(int) % 3 == 0)
        
        annotated_chunks.append(merged_precise)


    # --- PROCESS IMPRECISE READS (R2) ---
    imprecise_mask = reads_chrom['category'] == 'reverse_only'
    imprecise_reads = reads_chrom[imprecise_mask].copy()
    
    if not imprecise_reads.empty:
        imprecise_reads = imprecise_reads.sort_values(by='end')
        
        merged_imprecise = pd.merge_asof(
            imprecise_reads,
            sgrna_chrom[['cut_site', 'gene', 'accession']],
            left_on='end',
            right_on='cut_site',
            tolerance=1000,
            direction='nearest'
        )
        
        merged_imprecise['onTarget'] = merged_imprecise['gene'].notna()
        
        # Fill NaNs
        for col in ['bp_loss', 'aa_removed', 'is_in_frame', 'dist_to_cut_site']:
            merged_imprecise[col] = np.nan
            
        annotated_chunks.append(merged_imprecise)


# 4. EXPORT
# ---------------------------------------------------------
print("Concatenating and Saving...")
if annotated_chunks:
    final_df = pd.concat(annotated_chunks, ignore_index=True)
    
    final_df['gene'] = final_df['gene'].fillna("None")
    final_df['accession'] = final_df['accession'].fillna("None")

    output_path = os.path.join(TERMINUS, "results", f"{TERMINUS}_annotated_onTarget_v4.csv")
    final_df.to_csv(output_path, index=False)
    print(f"Done. Saved to {output_path}")
else:
    print("Warning: No annotated reads found. CSV not created.")
