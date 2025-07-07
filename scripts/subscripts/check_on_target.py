# Xavier Castellanos-Girouard
# 
# Assign bowtie alignments from nested PCR to sgRNAs
#
# Date First Created: April 5th 2025
# Date Last Modified: June 28th 2025


########## Import libraries ##########

import pandas as pd
import numpy as np
import sys
from tqdm import tqdm



########## Import Data ##########

## First Argument should be bowtie2 alignments formatted to clean csv files
annotated_reads_DF = pd.read_csv(sys.argv[1], index_col = False)

## Second Argument should be reference sgRNA library
sgRNA_file_DF = pd.read_csv(sys.argv[2], sep = "\t")
#sgRNA_file_DF = pd.read_csv("../../common_data/NTERM-LIBRARY-TABLE.tsv", sep = "\t")

## Get what terminus is targeted in the PRISM sample (third argument)
TERMINUS = sys.argv[3] # Should be either "Nterm" or "Cterm"


   

########## Format Data ##########

### Format read table

## Make a new column for on-target results. starts as False
annotated_reads_DF['gene_name'] = "None" # Also make column for gene name
annotated_reads_DF['accession_ID'] = "None" # Also make column for accession ID
annotated_reads_DF['onTarget'] = False

## Convert positions to integer type
annotated_reads_DF['start'] = annotated_reads_DF['start'].astype(int)
annotated_reads_DF['end'] = annotated_reads_DF['end'].astype(int)


### Format sgRNA table

## Split genomic location into chromosome and nucleotide position
sgRNA_file_DF[["chromosome", "position"]] = sgRNA_file_DF["genomic_location"].str.split(":",expand=True)

## Nucleotide position should be numeric
sgRNA_file_DF["position"] = sgRNA_file_DF["position"].astype(int)

## Change nucleotide position to cut site
# Note: Cut site is at 17th base pair of sgRNA, this is (3 nucleotides 
# 	before start of PAM site). Nucleotide position currently starts 
# 	at first sgRNA nucleotide.
sgRNA_file_DF["position"] = sgRNA_file_DF.apply(lambda x: x.position + 16 if x.strand=="+" else x.position, axis=1)
sgRNA_file_DF["position"] = sgRNA_file_DF.apply(lambda x: x.position - 17 if x.strand=="-" else x.position, axis=1)




########## Separate sgRNAs and reads according to Chromosome ##########

# Note: Separating the dataframes into chromosomes speends up the flow loop below

## Make a list of unique chromosomes targeted by our sgRNAs
unique_chr_ls = list(set(x for x in sgRNA_file_DF["chromosome"].tolist()))

## Make a dictionary of DFs, each DF is for a chromosome
sgRNA_DF_dict = {}
read_DF_dict = {}

# Sort the sgRNAs into dataframes according to their chromosome
for chrom in unique_chr_ls:
    sgRNA_DF_dict[chrom] = sgRNA_file_DF[sgRNA_file_DF["chromosome"]==chrom].reset_index(drop = True)

# Sort the reads into dataframes according to their chromosome
for chrom in unique_chr_ls:
    read_DF_dict[chrom] = annotated_reads_DF[annotated_reads_DF["chromosome"]==chrom].reset_index(drop = True)

# Delete original read dataframe to free up space
del annotated_reads_DF

########## Determine which bowtie alignments are onTarget ##########

# Iterate through each chromosome
for chrom in tqdm(unique_chr_ls):
    
    ## Get the sgRNA and read DFs for the current chromosome
    sgRNA_DF = sgRNA_DF_dict[chrom]
    read_DF = read_DF_dict[chrom]
    
    ## Get the positions of the sgRNAs, the associated gene, and accession ID
    sgRNA_pos_arr = sgRNA_DF["position"].values.astype(np.uint32)
    sgRNA_gene_arr = sgRNA_DF["gene"].values
    sgRNA_acces_arr = sgRNA_DF["accession"].values
    
    ## Get the start position of the read
    read_start_arr = read_DF['start'].values.astype(np.uint32)
    read_end_arr = read_DF['end'].values.astype(np.uint32)
    
    
    ##### CONCORDANT AND FORWARD-ONLY READS #####
    
    ## Determine lower and upper bounds for "on-target sites" 
    # For this run, we use (+/- 50 nt from sgRNA cut site)
    min_pos_arr = sgRNA_pos_arr - 50
    max_pos_arr = sgRNA_pos_arr + 50
    
    ## Check if read position is within any of the sgRNA cut-site intervals
    # Positions are compared to all lower and upper bounds.
    # Read position should be higher than min (cut site -50nt), 
    # but lower than max (cut site +50nt)
    # is_within_interval_arr is a 2D array where rows are reads, cols are
    # sgRNAs, and each individual element is TRUE/FALSE depending on whether
    # it is within the sgRNA cut-site interval. 
    is_within_min_arr = read_start_arr[:, np.newaxis] >= min_pos_arr
    is_within_max_arr = read_start_arr[:, np.newaxis] <= max_pos_arr
    is_within_intervals_arr = is_within_min_arr & is_within_max_arr
    
    # For every read, check whether it is in any sgRNA cut-site interval
    is_in_cutsite_arr = np.any(is_within_intervals_arr, axis=1)
    
    ## Make boolean indexes according to read category (proper pair or forward)
    proper_mask_arr = (read_DF['category'] == 'proper_pairs') & (is_in_cutsite_arr)
    forward_mask_arr = (read_DF['category'] == 'forward_only') & (is_in_cutsite_arr)
    
    ## Retrieve gene name and accession IDs for each read
    # This exploites the matrix structure of is_within_interval_arr
    # Where rows == read, col == sgRNA
    # np.argwhere gives you the position where the read and 
    # sgRNA match (True value in matrix)
    matches = np.array([(read_idx, sgRNA_idx) for read_idx, sgRNA_idx in np.argwhere(is_within_intervals_arr==True)])
    
    for (read_idx, sgRNA_idx) in matches:
        #read_idx = match[0]
        #sgRNA_idx = match[1]
        read_DF.loc[read_idx, "gene_name"] = sgRNA_gene_arr[sgRNA_idx]
        read_DF.loc[read_idx, "accession_ID"] = sgRNA_acces_arr[sgRNA_idx]
    
    
    ##### REVERSE-ONLY READS #####
    
        
    ## Determine lower and upper bounds for "on-target sites" 
    # For reverse only reads, we use (+/- 1000 nt from sgRNA cut site)
    min_pos_arr = sgRNA_pos_arr - 1000
    max_pos_arr = sgRNA_pos_arr + 1000
    
    ## Check if read position is within any of the sgRNA cut-site intervals
    # Note: the end position of the sequence is used, because it is closest
    # 	    to the junction.
    # Positions are compared to all lower and upper bounds.
    # Read position should be higher than min (cut site -1000nt), 
    # but lower than max (cut site +1000nt)
    # is_within_interval_arr is a 2D array where rows are reads, cols are
    # sgRNAs, and each individual element is TRUE/FALSE depending on whether
    # it is within the sgRNA cut-site interval. 
    is_within_min_arr = read_end_arr[:, np.newaxis] >= min_pos_arr
    is_within_max_arr = read_end_arr[:, np.newaxis] <= max_pos_arr
    is_within_intervals_arr = is_within_min_arr & is_within_max_arr
    
    # For every read, check whether it is in any sgRNA cut-site interval
    is_in_cutsite_arr = np.any(is_within_intervals_arr, axis=1)
    
    ## Make boolean indexes according to read category (proper pair or forward)
    reverse_mask_arr = (read_DF['category'] == 'reverse_only') & (is_in_cutsite_arr)
    
    ## Retrieve gene name and accession IDs for each read
    # This exploites the matrix structure of is_within_interval_arr
    # Where rows == read, col == sgRNA
    # np.argwhere gives you the position where the read and 
    # sgRNA match (True value in matrix)
    matches = np.array([(read_idx, sgRNA_idx) for read_idx, sgRNA_idx in np.argwhere(is_within_intervals_arr==True)])
    
    for (read_idx, sgRNA_idx) in matches:
        #read_idx = match[0]
        #sgRNA_idx = match[1]
        if read_DF.loc[read_idx, "category"] == "reverse_only":
            read_DF.loc[read_idx, "gene_name"] = sgRNA_gene_arr[sgRNA_idx]
            read_DF.loc[read_idx, "accession_ID"] = sgRNA_acces_arr[sgRNA_idx]
    
    
    ##### Assign on/off-target annotation to reads
    read_DF.loc[proper_mask_arr | forward_mask_arr | reverse_mask_arr, 'onTarget'] = True
    


########## Export ##########

annotated_reads_DF = pd.concat(list(read_DF_dict.values()))
output_filename = TERMINUS+"/results/"+TERMINUS+"_annotated_onTarget.csv"
annotated_reads_DF.to_csv(output_filename)


