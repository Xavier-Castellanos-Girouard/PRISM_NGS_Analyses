# Xavier Castellanos-Girouard
# 
# Make a clean table from annotations created with HTseq
#
# Date First Created: July 1st 2025
# Date Last Modified: July 5th 2025

########## Import libraries ##########

import pandas as pd
import pysam
import sys
from tqdm import tqdm


########## Import Data ##########

### First argument should be path with Bowtie2 alignment
file_path = sys.argv[1]
#file_path = "./Nterm/htseq_annotation/Nterm_annotated_R1.tsv"

## Get what terminus is targeted in the PRISM sample (third argument)
TERMINUS = sys.argv[2] # Should be either "Nterm" or "Cterm"

########## Process BAM file ##########

# Open sorted BAM file
bamfile = pysam.AlignmentFile(file_path, "rb")

# Initiate the lists of reads
proper_pairs = []
forward_only = []
reverse_only = []

#
for read in tqdm(bamfile.fetch(until_eof=True)):
    
    # Remove reads that are unmapped or failed the quality control filter
    if read.is_unmapped or read.is_qcfail:# or read.is_secondary or read.is_supplementary:
        continue
    
    ## Adjust the position column of the sam file.
    # In bowtie2, the position is always the left-most genomic location
    # of a read. We want to change this: the position should be the genomic
    # location of the first nucleotide of the read that aligns to the genome.
    # Determine the 5' position on the genome
    read_positions = read.get_reference_positions()
    reference_span = read.reference_length  # Sum of M, D, N, =, X in CIGAR

    if read.is_reverse:
        # 5' end is at the right-most mapped position
        start = read_positions[-1] if read_positions else read.reference_start
        end = start - reference_span
    else:
        # 5' end is at the left-most mapped position
        start = read.reference_start
        end = start + reference_span
    
    
    data = {
        "query_name": read.query_name,
        "chromosome": read.reference_name,
        "start": start,
        "end": end,
        "mapping_quality": read.mapping_quality,
        "is_read1": read.is_read1,
        "is_read2": read.is_read2,
        "flag": read.flag
    }

    ## Sort into categories according to concordance of FWD (R1) and REV (R2) reads.
    # A proper_pair means that both Read 1 (R1) and Read 2 (R2) align and
    # are facing each other in the correct direction (concordant).
    # forward_only means that the FWD read properly aligned, but not the REV or
    # that the FWD and REV reads were discordant. 
    # reverse_only means that the REV read properly aligned, but not the FWD or
    # that the FWD and REV reads were discordant. 
    if read.is_proper_pair:
        proper_pairs.append(data)
    elif read.is_read1 and (read.mate_is_unmapped or not read.is_proper_pair):
        forward_only.append(data)
    elif read.is_read2 and (read.mate_is_unmapped or not read.is_proper_pair):
        reverse_only.append(data)

bamfile.close()

# Convert to DataFrames
proper_pairs_df = pd.DataFrame(proper_pairs)
forward_only_df = pd.DataFrame(forward_only)
reverse_only_df = pd.DataFrame(reverse_only)


# Add category labels
proper_pairs_df["category"] = "proper_pairs"
forward_only_df["category"] = "forward_only"
reverse_only_df["category"] = "reverse_only"

# Concatenate all into one big DataFrame
all_reads_df = pd.concat([proper_pairs_df, forward_only_df, reverse_only_df], ignore_index=True)


########## Export data ##########

# Optional: save the combined table
all_reads_df.to_csv(TERMINUS+"/results/"+TERMINUS+"_mapped_clean.csv", index=False)

#proper_pairs_df.to_csv("Nterm/results/"+TERMINUS+"_proper_pairs.csv")
#forward_only_df.to_csv("Nterm/results/"+TERMINUS+"_forward_only.csv")
#reverse_only_df.to_csv("Nterm/results/"+TERMINUS+"_reverse_only.csv")


