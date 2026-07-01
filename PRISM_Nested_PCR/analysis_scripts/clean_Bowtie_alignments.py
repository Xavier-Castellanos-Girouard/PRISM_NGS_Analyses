# Xavier Castellanos-Girouard
# 
# Make a clean table from Bowtie2 alignments
#
# Date First Created: July 1st 2025
# Date Last Modified: November 23rd 2025

########## Import libraries ##########

import pysam
import sys
import csv
import os
from tqdm import tqdm

########## Import Data ##########

file_path = sys.argv[1]
TERMINUS = sys.argv[2] 

# Create output directory if it doesn't exist
output_dir = os.path.join(TERMINUS, "results")
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, f"{TERMINUS}_mapped_clean.csv")

########## Process BAM file ##########

bamfile = pysam.AlignmentFile(file_path, "rb")

# Open the CSV file for writing
with open(output_file, mode='w', newline='') as csvfile:
    fieldnames = [
        "query_name", "chromosome", "start", "end", 
        "mapping_quality", "is_read1", "is_read2", "flag", "category"
    ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for read in tqdm(bamfile.fetch(until_eof=True)):
        
        # Filter unmapped or QC fail
        if read.is_unmapped or read.is_qcfail:
            continue
        
        ## Coordinate calculation
        reference_span = read.reference_length

        if read.is_reverse:
            # 5' end is the right-most mapped position
            # read.reference_end points to one base AFTER the alignment
            start = read.reference_end - 1 
            end = start - reference_span
        else:
            # 5' end is the left-most mapped position
            start = read.reference_start
            end = start + reference_span
        
        ## Categorize reads according to paired or unpaired alignments
        category = "unknown"
        if read.is_proper_pair:
            category = "proper_pairs"
        elif read.is_read1 and (read.mate_is_unmapped or not read.is_proper_pair):
            category = "forward_only"
        elif read.is_read2 and (read.mate_is_unmapped or not read.is_proper_pair):
            category = "reverse_only"
        
        ## Add row to Dataframe directly on file
        writer.writerow({
            "query_name": read.query_name,
            "chromosome": read.reference_name,
            "start": start,
            "end": end,
            "mapping_quality": read.mapping_quality,
            "is_read1": read.is_read1,
            "is_read2": read.is_read2,
            "flag": read.flag,
            "category": category
        })

bamfile.close()
print(f"Processing complete. Saved to {output_file}")
