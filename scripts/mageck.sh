#!/bin/bash

## Xavier Castellanos-Girouard
#
# Counting the reads of PRISM sgRNA libraries
#
# Date First Created: January 18th 2025
# Date Last Modified: July 4th 2025

## Change directory to upper-most directory of the GitHub repo.
cd ~/Desktop/Mammalian_CRISPR/PRISM_NGS_Analysis/

## Make the directory for MAGECK output
mkdir ./results/MAGECK_Ouput/

## Make an array containing info on the NGS samples.
# Each line: sample_name, lib_type, sub_directory, sample_id
samples=(
    "Hela_Nlib NTERM SampleA1-Hela-Nlib S1"
    "Hela_Clib CTERM SampleA2-Hela-Clib S2"
    "HCT116_Nlib NTERM SampleA3-HCT116-Nlib S3"
    "HCT116_Clib CTERM SampleA4-HCT116-Clib S4"
    "Hap1_Nlib NTERM SampleA5-Hap1-Nlib S5"
    "Hap1_Clib CTERM SampleA6-Hap1-Clib S6"
    "Vero_Clib CTERM SampleA7-Vero-Clib S7"
    "Hek293_old_Nlib NTERM SampleA8-Hek293-Nlib-old S8"
    "Hek293_old_Clib CTERM SampleA9-Hek293-Clib-old S9"
    "Hek293_new_Nlib NTERM SampleA10-Hek293-Nlib-new S10"
    "Hek293_new_Clib CTERM SampleA11-Hek293-Clib-new S11"
    "Viral_Nlib NTERM SampleA12-Nlib-Viral S12"
    "Viral_Clib CTERM SampleB1-Clib-Viral S13"
    "Infected_Nlib NTERM SampleB2-Nlib-Infected S14"
    "Infected_Clib CTERM SampleB3-Clib-Infected S15"
    "Bac_Nlib NTERM SampleB4-Nlib-Bac S16"
    "Bac_Clib CTERM SampleB5-Clib-Bac S17"
)

## Loop through the samples and execute mageck-count (counts sgRNAs)
for entry in "${samples[@]}"; do
    read name lib subdir s_number <<< "$entry"
    mageck count \
        -l "./data/sgRNA_library_reference_files/${lib}_lib.csv" \
        --fastq "./data/DSP1645/HGKNVAFX5/fastq/Sample_${subdir}/${subdir}_${s_number}_R1_001.fastq.gz" \
        -n "./results/MAGECK_Output/Output_${name}"
done

