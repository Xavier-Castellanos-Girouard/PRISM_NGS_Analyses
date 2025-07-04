#!/bin/bash

## Xavier Castellanos-Girouard
#
# Counting the reads of PRISM sgRNA libraries
#
# Date First Created: January 18th 2025
# Date Last Modified: July 4th 2025

## Change directory to upper-most directory of the GitHub repo.
cd ~/Desktop/Mammalian_CRISPR/PRISM_NGS_Analysis/

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



#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA1-Hela-Nlib/SampleA1-Hela-Nlib_S1_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hela_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA2-Hela-Clib/SampleA2-Hela-Clib_S2_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hela_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA3-HCT116-Nlib/SampleA3-HCT116-Nlib_S3_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_HCT116_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA4-HCT116-Clib/SampleA4-HCT116-Clib_S4_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_HCT116_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA5-Hap1-Nlib/SampleA5-Hap1-Nlib_S5_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hap1_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA6-Hap1-Clib/SampleA6-Hap1-Clib_S6_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hap1_Clib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA7-Vero-Clib/SampleA7-Vero-Clib_S7_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Vero_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA8-Hek293-Nlib-old/SampleA8-Hek293-Nlib-old_S8_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hek293_old_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA9-Hek293-Clib-old/SampleA9-Hek293-Clib-old_S9_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hek293_old_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA10-Hek293-Nlib-new/SampleA10-Hek293-Nlib-new_S10_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hek293_new_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA11-Hek293-Clib-new/SampleA11-Hek293-Clib-new_S11_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Hek293_new_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleA12-Nlib-Viral/SampleA12-Nlib-Viral_S12_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Viral_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleB1-Clib-Viral/SampleB1-Clib-Viral_S13_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Viral_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleB2-Nlib-Infected/SampleB2-Nlib-Infected_S14_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Infected_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleB3-Clib-Infected/SampleB3-Clib-Infected_S15_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Infected_Clib

#mageck count -l ./data/sgRNA_library_reference_files/NTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleB4-Nlib-Bac/SampleB4-Nlib-Bac_S16_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Bac_Nlib

#mageck count -l ./data/sgRNA_library_reference_files/CTERM_lib.csv --fastq ./data/DSP1645/HGKNVAFX5/fastq/Sample_SampleB5-Clib-Bac/SampleB5-Clib-Bac_S17_R1_001.fastq.gz -n ./results/MAGECK_Output/Output_Bac_Clib
