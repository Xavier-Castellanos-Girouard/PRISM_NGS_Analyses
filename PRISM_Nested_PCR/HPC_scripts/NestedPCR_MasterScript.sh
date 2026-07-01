#! /bin/bash

# Xavier Castellanos-Girouard
#
# v3 of NGS pipeline for nested PCR data. 
# In this version. Fwd and rev reads are treated separate whenever possible
#
# Date First Created: November 23rd 2025
# Date Last Modified: November 23rd 2025

## Note: Information in the fastq file sequence identifier follows the following format:
# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>

#### Libraries, modules and other software ####

# Cutadapt version 5.0
# Samtools version 1.9
# Bowtie2 version 2.3.5.1
# Bowtie version 1.2.2

########## Initiate directories ##########

# Designate main directory. This should normally be the directory in which 
# "MasterScript.sh" is found. 
# main_dir="/home/xavier/Desktop/Mammalian_CRISPR/PRISM_Nested_PCR/Xaviers_Pipeline/"

# Go to main_dir
#cd $main_dir

##### Directories for N-terminal analysis #####

mkdir Nterm

## For bowtie and htseq mapping/annotation files
mkdir Nterm/Bowtie2_Mapping # Keep bowtie2 mappings here

## For Log/summary files created at every processing step
mkdir Nterm/Logs

## Keep final results here
mkdir Nterm/results


##### Directories for C-terminal analysis #####

mkdir Cterm

## For bowtie and htseq mapping/annotation files
mkdir Cterm/Bowtie2_Mapping # Keep bowtie2 mappings here

## For Log/summary files created at every processing step
mkdir Cterm/Logs

## Keep final results here
mkdir Cterm/results


########## Trim flanking regions of reads ##########

bash Nterm_Trimming.sh
bash Cterm_Trimming.sh


########## Aligning reads to the genome using Bowtie2 ##########

echo "Alignment of reads to the genome using Bowtie2"

# Note: We do the mapping for Forward and Reverse strands separately,
# 	this is to avoid downstream issues with annotation in htseq
# Note: We are not using --very-sensitive-local for now.

##### Alignment of forward and reverse reads
# Note: x argument is name of the file for human genome bowtie2 index
#  	 U argument is unpaired reads to be aligned
# 	 --local argument is used as opposed to "end-to-end" because our trimming is imperfect
# 	 -p argument is for number of cores to use
# 	 -S write to SAM alignments
# 	 -- dovetail allows reads to extend past each other
#       -I and -X indicates the minimum and maximum distance between read pairs



##### Alignment of forward and reverse for Nterminal sample #####

bowtie2 -x ../Human_Genome_hg38/Bowtie2Index/genome \
  -1 Nterm/Intermediate/polyG_removal/noPolyG_NLib_NestedPCR_R1.fastq.gz \
  -2 Nterm/Intermediate/polyG_removal/noPolyG_NLib_NestedPCR_R2.fastq.gz \
  --local \
  --very-sensitive-local \
  -p 8 \
  -S ./Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped.sam \
  --no-unal \
  -X 1000 \
  2> ./Nterm/Logs/Nterm_bowtie2_mapped.txt

## Convert sam files to bam files and sort according to position
samtools view -bS Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped.sam | samtools sort -o Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped_sorted.bam


##### Alignment of forward and reverse for Cterminal sample #####

bowtie2 -x ../Human_Genome_hg38/Bowtie2Index/genome \
  -1 Cterm/Intermediate/polyG_removal/noPolyG_CLib_NestedPCR_R1.fastq.gz \
  -2 Cterm/Intermediate/polyG_removal/noPolyG_CLib_NestedPCR_R2.fastq.gz \
  --local \
  --very-sensitive-local \
  -p 8 \
  -S ./Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped.sam \
  --no-unal \
  -X 1000 \
  2> ./Cterm/Logs/Cterm_bowtie2_mapped.txt

## Convert sam files to bam files and sort according to position
samtools view -bS Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped.sam | samtools sort -o Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped_sorted.bam


########## On-target analysis and frame-shift ##########


##### N-terminal analysis #####

## Clean the bowtie2 alignment files
python3.10 clean_Bowtie_alignments.py Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped_sorted.bam Nterm

## Find the reads that are on-target with sgRNA library
python3.10 check_on_target.py ./Nterm/results/Nterm_mapped_clean.csv ../../common_data/NTERM-LIBRARY-TABLE.tsv ../Human_Genome_hg38/hg38.gene_table.tsv Nterm

## Compute the nucleotide shift compared to gene start site
#python3.10 check_frame.py ./Nterm/results/Nterm_counted_onTarget.csv ../Human_Genome_hg38/hg38.gene_table.tsv Nterm


##### C-terminal analysis #####

## Clean the bowtie2 alignment files
python3.10 clean_Bowtie_alignments.py Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped_sorted.bam Cterm

## Find the reads that are on-target with sgRNA library
python3.10 check_on_target.py ./Cterm/results/Cterm_mapped_clean.csv ../../common_data/CTERM-LIBRARY-TABLE.tsv ../Human_Genome_hg38/hg38.gene_table.tsv Cterm

## Compute the nucleotide shift compared to gene start site
#python3.10 check_frame.py ./Cterm/results/Cterm_counted_onTarget.csv ../Human_Genome_hg38/hg38.gene_table.tsv Cterm




