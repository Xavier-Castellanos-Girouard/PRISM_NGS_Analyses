#! /bin/bash

# Xavier Castellanos-Girouard
#
# v2 of NGS pipeline for nested PCR data. 
# In this version. Fwd and rev reads are treated separate whenever possible
#
# Date First Created: April 1st 2025
# Date Last Modified: June 26th 2025

## Note: Information in the fastq file sequence identifier follows the following format:
# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>

#### Libraries, modules and other software ####

# Cutadapt version 5.0
# Samtools version 1.9
# Bowtie2 version 2.3.5.1
# Bowtie version 1.2.2 (NO LONGER RELEVANT)
# htseq-count version 2.0.9 (NO LONGER RELEVANT)

########## Initiate directories ##########

# Designate main directory. This should normally be the directory in which 
# "MasterScript.sh" is found. 
# main_dir="/home/xavier/Desktop/Mammalian_CRISPR/PRISM_Nested_PCR/Xaviers_Pipeline/"

# Go to main_dir
#cd $main_dir

##### Directories for N-terminal analysis #####

mkdir ../results/Nterm

## For bowtie and htseq mapping/annotation files
mkdir ../results/Nterm/Bowtie2_Mapping # Keep bowtie2 mappings here

## For Log/summary files created at every processing step
mkdir ../results/Nterm/Logs

## Keep final results here
mkdir ../results/Nterm/results


##### Directories for C-terminal analysis #####

mkdir ../results/Cterm

## For bowtie and htseq mapping/annotation files
mkdir ../results/Cterm/Bowtie2_Mapping # Keep bowtie2 mappings here

## For Log/summary files created at every processing step
mkdir ../results/Cterm/Logs

## Keep final results here
mkdir ../results/Cterm/results


########## Trim flanking regions of reads ##########

bash ./subscripts/Nterm_Trimming.sh
bash ./subscripts/Cterm_Trimming.sh


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

bowtie2 -x ../data/Human_Genome_hg38/Bowtie2Index/genome \
  -1 ../results/Nterm/Intermediate/polyG_removal/noPolyG_NLib_NestedPCR_R1.fastq.gz \
  -2 ../results/Nterm/Intermediate/polyG_removal/noPolyG_NLib_NestedPCR_R2.fastq.gz \
  --very-sensitive-local \
  --dovetail \
  -p 4 \
  -S ../results/Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped.sam \
  -I 10 \
  -X 1000 \
  2> ../results/Nterm/Logs/Nterm_bowtie2_mapped.txt

## Convert sam files to bam files and sort according to position
samtools view -bS ../results/Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped.sam | samtools sort -o ../results/Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped_sorted.bam


##### Alignment of forward and reverse for Cterminal sample #####

bowtie2 -x ../data/Human_Genome_hg38/Bowtie2Index/genome \
  -1 ../results/Cterm/Intermediate/polyG_removal/noPolyG_CLib_NestedPCR_R1.fastq.gz \
  -2 ../results/Cterm/Intermediate/polyG_removal/noPolyG_CLib_NestedPCR_R2.fastq.gz \
  --very-sensitive-local \
   --dovetail \
  -p 4 \
  -S ../results/Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped.sam \
  -I 10 \
  -X 1000 \
  2> ../results/Cterm/Logs/Cterm_bowtie2_mapped.txt

## Convert sam files to bam files and sort according to position
samtools view -bS ../results/Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped.sam | samtools sort -o ../results/Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped_sorted.bam


########## On-target analysis and frame-shift ##########


##### N-terminal analysis #####

## Clean the bowtie2 alignment files
# Arg1 is the input BAM file containing sorted bowtie2 alignments. Arg2 is the terminus of the sample
python3.10 ./subscripts/clean_Bowtie_alignments.py ../results/Nterm/Bowtie2_Mapping/Nterm_bowtie2_mapped_sorted.bam Nterm

## Find the reads that are on-target with sgRNA library
python3.10 ./subscripts/check_on_target.py ../results/Nterm/results/Nterm_mapped_clean.csv ../data/sgRNA_library/NTERM-LIBRARY-TABLE.tsv Nterm

## Compute the nucleotide shift compared to gene start site
#python3.10 check_frame.py ./Nterm/results/Nterm_counted_onTarget.csv ../Human_Genome_hg38/hg38.gene_table.tsv Nterm


##### C-terminal analysis #####

## Clean the bowtie2 alignment files
# Arg1 is the input BAM file containing sorted bowtie2 alignments. Arg2 is the terminus of the sample
python3.10 ./subscripts/clean_Bowtie_alignments.py ../results/Cterm/Bowtie2_Mapping/Cterm_bowtie2_mapped_sorted.bam Cterm

## Find the reads that are on-target with sgRNA library
python3.10 ./subscripts/check_on_target.py ../results/Cterm/results/Cterm_mapped_clean.csv ../data/CTERM-LIBRARY-TABLE.tsv Cterm

## Compute the nucleotide shift compared to gene start site
#python3.10 check_frame.py ./Cterm/results/Cterm_counted_onTarget.csv ../Human_Genome_hg38/hg38.gene_table.tsv Cterm




