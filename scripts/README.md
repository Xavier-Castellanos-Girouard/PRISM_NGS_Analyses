# Instructions on scripts #

## Counting sgRNA reads ##
**mageck.sh** is a script that counts sgRNAs from NGS reads using the MAGeCK programme (Li et al., 2014).
The script currently runs on the NGS data for the PRISM paper (Karimi et al., 2025).
To use the script as-is, use the same file architecture as the GitHub repository. You will need to change the first line of the script (cd command) to the path towards the analysis file on your computer.
Use the fetchFastQ scripts (see data directory of this repo) to download the NGS data.

## Nested PCR Analyses ##

**NestedPCR_MasterScript.sh** is a master script to execute all the analyses related to the NestedPCR data (Figure 4 of Karimi et al., 2025).
First, the master script will run Cterm_Trimming.sh and Nterm_Trimming.sh, which use Cutadapt to remove unwanted reads (plasmid) and trim flanking regions (DHFR and linker).
Second, the master script will align the reads to a reference genome using Bowtie2. The alignments will be sorted and converted to BAM format using samtools.
Third, the master script will filter the Bowtie2 reads using the clean_Bowtie_alignments.py script. This removes the SAM header, filters unaligned reads, etc.
Fourth, the master script will check if the aligned reads are consistent with on-Target integration of the donor DNA. A gene name and accession ID will also be assigned to each read if it is onTarget.

The Final_Analysis_and_Figure.R script is meant to be opened with RStudio to view the data. We currently using a Mapping Quality Threshold (MAPQ) of 20 and use concordant and reverse-only reads to estimate on-target ratios.

## REFERENCES ##

Li, W., Xu, H., Xiao, T. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol 15, 554 (2014). 
https://doi.org/10.1186/s13059-014-0554-4


Karimi M. et al. A strategy for genome-wide seamless tagging of human protein-coding genes. BiorXiv (2025).
doi: https://doi.org/10.1101/2025.03.04.641506 
