## Nested PCR Analyses ##

**NestedPCR_MasterScript.sh** is a master script to execute all the analyses related to the NestedPCR data (Figure 4 of Karimi et al., 2025). <br>
- First, the master script will run Cterm_Trimming.sh and Nterm_Trimming.sh, which use Cutadapt to remove unwanted reads (plasmid) and trim flanking regions (DHFR and linker). <br>
- Second, the master script will align the reads to a reference genome using Bowtie2. The alignments will be sorted and converted to BAM format using samtools. <br>
- Third, the master script will filter the Bowtie2 reads using the clean_Bowtie_alignments.py script. This removes the SAM header, filters unaligned reads, etc. <br>
- Fourth, the master script will check if the aligned reads are consistent with on-Target integration of the donor DNA. A gene name and accession ID will also be assigned to each read if it is onTarget. <br>
<br>
The Final_Analysis_and_Figure.R script is meant to be opened with RStudio to view the data and make figures. We currently using a Mapping Quality Threshold (MAPQ) of 20 and use concordant and reverse-only reads to estimate on-target ratios. <br>
 <br>
The scripts executed by NestedPCR_MasterScript.sh are located in the subscripts directory.
