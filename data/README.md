# Data for the PRISM analyses

## 1. Counting sgRNAs using MAGeCK

The **DSP1645** directory contains the fetchFASTQ and fetchBams scripts which allow you to download the raw NGS data. <br>
<br>
The following two .csv files contain the reference sgRNA sequences for the C-terminal and N-terminal sgRNA libraries. These are formatted such that the content is recognized by MAGeCK. <br>
./sgRNA_library_reference_files/CTERM_lib.csv and ./sgRNA_library_reference_files/CTERM_lib.csv <br>

## 2. On-Target analyses using Nested PCR data

The **DSP1374** (N-terminal sample) and **DSP1472** (C-terminal sample) directories contain the fetchFASTQ and fetchBams scripts which allow you to download the raw NGS data. <br>
<br>
The following two .csv files contain the reference sgRNA sequences for the C-terminal and N-terminal sgRNA libraries. These also contain other information pertaining the the sgRNAs, such as the target gene name and accession number. <br>
./sgRNA_library_reference_files/CTERM-LIBRARY-TABLE.tsv and ./sgRNA_library_reference_files/NTERM-LIBRARY-TABLE.tsv <br>

The Human_Genome_hg38 directory should contain genome assemblies for Bowtie2 alignements. The .bt2 index files can be found at Illumina iGenomes: https://support.illumina.com/sequencing/sequencing_software/igenome.html <br>
You should specifically download the file Homo_sapiens_UCSC_hg38.tar.gz (Homo sapiens ; UCSC ; hg38).

The hg38.gene_table.tsv table in Human_Genome_hg38 was taken from Louis Gauthier's GitHub HumanInteractome repository (in the annotations directory). This contains the gene information that was used for the initial library design. We use it to annotated our reads in the check_ontarget.py python script.
