# Instructions on scripts #

## Counting sgRNA reads ##
**mageck.sh** is a script that counts sgRNAs from NGS reads using the MAGeCK programme (Li et al., 2014). <br>
<br>
The script currently runs on the NGS data for the PRISM paper (Karimi et al., 2025). <br>
<br>
To use the script as-is, use the same file architecture as the GitHub repository. You will need to change the first line of the script (cd command) to the path towards the PRISM_NGS_Analyses directory on your computer. <br>
<br>
Use the fetchFastQ scripts (see ./data directory of the repo) to download the NGS data. <br>

## REFERENCES ##

Li, W., Xu, H., Xiao, T. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol 15, 554 (2014). 
https://doi.org/10.1186/s13059-014-0554-4


Karimi M. et al. A strategy for genome-wide seamless tagging of human protein-coding genes. BiorXiv (2025).
doi: https://doi.org/10.1101/2025.03.04.641506 
