# PRISM_NGS_Analyses
Computational Pipeline for the PRISM paper

## References ##

**PRISM PAPER** Karimi M. et al. A strategy for genome-wide seamless tagging of human protein-coding genes. BiorXiv (2025).
https://doi.org/10.1101/2025.03.04.641506 

**MAGeCK** Li, W., Xu, H., Xiao, T. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol 15, 554 (2014). 
https://doi.org/10.1186/s13059-014-0554-4

**Cutadapt** Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. embnet.journal 17, (2011).
https://doi.org/10.14806/ej.17.1.200

**Bowtie2** Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). 
https://doi.org/10.1038/nmeth.1923

**Samtools** Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics  25,  2078–2079 (2009).

## Dependencies ##

- Cutadapt version 5.0
- Samtools version 1.9
- Bowtie2 version 2.3.5.1
- Python version 3.10
  - pandas 2.2.3
  - numpy 2.2.4
  - tqdm 4.67.1
- R version 4.4.0
  - dplyr 1.1.4
  - tidyr 1.3.1
  - data.table 1.15.4
  - ggplot2 3.5.1
