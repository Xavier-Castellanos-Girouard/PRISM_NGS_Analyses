The Human_Genome_hg38 directory should contain genome assemblies for Bowtie2 alignements. The .bt2 index files can be found at Illumina iGenomes: https://support.illumina.com/sequencing/sequencing_software/igenome.html <br>
You should specifically download the file Homo_sapiens_UCSC_hg38.tar.gz (Homo sapiens ; UCSC ; hg38). <br>
<br>
The hg38.gene_table.tsv table in Human_Genome_hg38 was taken from Louis Gauthier's GitHub HumanInteractome repository (in the annotations directory). This contains the gene information that was used for the initial library design. We use it to annotated our reads in the check_ontarget.py python script.
