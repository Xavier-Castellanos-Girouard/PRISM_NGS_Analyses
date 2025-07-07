#! /bin/bash

# Xavier Castellanos-Girouard
#
# Trimming Cterminal-tagged nested PCR reads using CutADAPT
#
# Date First Created: June 25th 2025
# Date Last Modified: July 4th 2025

## Note: Information in the fastq file sequence identifier follows the following format:
# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>


########## Initiate directories ##########

## Keep discarded reads here
mkdir Cterm/Discarded
mkdir Cterm/Discarded/plasmid_removal
mkdir Cterm/Discarded/DHFR_removal
mkdir Cterm/Discarded/linker_removal

## For files that are used as input for downstream in the script
mkdir Cterm/Intermediate
mkdir Cterm/Intermediate/plasmid_removal
mkdir Cterm/Intermediate/DHFR_removal
mkdir Cterm/Intermediate/linker_removal
mkdir Cterm/Intermediate/polyG_removal




########## Filter reads that contain plasmid sequence ##########

echo "Trimming Cterminal flanking sequences..."

# Note: Some forward reads contain the plasmid sequence which harbored the
# 	donor DNA. We remove these reads as they are not informative ; we 
# 	are looking for donor DNAs which were inserted in the genome.
# 	Important: If the sequence appears in the forward read, we also 
# 	remove the corresponding reverse read and vice versa (-p argument).
#
# Note: Plasmid sequence: "CCATCGCGA" and "CCATCGCGAGTTTTAGAGCTAGAAATAGCAAG". 
# 	First sequence is a match to miniplasmid precursor (backbone).
#	Second sequence is a partial match to 5'flanking regions for sgRNAs
#
# Note: We filter these unwanted sequences by treating them as an "adapter" 
# 	that occurs at the end (flexible) of the read (-a argument of cutadapt). 
# 	In this case, untrimmed reads will contain our wanted reads, so we
# 	specify the output using the --untrimmed-output argument.


### First filter on forward primer (full length sequence)
#
# Note: We tolerate a 25% error rate (-e argument) and require an overlap
# 	of minimum 18 nucleotides (-O argument).
cutadapt -a CCATCGCGAGTTTTAGAGCTAGAAATAGCAAGX `# Sequence we want to remove.` \
  -O 24 `# Minimum overlap` \
  -e 0.1 `# error rate` \
  --cores 5 \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Cterm/Intermediate/plasmid_removal/noplasmidFWD_CLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads not containing plasmid sequence` \
  --untrimmed-paired-output Cterm/Intermediate/plasmid_removal/noplasmidFWD_CLib_NestedPCR_R2.fastq.gz `# Output for rev reads that pair with fwd reads not containing plasmid sequence` \
  ../data/DSP1472/HLTJ3BGXM/fastq/Sample_Human_Clib/Human_Clib_S1_R1_001.fastq.gz `# raw fwd reads input file` \
  ../data/DSP1472/HLTJ3BGXM/fastq/Sample_Human_Clib/Human_Clib_S1_R2_001.fastq.gz `# raw rev reads input file` \
  1> ./Cterm/Logs/Cterm_FWD_RemovePlasmid.txt

# If Discarded reads are of interest, use the following:
#  -o Cterm/Discarded/plasmid_removal/plasmidFWD_CLib_NestedPCR_R1.fastq.gz \
#  -p Cterm/Discarded/plasmid_removal/plasmidFWD_CLib_NestedPCR_R2.fastq.gz \



### Second filter on forward primer (partial length sequence)
#
# Note: Some reads contain the plasmid sequence only partially, 
# 	around the first dozen nucleotides. We can remove this with 
# 	a more specific filter that is strict enough to not remove
# 	other important sequences.

cutadapt -a CCATCGCGAGTTTTA `# Sequence we want to remove.` \
  -e 1 `# error rate` \
  -O 12 `# Minimum overlap` \
  --cores 5 \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Cterm/Intermediate/plasmid_removal/noplasmidFWD2_CLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads not containing partial plasmid sequence` \
  --untrimmed-paired-output Cterm/Intermediate/plasmid_removal/noplasmidFWD2_CLib_NestedPCR_R2.fastq.gz `# Output for rev reads that pair with fwd reads not containing partial plasmid sequence` \
  Cterm/Intermediate/plasmid_removal/noplasmidFWD_CLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  Cterm/Intermediate/plasmid_removal/noplasmidFWD_CLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Cterm/Logs/Cterm_FWD2_RemovePlasmid.txt

# If Discarded reads are of interest, use the following:
#  -o Cterm/Discarded/plasmid_removal/plasmidFWD2_CLib_NestedPCR_R1.fastq.gz \
#  -p Cterm/Discarded/plasmid_removal/plasmidFWD2_CLib_NestedPCR_R2.fastq.gz \



### Third filter on forward primer (alternative plasmid sequence?)
cutadapt -a AGACGTTCTGCTTCACTCTCCCCATCT `# Sequence we want to remove.` \
  -e 0.1 `# error rate` \
  -O 24 `# Minimum overlap` \
  --cores 5 \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Cterm/Intermediate/plasmid_removal/noplasmidFWD3_CLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads not containing partial plasmid sequence` \
  --untrimmed-paired-output Cterm/Intermediate/plasmid_removal/noplasmidFWD3_CLib_NestedPCR_R2.fastq.gz `# Output for rev reads that pair with fwd reads not containing partial plasmid sequence` \
  Cterm/Intermediate/plasmid_removal/noplasmidFWD2_CLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  Cterm/Intermediate/plasmid_removal/noplasmidFWD2_CLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Cterm/Logs/Cterm_FWD3_RemovePlasmid.txt

# If Discarded reads are of interest, use the following:
#  -o Cterm/Discarded/plasmid_removal/plasmidFWD3_CLib_NestedPCR_R1.fastq.gz \
#  -p Cterm/Discarded/plasmid_removal/plasmidFWD3_CLib_NestedPCR_R2.fastq.gz \



### To be thorough, also remove plasmid sequence if it appears in reverse primer
#
# Note: Reverse complement of plasmid sequence is the following: CTTGCTATTTCTAGCTCTAAAACTCGCGATGG
# Note: We tolerate a 25% error rate (-e argument) and require an overlap
# 	of minimum 18 nucleotides.

cutadapt -G CTTGCTATTTCTAGCTCTAAAACTCGCGATGG `# Sequence we want to remove.` \
  -e 0.1 `# error rate` \
  -O 24 `# Minimum overlap` \
  --cores 5 \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Cterm/Intermediate/plasmid_removal/noplasmidREV_CLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads that pair with rev reads not containing plasmid sequence` \
  --untrimmed-paired-output Cterm/Intermediate/plasmid_removal/noplasmidREV_CLib_NestedPCR_R2.fastq.gz `# Output for rev reads not containing plasmid sequence` \
  ./Cterm/Intermediate/plasmid_removal/noplasmidFWD3_CLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  ./Cterm/Intermediate/plasmid_removal/noplasmidFWD3_CLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Cterm/Logs/Cterm_REV_RemovePlasmid.txt

# If Discarded reads are of interest, use the following:
#  -o Cterm/Discarded/plasmid_removal/plasmidREV_CLib_NestedPCR_R1.fastq.gz \
#  -p Cterm/Discarded/plasmid_removal/plasmidREV_CLib_NestedPCR_R2.fastq.gz \



### Also remove short plasmid sequence on reverse primer (ACTCGCGATGG)
cutadapt -G TAAAACTCGCGATGG `# Sequence we want to remove.` \
  -e 0.1 `# error rate` \
  -O 15 `# Minimum overlap` \
  --cores 5 \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Cterm/Intermediate/plasmid_removal/noplasmidREV2_CLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads that pair with rev reads not containing partial plasmid sequence` \
  --untrimmed-paired-output Cterm/Intermediate/plasmid_removal/noplasmidREV2_CLib_NestedPCR_R2.fastq.gz `# Output for rev reads not containing partial plasmid sequence` \
  Cterm/Intermediate/plasmid_removal/noplasmidREV_CLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  Cterm/Intermediate/plasmid_removal/noplasmidREV_CLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Cterm/Logs/Cterm_REV2_RemovePlasmid.txt

# If Discarded reads are of interest, use the following:
#  -o Cterm/Discarded/plasmid_removal/plasmidREV2_CLib_NestedPCR_R1.fastq.gz \
#  -p Cterm/Discarded/plasmid_removal/plasmidREV2_CLib_NestedPCR_R2.fastq.gz \



########## Trim DHFR sequence from reads ##########

## Trim DHFR from 5' flanking region of forward read

# *UPDATE THIS NOTE*Note: Anchored trimming is used (^), because DHFR is always at the 5' end.
# 	and we are looking for instances where the donor DNA (DHFR-Linker)
# 	is inserted without without indels in the donor. cutadapt -g ^ADAPTER
# *UPDATE THIS NOTE*Note: We tolerate a high error rate because practically all reads
#   contain the DHFR sequence (FWD primer sits on DHFR)
#
# Note: Upstream reverse complement DHFR-F[1,2] sequence: "gttcaatggtcgaaccat"
#	The 3'-5' translation: MVRPLN (Nterminal of DHFR-F[1,2])
cutadapt -g ^GTTCAATGGTCGAACCAT `# Sequence we want to trim.` \
  -e 0.25 `# error rate` \
  -o Cterm/Intermediate/DHFR_removal/noDHFR_FWD_CLib_NestedPCR_R1.fastq.gz \
  -p Cterm/Intermediate/DHFR_removal/noDHFR_FWD_CLib_NestedPCR_R2.fastq.gz \
  --untrimmed-output Cterm/Discarded/DHFR_removal/DHFR_FWD_CLib_NestedPCR_R1.fastq.gz \
  --untrimmed-paired-output Cterm/Discarded/DHFR_removal/DHFR_FWD_CLib_NestedPCR_R2.fastq.gz \
  Cterm/Intermediate/plasmid_removal/noplasmidREV2_CLib_NestedPCR_R1.fastq.gz \
  Cterm/Intermediate/plasmid_removal/noplasmidREV2_CLib_NestedPCR_R2.fastq.gz \
  1> ./Cterm/Logs/Cterm_FWD_RemoveDHFR.txt


########## Trim linker sequence from reads ##########
# Note: Anchored trimming is used to remove linker. (cutadapt -g ^ADAPTER)
# Note: Linker information: "AGACCCACCGCCTCCTGATCCGCCACCGCC" (length 30)
# 	Amino-acid translation (3'-5'): "GGGGSGGGGS" (length 10)
# Note: Because linker can lose an amino acid, we sequentially remove shorter 
# 	and shorter versions of the linker to recoup as many reads as possible.
# Note: We tolerate only one error, otherwise you may have instances where 
# 	we trim the stop codon.


##### Sequential trimming #####

i=1
curr_input_R1="./Cterm/Intermediate/DHFR_removal/noDHFR_FWD_CLib_NestedPCR_R1.fastq.gz"
curr_input_R2="./Cterm/Intermediate/DHFR_removal/noDHFR_FWD_CLib_NestedPCR_R2.fastq.gz"
linker="AGACCCACCGCCTCCTGATCCGCCACCGCC"

intermediates_dir="./Cterm/Intermediate/linker_removal"
discarding_dir="./Cterm/Discarded/linker_removal"
  
#echo $intermediates_dir

while (( ${#linker} >= 15 )); do
  
  echo $linker
  
  cutadapt -g "^$linker"  `# Sequence we want to trim.` \
    -e 1 `# number of tolerated errors` \
    -o "$intermediates_dir/noLinker${i}_CLib_NestedPCR_R1.fastq.gz" \
    -p "$intermediates_dir/noLinker${i}_CLib_NestedPCR_R2.fastq.gz" \
    --untrimmed-output "$discarding_dir/noLinker${i}_CLib_NestedPCR_R1.fastq.gz" \
    --untrimmed-paired-output "$discarding_dir/noLinker${i}_CLib_NestedPCR_R2.fastq.gz" \
    "$curr_input_R1" \
    "$curr_input_R2" \
    1> "./Cterm/Logs/Cterm_RemoveLinker${i}.txt"
  
  curr_input_R1="$discarding_dir/noLinker${i}_CLib_NestedPCR_R1.fastq.gz"
  curr_input_R2="$discarding_dir/noLinker${i}_CLib_NestedPCR_R2.fastq.gz"
  ((i++))
  
  # Trim 3 nt from linker until 15 nt (if statement avoids preemptive break)
  if (( ${#linker} > 15 )); then
    linker="${linker:0:${#linker}-3}"
  else
    break
  fi

done


##### Short fragment trimming #####
# Note: Certain 5' sequences will not be aligned perfectly.
# 	Use a stringent filter with 3' end of the linker to
# 	accomodate for this. 

cutadapt -g "TCCGCCACCGCC"  `# Sequence we want to trim.` \
    -O 11 \
    -e 1 `# number of tolerated errors` \
    -o "$intermediates_dir/noLinker7_CLib_NestedPCR_R1.fastq.gz" \
    -p "$intermediates_dir/noLinker7_CLib_NestedPCR_R2.fastq.gz" \
    --untrimmed-output "$discarding_dir/noLinker7_CLib_NestedPCR_R1.fastq.gz" \
    --untrimmed-paired-output "$discarding_dir/noLinker7_CLib_NestedPCR_R2.fastq.gz" \
    "./Cterm/Discarded/linker_removal/noLinker6_CLib_NestedPCR_R1.fastq.gz" \
    "./Cterm/Discarded/linker_removal/noLinker6_CLib_NestedPCR_R2.fastq.gz" \
    1> "./Cterm/Logs/Cterm_RemoveLinker7.txt"


### Combine different trimmings into a single fastq.gz file
cat Cterm/Intermediate/linker_removal/*_R1.fastq.gz > Cterm/Intermediate/linker_removal/noLinker_CLib_NestedPCR_R1.fastq.gz
cat Cterm/Intermediate/linker_removal/*_R2.fastq.gz > Cterm/Intermediate/linker_removal/noLinker_CLib_NestedPCR_R2.fastq.gz




########## Remove PolyG tails ##########

cutadapt -G "G{75}" `# Sequence we want to trim.` \
    -q 10 `# Phred score cutoff` \
    --trim-n `# Trim uncalled nucleotides from read extremities` \
    --minimum-length 10 `# Minimum length of the read` \
    -o "Cterm/Intermediate/polyG_removal/noPolyG_CLib_NestedPCR_R1.fastq.gz" \
    -p "Cterm/Intermediate/polyG_removal/noPolyG_CLib_NestedPCR_R2.fastq.gz" \
    "$intermediates_dir/noLinker_CLib_NestedPCR_R1.fastq.gz" \
    "$intermediates_dir/noLinker_CLib_NestedPCR_R2.fastq.gz" \
    1> "./Cterm/Logs/Cterm_RemovePolyG.txt"



