#! /bin/bash

# Xavier Castellanos-Girouard
#
# Trimming Nterminal-tagged nested PCR reads using CutADAPT
#
# Date First Created: Nov 23rd 2025
# Date Last Modified: Nov 23rd 2025

## Note: Information in the fastq file sequence identifier follows the following format:
# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>


########## Initiate directories ##########

## Keep discarded reads here
mkdir Nterm/Discarded
mkdir Nterm/Discarded/plasmid_removal
mkdir Nterm/Discarded/DHFR_removal
mkdir Nterm/Discarded/linker_removal

## For files that are used as input downstream in the script
mkdir Nterm/Intermediate
mkdir Nterm/Intermediate/plasmid_removal
mkdir Nterm/Intermediate/DHFR_removal
mkdir Nterm/Intermediate/linker_removal
mkdir Nterm/Intermediate/polyG_removal




########## Filter reads that contain plasmid sequence ##########

echo "Trimming Nterminal flanking sequences..."

# Note: Some forward reads contain the plasmid sequence which harbored the
# 	donor DNA. We remove these reads as they are not informative ; we 
# 	are looking for donor DNAs which were inserted in the genome.
# 	Important: If the sequence appears in the forward read, we also 
# 	remove the corresponding reverse read and vice versa (-p argument).
#
# Note: Plasmid sequence: "TCCACGCGTAACTAAGTGGGC" (most common).
# 	Full sequence is likely something like: TCCACGCGTAACTAAGTGGGCCCGCCCCAACTGGGGTATCCTTT
# 	Partial match to miniplasmid precursor.
#
# Note: We filter these unwanted sequences by treating them as an "adapter" 
# 	that can occur anywhere within the read (-g argument of cutadapt). 
# 	In this case, untrimmed reads will contain our wanted reads, so we
# 	specify the output using the --untrimmed-output argument.


## First filter on forward primer (full length sequence)
#
# Note: We tolerate a 10% error rate (-e argument) and require an overlap
# 	of minumum 12 nucleotides (-O argument).
cutadapt -a TCCACGCGTAACTAAGTGGGCX `# Sequence we want to remove.` \
  -j 0 \
  -O 12 `# Minimum overlap` \
  -e 0.1 `# error rate` \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Nterm/Intermediate/plasmid_removal/noplasmidFWD_NLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads not containing plasmid sequence` \
  --untrimmed-paired-output Nterm/Intermediate/plasmid_removal/noplasmidFWD_NLib_NestedPCR_R2.fastq.gz `# Output for rev reads that pair with fwd reads not containing plasmid sequence` \
  ../data/DSP1374/HKCLLBGXL/fastq/Sample_N-Lib-Nested-PCR/N-Lib-Nested-PCR_S1_R1_001.fastq.gz `# raw fwd reads input file` \
  ../data/DSP1374/HKCLLBGXL/fastq/Sample_N-Lib-Nested-PCR/N-Lib-Nested-PCR_S1_R2_001.fastq.gz `# raw rev reads input file` \
  1> ./Nterm/Logs/Nterm_FWD_RemovePlasmid.txt



## Second filter on forward primer (partial length sequence)
#
# Note: Some reads contain the plasmid sequence only partially, 
# 	around the first dozen nucleotides. We can remove this with 
# 	a more specific filter that is strict enough to not remove
# 	other important sequences.

cutadapt -a TCCACGCGTAAG `# Sequence we want to remove.` \
  -j 0 \
  -e 1 `# error rate` \
  -O 12 `# Minimum overlap` \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Nterm/Intermediate/plasmid_removal/noplasmidFWD2_NLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads not containing partial plasmid sequence` \
  --untrimmed-paired-output Nterm/Intermediate/plasmid_removal/noplasmidFWD2_NLib_NestedPCR_R2.fastq.gz `# Output for rev reads that pair with fwd reads not containing partial plasmid sequence` \
  Nterm/Intermediate/plasmid_removal/noplasmidFWD_NLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  Nterm/Intermediate/plasmid_removal/noplasmidFWD_NLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Nterm/Logs/Nterm_FWD2_RemovePlasmid.txt


## To be thorough, also remove plasmid sequence if it appears in reverse primer
#
# Note: Reverse complement of plasmid sequence is the following: GCCCACTTAGTTACGCGTGGA
# Note: We tolerate a 10% error rate (-e argument) and require an overlap
# 	of minumum 20 nucleotides.

cutadapt -G GCCCACTTAGTTACGCGTGGA `# Sequence we want to remove.` \
  -j 0 \
  -e 1 `# error rate` \
  -O 20 `# Minimum overlap` \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Nterm/Intermediate/plasmid_removal/noplasmidREV_NLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads that pair with rev reads not containing plasmid sequence` \
  --untrimmed-paired-output Nterm/Intermediate/plasmid_removal/noplasmidREV_NLib_NestedPCR_R2.fastq.gz `# Output for rev reads not containing plasmid sequence` \
  ./Nterm/Intermediate/plasmid_removal/noplasmidFWD2_NLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  ./Nterm/Intermediate/plasmid_removal/noplasmidFWD2_NLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Nterm/Logs/Nterm_REV_RemovePlasmid.txt



## Also remove short plasmid sequence on reverse primer (TTACGCGTGGA)
cutadapt -G CTTACGCGTGGA `# Sequence we want to remove.` \
  -j 0 \
  -e 1 `# error rate` \
  -O 12 `# Minimum overlap` \
  --cores 5 \
  -o /dev/null `# Discard fwd reads with plasmid` \
  -p /dev/null `# Discard rev reads that pair with fwd reads containing plasmid sequence` \
  --untrimmed-output Nterm/Intermediate/plasmid_removal/noplasmidREV2_NLib_NestedPCR_R1.fastq.gz  `# Output for fwd reads that pair with rev reads not containing partial plasmid sequence` \
  --untrimmed-paired-output Nterm/Intermediate/plasmid_removal/noplasmidREV2_NLib_NestedPCR_R2.fastq.gz `# Output for rev reads not containing partial plasmid sequence` \
  Nterm/Intermediate/plasmid_removal/noplasmidREV_NLib_NestedPCR_R1.fastq.gz `# fwd reads input file` \
  Nterm/Intermediate/plasmid_removal/noplasmidREV_NLib_NestedPCR_R2.fastq.gz `# rev reads input file` \
  1> ./Nterm/Logs/Nterm_REV2_RemovePlasmid.txt



########## Trim DHFR sequence from reads ##########

## Trim DHFR from 5' of forward read

# Note: Anchored trimming is used (^), because DHFR is always at the 5' end.
# 	and we are looking for instances where the donor DNA (DHFR-Linker)
# 	is inserted without without indels in the donor. cutadapt -g ^ADAPTER
# Note: We tolerate a high error rate because practically all reads
#   contain the DHFR sequence (FWD primer sits on DHFR)
#
# Note: Upstream DHFR-F[1,2] sequence: "CAACCGGAATTGGGTACC"
#	The 5'-3' translation: QPELGT (Cterminal of DHFR-F[1,2])
cutadapt -g ^CAACCGGAATTGGGTACC `# Sequence we want to trim.` \
  -e 0.25 `# error rate` \
  -m 15 `# minimum sequence length after trim` \
  -j 0 \
  -o Nterm/Intermediate/DHFR_removal/noDHFR_FWD_NLib_NestedPCR_R1.fastq.gz \
  -p Nterm/Intermediate/DHFR_removal/noDHFR_FWD_NLib_NestedPCR_R2.fastq.gz \
  --discard-untrimmed \
  Nterm/Intermediate/plasmid_removal/noplasmidREV2_NLib_NestedPCR_R1.fastq.gz \
  Nterm/Intermediate/plasmid_removal/noplasmidREV2_NLib_NestedPCR_R2.fastq.gz \
  1> ./Nterm/Logs/Nterm_FWD_RemoveDHFR.txt



########## Trim linker sequence from reads ##########
# Note: Anchored trimming is used to remove linker. (cutadapt -g ^ADAPTER)
# Note: Linker information: "GGTGGCGGTGGCTCTGGAGGTGGTGGGTCCCCATCG" (length 36)
# 	Amino-acid translation: "GGGGSGGGGSPS" (length 12)
# Note: Because linker can lose an amino acid, we sequentially remove shorter 
# 	and shorter versions of the linker to recoup as many reads as possible.



LINKER_SEQ="GGTGGCGGTGGCTCTGGAGGTGGTGGGTCCCCATCG"
FULL_LEN=${#LINKER_SEQ}

# Input Files
IN_R1="Nterm/Intermediate/DHFR_removal/noDHFR_FWD_NLib_NestedPCR_R1.fastq.gz"
IN_R2="Nterm/Intermediate/DHFR_removal/noDHFR_FWD_NLib_NestedPCR_R2.fastq.gz"

# Final Output Files (We will append to these)
OUT_R1="Nterm/Intermediate/linker_removal/noLinker_NLib_NestedPCR_R1.fastq.gz"
OUT_R2="Nterm/Intermediate/linker_removal/noLinker_NLib_NestedPCR_R2.fastq.gz"

# Clear output files if they exist (to avoid appending to old runs)
> "$OUT_R1"
> "$OUT_R2"

# Temporary files for the "Next Round"
TEMP_UNTRIMMED_R1="temp_untrimmed_R1.fastq.gz"
TEMP_UNTRIMMED_R2="temp_untrimmed_R2.fastq.gz"

# Initialize: Copy input to the first "Untrimmed" bucket
cp "$IN_R1" "$TEMP_UNTRIMMED_R1"
cp "$IN_R2" "$TEMP_UNTRIMMED_R2"

echo "Starting Cascading Trim..."

# Loop from 36 down to 15
for (( i=$FULL_LEN; i>=15; i-- )); do
    
    # Define current adapter (Prefix)
    SUBSEQ="${LINKER_SEQ:0:$i}"
    
    echo "  Filtering for Length $i..."
    
    # 1. Run Cutadapt on the CURRENT batch of untrimmed reads
    #    - Matched reads -> Appended to FINAL Output
    #    - Untrimmed reads -> Saved to TEMP Output for next loop
    cutadapt \
        -g "X${SUBSEQ}" \
        -j 0 -e 0.15 \
        --untrimmed-output "next_round_R1.fastq.gz" \
        --untrimmed-paired-output "next_round_R2.fastq.gz" \
        -o "trimmed_chunk_R1.fastq.gz" \
        -p "trimmed_chunk_R2.fastq.gz" \
        "$TEMP_UNTRIMMED_R1" "$TEMP_UNTRIMMED_R2" \
        > /dev/null

    # 2. Append the successfully trimmed reads to our Final Pile
    cat "trimmed_chunk_R1.fastq.gz" >> "$OUT_R1"
    cat "trimmed_chunk_R2.fastq.gz" >> "$OUT_R2"

    # 3. Move the untrimmed reads to be the input for the next loop
    mv "next_round_R1.fastq.gz" "$TEMP_UNTRIMMED_R1"
    mv "next_round_R2.fastq.gz" "$TEMP_UNTRIMMED_R2"
done

# Cleanup
rm "trimmed_chunk_R1.fastq.gz" "trimmed_chunk_R2.fastq.gz" "$TEMP_UNTRIMMED_R1" "$TEMP_UNTRIMMED_R2"



########## Remove PolyG tails and pass standard quality filters ##########

cutadapt -G "G{55}"  `# Sequence we want to trim.` \
    --trim-n `# Trim uncalled nucleotides from read extremities` \
    -m 20 `# Minimum length of the read` \
    -o "Nterm/Intermediate/polyG_removal/noPolyG_NLib_NestedPCR_R1.fastq.gz" \
    -p "Nterm/Intermediate/polyG_removal/noPolyG_NLib_NestedPCR_R2.fastq.gz" \
    Nterm/Intermediate/linker_removal/noLinker_NLib_NestedPCR_R1.fastq.gz \
    Nterm/Intermediate/linker_removal/noLinker_NLib_NestedPCR_R2.fastq.gz \
    1> "./Nterm/Logs/Nterm_RemovePolyG.txt"


