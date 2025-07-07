# Xavier Castellanos-Girouard
#
# Final Analysis and Figures for on-Targeting and insert shifts
#
# Date First Created: May 22nd 2025
# Date Last Modified: July 5th 2025

#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

#### Import data ####

#TEST <- fread("./Nterm/results/Nterm_mapped_clean.csv")

Nterm_annotated_DF <- fread("./Nterm/results/Nterm_annotated_onTarget.csv")
Nterm_annotated_DF <- Nterm_annotated_DF[Nterm_annotated_DF$mapping_quality >= 20,]
Nterm_annotated_DF <- Nterm_annotated_DF %>% select(-c("V1"))
Nterm_annotated_DF <- Nterm_annotated_DF[(Nterm_annotated_DF$category == "proper_pairs") | (Nterm_annotated_DF$category == "reverse_only")]

Cterm_FullLinker_DF <- fread("./Cterm/results/Cterm_annotated_onTarget.csv")
Cterm_FullLinker_DF <- Cterm_FullLinker_DF[Cterm_FullLinker_DF$mapping_quality >= 20,]
Cterm_FullLinker_DF <- Cterm_FullLinker_DF[(Cterm_FullLinker_DF$category == "proper_pairs") | (Cterm_FullLinker_DF$category == "reverse_only")]

#Nterm_FullLinker_DF <- read.csv("./Nterm/results/Nterm_FrameShift_R1.csv")
#Cterm_FullLinker_DF <- read.csv("./Cterm/results/Cterm_FrameShift_R1.csv")

# Masih results
#R1 <- read.csv("/home/xavier/Desktop/Mammalian_CRISPR/PRISM_Nested_PCR/res_n_library/plasmidLinkerPolyGQTrimmed.ReadFiltered.sorted.R1.genome_annotation.tsv", sep = "\t", header = FALSE)
#R1 <- R1[!(R1$V1 %in% c("__alignment_not_unique", "__ambiguous", "__not_aligned", "__too_low_aQual")),]

#R2 <- read.csv("/home/xavier/Desktop/Mammalian_CRISPR/PRISM_Nested_PCR/res_n_library/plasmidLinkerPolyGQTrimmed.ReadFiltered.sorted.R2.genome_annotation.tsv", sep = "\t", header = FALSE)
#R2 <- R2[!(R2$V1 %in% c("__alignment_not_unique", "__ambiguous", "__not_aligned", "__no_feature", "__too_low_aQual")),]


#### Format Data ####

## Make a user-defined function to extract gene names and accession IDs
get_unique_gene_ID <- function(x){ # Sometimes on of the reads does not have an assigned gene/ID
  if(length(unique(x[which(x!="None")]))==0){return("None")} # If no read has a gene, return None
  else{return(unique(x[which(x!="None")]))} # If at least one of the reads has a gene, return gene name 
}

## Format proper pair reads, such that only one entry is retained
Nterm_annotated_proper_pair_DF <- 
  Nterm_annotated_DF %>%
  filter(category == "proper_pairs") %>%
  group_by(query_name) %>%
  summarise(query_name = unique(query_name),
            chromosome = unique(chromosome),
            #start = start[which(is_read1 == TRUE)], # First nt of read1 is always start ***NOT TRUE, FIX THIS
            #end = end[which(is_read2 == TRUE)], # Last nt of read2 is always the end
            mapping_quality = unique(mapping_quality),
            is_read1 = TRUE,
            is_read2 = TRUE,
            flag = 0,
            category = unique(category),
            gene_name = get_unique_gene_ID(gene_name),
            accession_ID = get_unique_gene_ID(accession_ID),
            onTarget = any(onTarget == TRUE)) %>%
  ungroup()

Nterm_annotated_reverse_only_DF <-
  Nterm_annotated_DF %>%
  filter(category == "reverse_only")

#### Check on-Target ratio ####

Nterm_FullLink_onTar_DF <-
  Nterm_annotated_proper_pair_DF  %>%
  dplyr::group_by(onTarget) %>%
  summarize(read_count = n())

Nterm_FullLink_onTar_DF <-
  Nterm_annotated_DF  %>%
  filter(category == "proper_pairs") %>%
  dplyr::group_by(onTarget) %>%
  summarize(read_count = n())

Nterm_FullLink_onTar_DF <-
  Nterm_annotated_DF %>%
  dplyr::group_by(onTarget) %>%
  summarize(read_count = n())

Nterm_FullLink_onTar_DF <-
  Nterm_annotated_reverse_only_DF %>%
  dplyr::group_by(onTarget) %>%
  summarize(read_count = n())

Cterm_FullLink_onTar_DF <-
  Cterm_FullLinker_DF %>%
  dplyr::group_by(onTarget) %>%
  summarize(read_count = n())


#### Make Figures for Nterm ####


Nterm_FullLink_onTar_p <-
  ggplot(data = Nterm_FullLink_onTar_DF, 
         mapping = aes(x = onTarget, 
                       y = read_count/sum(read_count)*100)) +
  geom_col(fill = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black"),
    axis.line.y = element_line(color ="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black"),
    axis.text.y.left = element_text(color = "black"),
    axis.title.x.bottom = element_text(color = "black"),
    axis.title.y.left = element_text(color = "black")) +
  scale_y_continuous(limits = c(0,100),
                     breaks = seq(0,100, 10)) +
  xlab("Insert On Target") +
  ylab("Read percentage")

Nterm_FullLink_onTar_p

##
Nterm_FullLink_onTar_violin_p <-
  ggplot(data = Nterm_FullLinker_DF, 
         mapping = aes(x = onTarget, 
                       y = counts)) +
  geom_violin(fill = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black"),
    axis.line.y = element_line(color ="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black"),
    axis.text.y.left = element_text(color = "black"),
    axis.title.x.bottom = element_text(color = "black"),
    axis.title.y.left = element_text(color = "black")) +
  xlab("Insert On Target") +
  ylab("counts") +
  scale_y_log10()

Nterm_FullLink_onTar_violin_p



#### Make figures for Cterm ####


Cterm_FullLink_onTar_p <-
  ggplot(data = Cterm_FullLink_onTar_DF, 
         mapping = aes(x = onTarget, 
                       y = read_count/sum(read_count)*100)) +
  geom_col(fill = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black"),
    axis.line.y = element_line(color ="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black"),
    axis.text.y.left = element_text(color = "black"),
    axis.title.x.bottom = element_text(color = "black"),
    axis.title.y.left = element_text(color = "black")) +
  scale_y_continuous(limits = c(0,100),
                     breaks = seq(0,100, 10)) +
  xlab("Insert On Target") +
  ylab("Read percentage")

Cterm_FullLink_onTar_p

##
Cterm_FullLink_onTar_violin_p <-
  ggplot(data = Cterm_FullLinker_DF, 
         mapping = aes(x = onTarget, 
                       y = counts)) +
  geom_violin(fill = "black") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black"),
    axis.line.y = element_line(color ="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black"),
    axis.text.y.left = element_text(color = "black"),
    axis.title.x.bottom = element_text(color = "black"),
    axis.title.y.left = element_text(color = "black")) +
  xlab("Insert On Target") +
  ylab("counts") +
  scale_y_log10()

Cterm_FullLink_onTar_violin_p