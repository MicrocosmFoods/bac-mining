#! /usr/local/bin/R

library(tidyverse)

######################################################################################
# Merge peptide stats with metadata
# Merge results from DeepSig, peptides.py, and DIAMOND Blastp results against database
######################################################################################

# command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# input files
peptides_info_tsv <- args[1]
deep_sig_tsv <- args[2]
diamond_tsv <- args[3]
genome_metadata <- args[4]
output_tsv <- args[5]

# read in files
peptides_info <- read_tsv(peptides_info_tsv)
deepsig_info <- read_tsv(deep_sig_tsv, col_names = c("peptide_id", 
                                                    "tool",
                                                    "deepsig_feature",
                                                    "deepsig_feature_start",
                                                    "deepsig_feature_end",
                                                    "deepsig_feature_score",
                                                    "tmp1",
                                                    "tmp2",
                                                    "deepsig_description"))
diamond_blast_results <- read_tsv(diamond_tsv)  %>% 
    mutate(peptide_id = qseqid)  %>% 
    select(-qseqid)
genome_metadata <- read_tsv(genome_metadata)  %>% 
    select(mag_id, substrate, species, group)

# merge peptides info with deepsig results
# merge with diamond blastp results
# create mag_id column for merging with metadata
# merge with metadata
all_peptide_info_metadata <- left_join(deepsig_info, peptides_info)  %>% 
    left_join(diamond_blast_results)  %>% 
    mutate(mag_id = str_extract(peptide_id, "^.*?(?=_id_)"))  %>% 
    left_join(genome_metadata)  %>% 
    select(mag_id, peptide_id, everything())

# write to tsv
write_tsv(all_peptide_info_metadata, output_tsv)

