#! /usr/local/bin/R

library(tidyverse)

#################################################
# Join metadata with BGC info
#################################################

# command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# input files
genome_metadata_file <- args[1]
bgc_metadata_file <- args[2]
stb_tsv <- args[3]
output_tsv <- args[4]

# read in files
genome_metadata <- read_tsv(genome_metadata_file)
bgc_metadata <- read_tsv(bgc_metadata_file)
stb_tsv <- read_tsv(stb_tsv)

colnames(bgc_metadata) <- c("bgc_id", "scaffold_id", "description", "product", "bigscape_class", "organism", "taxonomy")

# join bgc metadata with STB genome info for corresponding scaffolds to genome ID
# then join with genome metadata
bgc_metadata <- left_join(bgc_metadata, stb_tsv, by = "scaffold_id")  %>% 
    select(mag_id, bgc_id, product, bgc_class)  %>% 
    left_join(genome_metadata)  %>% 
    select(mag_id, bgc_id, product, bgc_class, substrate, species, group)

write_tsv(bgc_metadata, output_tsv)