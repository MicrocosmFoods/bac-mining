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
output_annotation_tsv <- args[4]
output_substrate_counts_tsv <- args[5]
output_phylo_group_counts_tsv <- args[6]

# read in files
genome_metadata <- read_tsv(genome_metadata_file)
bgc_info <- read_tsv(bgc_metadata_file)
stb_tsv <- read_tsv(stb_tsv)

colnames(bgc_info) <- c("bgc_id", "scaffold_id", "description", "product", "bgc_class", "organism", "taxonomy")
colnames(stb_tsv) <- c("mag_id", "scaffold_id")
# join bgc metadata with STB genome info for corresponding scaffolds to genome ID
# then join with genome metadata
bgc_metadata <- left_join(bgc_info, stb_tsv, by = "scaffold_id")  %>% 
    filter(!str_starts(bgc_id, "BGC"))  %>% 
    select(mag_id, bgc_id, product, bgc_class)  %>% 
    left_join(genome_metadata)  %>% 
    select(mag_id, bgc_id, product, bgc_class, substrate, species, group)

bgc_substrate_type_counts <- bgc_metadata  %>% 
    group_by(substrate, bgc_class)  %>% 
    count()  %>% 
    arrange(desc(n))

bgc_phylo_group_type_counts <- bgc_metadata  %>% 
    group_by(group, bgc_class)  %>% 
    count()  %>% 
    arrange(desc(n))

write_tsv(bgc_metadata, output_annotation_tsv)
write_tsv(bgc_substrate_type_counts, output_substrate_counts_tsv)
write_tsv(bgc_phylo_group_type_counts, output_phylo_group_counts_tsv)