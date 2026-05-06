# genomad
library(purrr)
library(readr)
library(dplyr)
library(stringr)

# metadata 
metadata <- read.table('data/metadata_genomes.tsv', sep = '\t', header = TRUE)


# list all summary files
f_viruses <- list.files(
  path = '012_genomad',
  pattern = "_virus_summary.tsv$",
  recursive = TRUE,
  full.names = TRUE)

# read and combine
viruses <- map_dfr(f_viruses, function(f) {
  genome_name <- str_extract(f, "H[0-9]+_[0-9]+")
  
  read_tsv(f, col_types = cols(.default = "c"),  # force all columns as character
    progress = FALSE,
    show_col_types = FALSE) %>%
    mutate(genome = genome_name) }) %>% 
  left_join(metadata, by = join_by('genome' == 'samples'))

# read plasmids into R 
f_plasmids <- list.files(
  path = '012_genomad',
  pattern = "_plasmid_summary.tsv$",
  recursive = TRUE,
  full.names = TRUE)

plasmids <- map_dfr(f_plasmids, function(f) {
  genome_name <- str_extract(f, "H[0-9]+_[0-9]+")
  
  read_tsv(f, col_types = cols(.default = "c"),  # force all columns as character
           progress = FALSE,
           show_col_types = FALSE) %>%
    mutate(genome = genome_name)
}) %>% 
  left_join(metadata, by = join_by('genome' == 'samples')) %>% 
  filter(plasmid_score > 0.9)

# viruses through time? 
viruses %>% 
  group_by(seq_name) %>% 
  reframe(n_genomes = n_distinct(genome))

