# genomad
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 12),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 12)))

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


# Virus clusters based on MMseq2 (identity 90%, coverage 80%)
virus_clusters <- read_tsv('012_genomad/virus_mmseq2/virus_clusters.tsv', col_names = c('cluster', 'member')) %>%
  mutate(genome1 = str_extract(member, "H[0-9]+_[0-9]+"), 
         genome2 = str_extract(cluster, "H[0-9]+_[0-9]+"), 
         seq_name = str_remove(cluster, "^[^|]*\\|")) %>% 
  filter(!genome1 %in% c("H005_126", "H005_94"), 
         !genome2 %in% c("H005_126", "H005_94")) %>% 
  left_join(metadata, by = join_by('genome1' == 'samples')) %>% 
  left_join(select(viruses, seq_name, topology, virus_score), by = 'seq_name') %>% 
  filter(virus_score > 0.9)
  

# Count shared clusters and who has them?
virus_clusters %>% 
  group_by(cluster) %>% 
  reframe(n_genomes = n_distinct(genome1), 
          name_genomes = list(genome1))

virus_clusters %>% 
  ggplot(aes(x = seq_name, fill = as.factor(STRAINS))) +
  geom_bar() +
  coord_flip() +
  facet_wrap(~topology, scales = 'free', nrow = 3) +
  labs(x = 'Number of clusters detected', y = 'Viral cluster', fill = 'Strain')
ggsave('out/genomad/number_clusters_genome.png')

# How do viral clusters move across strains over time?
virus_clusters %>%
  ggplot(aes(x = as.factor(TIMEPOINT), y = seq_name, fill = as.factor(STRAINS))) +
  geom_tile() + 
  facet_wrap(~topology, scales = 'free', nrow = 3) +
  labs(x = 'Time point', y = 'Viral cluster', fill = 'Strain')
ggsave('out/genomad/viral_heatmap.png')


# Plasmids 
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
    mutate(genome = genome_name) }) %>% 
  left_join(metadata, by = join_by('genome' == 'samples')) %>% 
  filter(plasmid_score > 0.9)

plasmid_clusters <- read_tsv('012_genomad/plasmid_mmseq2/plasmid_clusters.tsv', col_names = c('cluster', 'member')) %>%
  mutate(genome1 = str_extract(member, "H[0-9]+_[0-9]+"), 
         genome2 = str_extract(cluster, "H[0-9]+_[0-9]+"), 
         seq_name = str_remove(cluster, "^[^|]*\\|")) %>% 
  filter(!genome1 %in% c("H005_126", "H005_94"), 
         !genome2 %in% c("H005_126", "H005_94")) %>% 
  left_join(metadata, by = join_by('genome1' == 'samples')) %>% 
  left_join(select(plasmids, seq_name, topology, plasmid_score), by = 'seq_name') %>% 
  filter(plasmid_score > 0.9)

plasmid_clusters %>% 
  ggplot(aes(x = as.factor(TIMEPOINT), y = seq_name, fill = as.factor(STRAINS))) +
  geom_tile() + 
  facet_wrap(~topology, scales = 'free', nrow = 3) +
  labs(x = 'Time point', y = 'Plasmid cluster', fill = 'Strain')
ggsave('out/genomad/plasmid_heatmap.png')
