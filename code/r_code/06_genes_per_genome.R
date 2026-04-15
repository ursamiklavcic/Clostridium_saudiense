library(dplyr)
library(tidyr)
library(ggplot2)

# gene presence absence 

gene_pa <- read.table('004_prokka/roary/gene_presence_absence.csv', sep = ',', header = T) %>% 
  pivot_longer(names_to = 'genome', values_to = 'PA', cols = starts_with('H00')) %>%  
  mutate(PA_01 = ifelse(PA != '', 1, 0), 
         genome = str_remove_all(genome, ".fasta"), 
         genome = str_remove_all(genome, "_hybrid"))

metadata_genomes <- read.table('~/projects/C_saudiense/metadata_genomes.tsv', sep = '\t', header = T) %>%  
  rename('genome' = 'samples')
genes <- left_join(gene_pa, metadata_genomes, by = 'genome')

# Number of genes per genome 
n_genes <- genes %>%  filter(PA_01 == 1) %>%  
  group_by(genome) %>%  
  reframe(n_genes = n_distinct(Gene))

# Number of genes in core set for each strain? (99% genomes)
core_genes <- genes %>% 
  filter(No..isolates >= 0.99 * 42) %>%  
  distinct(Gene)

n_core <- genes %>%  
  filter(PA_01 == 1, Gene %in% core_genes$Gene) %>% 
  group_by(genome, STRAIN_GROUP, STRAINS) %>%  
  reframe(n_core_genes = n_distinct(Gene))

# Number of accessory genes per genome
n_accessory <- genes %>%  
  filter(PA_01 == 1, !Gene %in% core_genes$Gene) %>% 
  group_by(genome, STRAIN_GROUP, STRAINS) %>%  
  reframe(n_accessory_genes = n_distinct(Gene))

# Genes present in only one genome 
unique_genes <- genes %>%  
  filter(PA_01 == 1) %>%  
  group_by(Gene) %>% 
  reframe(n_genomes = n_distinct(genome)) %>%  
  filter(n_genomes <= 1)

unique_genes_genome <- genes %>% 
  filter(PA_01 == 1, Gene %in% unique_genes$Gene) 

# How many unique genes per genome? 
n_unique_genes <- genes %>% 
  filter(PA_01 == 1, Gene %in% unique_genes$Gene) %>% 
  group_by(genome) %>%  
  reframe(n_unique_genes = n_distinct(Gene))

# One table 
all <- n_genes %>%  
  full_join(n_core, by = 'genome') %>% 
  full_join(n_accessory, by = c('genome', 'STRAIN_GROUP', 'STRAINS')) %>%  
  full_join(n_unique_genes, by = 'genome')

write.csv(all, 'out/genes_genome.csv')
