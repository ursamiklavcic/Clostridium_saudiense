
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)

roary <- read.table('011_analysis_with_NCBI_genomes/data/gene_presence_absence.csv', sep = ',', header = TRUE) %>%  
  pivot_longer(names_to = 'genome', values_to = 'PA', cols = starts_with(c('H00', 'GC'))) %>%  
  mutate(PA_01 = ifelse(PA != '', 1, 0)) %>% 
  mutate(group = ifelse(grepl('H00', genome), 'Your_42', 'Public_37')) 

roary %>% 
  filter(No..isolates > 70) %>% 
  ggplot(aes(x = Gene, y = genome, fill = as.factor(PA_01))) +
  geom_tile()


roary %>% 
  ggplot(aes (y = genome, x = n_distinct(PA), fill = as.factor(PA_01))) +
  geom_col() +
  labs(x = 'Number of genes', y = '', fill = '1 = present\n0=absent')
ggsave('out/roary/all_genomes_present_absent_genes.png', dpi =400)

# Number of present genes per genome 
roary %>% 
  group_by(genome, group) %>% 
  reframe(present = sum(PA_01 == 1)) %>% 
  ggplot(aes(y = genome, x = present, fill = group)) +
  geom_col() +
  labs(x = 'Number of genes present per genome')
ggsave('out/roary/no_present_genes_genome.png', dpi= 400)

# Number of absent genes per genome 
roary %>% 
  group_by(genome, group) %>% 
  reframe(absent = sum(PA_01 == 0)) %>% 
  ggplot(aes(y = genome, x = absent, fill = group)) +
  geom_col() +
  labs(x = 'Number of genes absent per genome')
ggsave('out/roary/no_absent_genes_genome.png', dpi= 400)

# calculate core genome
length(unique(roary$genome))
#79

roary %>% 
  filter(PA_01 == 1) %>% 
  group_by(Gene) %>% 
  reframe(n_genomes = n_distinct(genome)) # %>% 
  #filter(n_genomes > as.integer(79*0.99)) 73
  #filter(n_genomes > as.integer(79*0.95)) 747
  #filter(n_genomes > as.integer(79*0.90)) 1565

# Geni prisotni v več kot 90% genomov, kateri bodo štrleli ven? 
roary %>% filter(No..isolates > as.integer(79*0.9)) %>% 
  ggplot(aes(x = Gene, y = genome, fill = as.factor(PA_01))) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, size = 2)) +
  labs(y = '', x = 'Genes', fill = 'Ali je gen prisoten?\n1=DA, 0=NE')


# Core genome by groups? 
roary %>% 
  filter(No..isolates > 10) %>% 
  group_by(Gene) %>% 
  reframe(my_genomes = sum(PA_01[grepl('H00', genome)]), 
          NCBI = sum(PA_01[grepl('GC', genome)])) %>%  
  pivot_longer(c(my_genomes, NCBI), names_to = 'group', values_to = 'presence') %>% 
  ggplot(aes(x = presence, fill = group)) +
  geom_histogram(bins = 20, position = 'identity') +
  facet_wrap(~group) +
  labs(x = 'Number of genes present in my genomes VS NCBI genomes', y = 'Number of genomes with x genes', fill = '')
ggsave('out/roary/genomes_genes_histogram.png', dpi=400)  
  
roary %>% 
  filter(No..isolates > 10) %>% 
  group_by(Gene) %>% 
  reframe(my_genomes = sum(PA_01[grepl('H00', genome)])/42 *100, 
          NCBI = sum(PA_01[grepl('GC', genome)])/37 * 100) %>%  
  pivot_longer(c(my_genomes, NCBI), names_to = 'group', values_to = 'presence') %>%  
  ggplot(aes(y = group, x = presence, fill = group)) +
  #geom_point(size = 3, position = position_dodge(width = 0.8))
  geom_violin()
ggsave('out/roary/violin_presence_genes.png', dpi=400)


##
# Correlation between the abundance of a strain and accessory genome 
metadata <- read.table('data/metadata_genomes.tsv', header = TRUE)
genes_metagenomes <- readRDS('006_quantify_genes_Lara/genes_norm_count.RDS') 

n_accessory_genes <- roary %>% filter(genome %in% metadata$samples) %>% 
  filter(No..isolates < 41, PA_01 == 1) %>% 
  group_by(genome) %>% 
  reframe(n_accessory_genes = n_distinct(Gene)) %>% 
  left_join(metadata, by = join_by('genome' == 'samples')) %>% 
  group_by(STRAINS) %>% 
  mutate(mean_accessory_genes = as.integer(mean(n_accessory_genes)))

n_accessory_genes %>% 
  group_by(STRAINS) %>% 
  mutate(genomes_per_strain = n_distinct(genome)) %>% 
  ggplot(aes(x = as.factor(STRAINS), y = n_accessory_genes, fill= as.factor(genomes_per_strain))) +
  geom_col() +
  labs(x = 'Strains', y = 'Number of accessory genomes', fill = 'Genomes\nper strain')
ggsave('out/roary/n_accesory_genes_strain.png', dpi=400)

# Correlation strain abundance and size accessory genome
genes_metagenomes %>% 
  filter(person == 'H', biota == 'untreated sample') %>% 
  left_join(n_accessory_genes, by = join_by('STRAIN' == 'STRAINS'), relationship = "many-to-many") %>% 
  select(STRAIN, mean_accessory_genes, sum_norm_count) %>%  
  distinct() %>% 
  ggplot(aes(x = mean_accessory_genes, y = sum_norm_count)) +
  geom_point(aes(color = as.factor(STRAIN))) +
  geom_smooth(method = 'lm') +
  scale_y_log10() +
  labs(x = 'Number of accessory genes', y = 'Normlaized read counts', color = 'STRAIN')
ggsave('out/roary/no_accesory_genes_abundance_strain.png', dpi=400)

# remove low coverage genomes: 
# GCA_036305605.1_ASM3630560v1_genomic in GCA_981806595.1_ERR14205218_bin_42192_genomic
roary %>% 
  filter(!genome %in% c('GCA_036305605.1_ASM3630560v1_genomic', 
                        'GCA_981806595.1_ERR14205218_bin_42192_genomic')) %>% 
  filter(PA_01 == 1) %>% 
  group_by(Gene) %>% 
  reframe(n_genomes = n_distinct(genome)) %>% 
  filter(n_genomes > as.integer(77*0.99))
  