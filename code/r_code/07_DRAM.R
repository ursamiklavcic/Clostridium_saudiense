# DRAM 


library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(stringr)

summaries <- read_xlsx('005_DRAM/distill/metabolism_summary.xlsx') %>%  
  pivot_longer(names_to = 'samples', values_to = 'PA', cols = starts_with('H00'))

metadata <- read.table('metadata_genomes.tsv', header = T)

strains <- full_join(summaries, metadata, by = 'samples') %>%  
  mutate(samples = str_remove(samples, '_hybrid'))

strains_genes <- strains %>% 
  group_by(module, samples, TIMEPOINT, ANTIBIOTIC_THERAPY, STRAIN_GROUP, STRAINS) %>% 
  reframe(PA = sum(as.numeric(PA), na.rm = T))

strains %>%  
  filter(!is.na(module), PA > 0) %>% 
  ggplot(aes(x = samples, y = module, fill = as.numeric(PA))) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('out/DRAM/modules_genome.svg', dpi=600)
