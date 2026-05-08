# DRAM 

library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(stringr)

metadata <- read.table('data/metadata_genomes.tsv', header = T)

strains <- read_xlsx('005_DRAM/distill/metabolism_summary.xlsx') %>%  
  pivot_longer(names_to = 'samples', values_to = 'PA', cols = starts_with('H00')) %>% 
  mutate(PA = as.integer(as.numeric(PA)), 
         samples = str_remove(samples, '_hybrid')) %>%  
  filter(!samples %in% c('H005_94', 'H005_126')) %>% 
  left_join(metadata, by = 'samples')
  

strains %>% 
  group_by(STRAINS, module) %>% 
  reframe(PA = ifelse(sum(PA) > 0, 1, 0)) %>% 
  filter(PA > 0) %>% 
  ggplot(aes(x = as.factor(STRAINS), y = module, fill = as.factor(PA))) +
  geom_point(fill = 'black', alpha = .4, size = 5) +
  labs(x = 'Strains', y = 'Module') +
  theme_bw()
ggsave('out/DRAM/modules_genome.svg', dpi=600)
