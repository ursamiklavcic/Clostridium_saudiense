# metaphlan results
library(stringr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(tidyr)

theme_set(theme_bw(base_size = 14))

metadata <- read.table('~/projects/longitudinal_shotgun/data/metadata.csv', header= TRUE, sep = ';') %>%
  mutate(date = dmy(date))

abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain'),
           sep="\\|")

## Cultivated species Clostridium saudiense
cs <- abund %>% filter(Species == 's__Clostridium_saudiense') %>%
 # filter(Strain == 't__SGB6178') %>% 
  pivot_longer(cols = -c(Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain)) 

cs %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(biota == 'bulk microbiota') %>% 
  filter(value > 0) %>% 
  ggplot(aes(x = day, y = value, color = person)) +
  geom_point(size = 3, show.legend = F) +
  #geom_line(linewidth = 1, show.legend = F) +
  scale_y_log10() +
  #scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  facet_wrap(~person, scales = 'free_y') +
  labs(x = 'Day', y = 'Relative abundance', color = 'Individual')
ggsave('out/strainphlan/relative_abundance_microbiota.png')

cs %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(biota == 'ethanol treated sample') %>% 
  filter(value > 0) %>% 
  ggplot(aes(x = day, y = value, color = person)) +
  geom_point() +
  geom_line(linewidth=2) +
  scale_y_log10() +
  #scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  facet_wrap(~person, scales = 'free_y') +
  labs(x = 'Day', y = 'Relative abundance', color = 'Individual')
ggsave('out/strainphlan/relative_abundance_ethanol.png')
