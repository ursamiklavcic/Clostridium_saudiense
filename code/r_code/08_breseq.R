# breseq results 
# code from https://github.com/barricklab/breseq/wiki/Tutorial-Clones


library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

theme_set(theme_bw(base_size=14))

compare <- read_csv('03_breseq_results/compare.csv')


count <- read_csv('03_breseq_results/count.csv') %>% 
  mutate(sample = substr(sample, 1, 8))

count %>% 
  ggplot(aes(x = as.factor(time), y = total)) +
  geom_boxplot() +
  geom_point(size = 3)

count_long <- count %>% 
  select(sample, STRAINGROUP, base_substitution, small_indel, 
         large_deletion, large_insertion, large_amplification, large_substitution, 
         mobile_element_insertion, gene_conversion, inversion, total) %>% 
  pivot_longer(cols = c(base_substitution, small_indel, large_deletion, large_insertion, large_amplification, large_substitution, mobile_element_insertion, gene_conversion, inversion),
               names_to = "what",
               values_to = "value" ) 
write_csv2(count_long %>%  filter(value > 0) %>% 
             pivot_wider(names_from = 'what', values_from = 'value'), '03_breseq_results/count_long.csv')

count_long %>% 
  filter(value > 0) %>% 
  ggplot(aes(x = value, y = sample, fill = what)) +
  geom_col() +
  facet_wrap(~STRAINGROUP, scales = 'free')

# Which proteins etc are SNPs in?
genes_strains <- compare %>% 
  mutate(sample = substr(title, 1, 8)) %>% 
  left_join(select(count, sample, STRAINGROUP), by = 'sample') %>% 
  select(gene_product, sample, STRAINGROUP) %>%  
  group_by(gene_product, STRAINGROUP) %>% 
  reframe(n = n_distinct(sample)) %>% 
  pivot_wider(names_from = 'gene_product', values_from = 'n')
write_csv2(genes_strains, '03_breseq_results/genes_per_strains.csv')  

long <- compare %>% 
  mutate(sample = substr(title, 1, 8)) %>% 
  left_join(select(count, sample, STRAINGROUP), by = 'sample') %>% 
  select(gene_product, sample, STRAINGROUP) 

long_genes <- long %>% 
  mutate(genes = ifelse(grepl("hypothetical", gene_product), 'Hypothetical protein', gene_product)) %>%
  group_by(genes) %>% 
  reframe(n = n_distinct(sample))

compare %>% 
  group_by(mutation_category) %>% 
  reframe(n = n_distinct(gene_name))

# SNPs 
