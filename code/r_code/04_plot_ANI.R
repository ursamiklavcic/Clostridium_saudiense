# ANI 

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 11),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 12)))

ANIm <- read.table('002_ANI/ANIm_identity.tsv', sep = '\t', header = TRUE) %>%  
  pivot_longer(names_to = 'sample2', values_to = 'ANIm', cols = starts_with('H00')) %>%  
  rename(sample1= X)

ANIb <- read.table('002_ANI/ANIb_identity.tsv', sep = '\t', header = TRUE) %>%  
  pivot_longer(names_to = 'sample2', values_to = 'ANIb', cols = starts_with('H00')) %>%  
  rename(sample1= X)

sourmash <- read.table('002_ANI/sourmash_identity.tsv', sep = '\t', header = TRUE) %>%  
  pivot_longer(names_to = 'sample2', values_to = 'sourmash', cols = starts_with('H00')) %>%  
  rename(sample1= X)

ANI <- full_join(ANIm, ANIb, by = c('sample1', 'sample2')) %>%  
  full_join(sourmash, by = c('sample1', 'sample2')) %>%  
  filter(sample1 != sample2) 

# % genome coverage with which ANI was calculated 
ANI_coverage <- full_join(read.table('002_ANI/ANIb_query_cov.tsv', sep = '\t', header = TRUE) %>% 
                            pivot_longer(names_to = 'sample2', values_to = 'coverage_ANIb', cols = starts_with('H0')) %>% 
                            rename('sample1' = 'X'), 
                          read.table('002_ANI/ANIm_query_cov.tsv', sep = '\t', header = TRUE) %>% 
                            pivot_longer(names_to = 'sample2', values_to = 'coverage_ANIm', cols = starts_with('H0')) %>% 
                            rename('sample1' = 'X'), by = c('sample1', 'sample2')) %>%  
  full_join(read.table('002_ANI/sourmash_query_cov.tsv', sep = '\t', header = TRUE) %>% 
              pivot_longer(names_to = 'sample2', values_to = 'coverage_sourmash', cols = starts_with('H0')) %>% 
              rename('sample1' = 'X'), by = c('sample1', 'sample2'))
                          

# STRAIN BREAK 
ANIb_break <- ANI %>% 
  ggplot(aes(x = ANIb)) +
  geom_histogram(binwidth = 0.0001) +
  labs(x = 'ANIb percentage identity [%]', y = 'Pairwise comparisons') +
  coord_cartesian(xlim = c(0.99, 1))
ANIb_break
ggsave('out/ANI/ANIb_strain_break.png', dpi = 600)

ANI %>% 
  ggplot(aes(x = ANIm)) +
  geom_histogram(binwidth = 0.0001) +
  labs(x = 'ANIm percentage identity [%]', y = 'Pairwise comparisons') +
  coord_cartesian(xlim = c(0.99, 1))
ggsave('out/ANI/ANIm_strain_break.png', dpi = 600)

ANI %>% 
  ggplot(aes(x = sourmash)) +
  geom_histogram(binwidth = 0.0001) +
  labs(x = 'ANI sourmash percentage identity [%]', y = 'Pairwise comparisons') +
  coord_cartesian(xlim = c(0.99, 1))
ggsave('out/ANI/sourmash_strain_break.png', dpi = 600)


# ANI pairwise compariosns with genome coverage 
ANIb_coverage <- ANI %>% full_join(ANI_coverage, by = c('sample1', 'sample2')) %>% 
  filter(ANIb > 0.99, sample1 != sample2) %>% 
  ggplot(aes(y = coverage_ANIb, x = ANIb)) + 
  geom_point() +
  labs(x = 'ANIb percentage identity [%]', y = 'Genome coverage') +
  coord_cartesian(xlim = c(0.99, 1))
ANIb_coverage

# Article plot 
ggarrange(ANIb_break + geom_vline(xintercept = 0.9972, color = 'skyblue', linewidth = 2)
          + theme(axis.title.x = element_blank(), 
                                axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank()) +
            labs(tag = 'A'), 
          ANIb_coverage + labs(tag = 'B'), 
          nrow = 2, align = 'v', heights = c(0.9, 1))
ggsave('out/Frontiers/genome_coverage_ANIb.svg')

# Correlation between ANI and coverage ? 
ANI_both <- full_join(ANI, ANI_coverage, by = c('sample1', 'sample2')) %>% 
  filter(ANIb > 0.99)

cor.test(ANI_both$ANIb, ANI_both$coverage_ANIb, method = 'spearman')

# Strain brak at 0.998 % ANI
threshold <- 0.9972  
ANI_filt <- filter(ANI, ANIb > threshold)

library(igraph)
# Create an undirected graph from the pairs
g <- graph_from_data_frame(ANI_filt[, c("sample1", "sample2")], directed = FALSE)

# Find connected groups
components <- components(g)

# Create a dataframe mapping each genome to its group
strain_groups <- data.frame(
  sample = names(components$membership),
  STRAIN_GROUP_ANIb = components$membership)

# Save 
write.csv(strain_groups, '002_ANI/ANIb_strain_groups.csv')

# Strains through time with ANIb 
metadata <- read.table('metadata_genomes.csv', sep = '\t', header = TRUE)
ANIb_strains <- left_join(metadata, strain_groups, by = join_by('samples' == 'sample'))

ANIb_strains %>%  
  ggplot(aes(x = as.factor(TIMEPOINT), fill = as.factor(STRAIN_GROUP_ANIb))) +
  geom_bar() +
  labs(x = 'Time point', y = '# isolate genomes', fill = 'Strain group\nANIb')
ggsave('out/ANI/ANIb_strain_timepoints.png')

# Strain groups based on SNPs 
ANIb_strains %>%  
  ggplot(aes(x = as.factor(TIMEPOINT), fill = as.factor(STRAINS))) +
  geom_bar() +
  labs(x = 'Time point', y = '# isolate genomes', fill = 'Strain group\ncore SNPs')
ggsave('out/ANI/SNPs_strain_timepoints.png')

ANIb_strains %>%  
  ggplot(aes(x = as.factor(TIMEPOINT), fill = as.factor(STRAIN_GROUP))) +
  geom_bar() +
  labs(x = 'Time point', y = '# isolate genomes', fill = 'Strain group\ncore SNPs')
ggsave('out/ANI/SNPs_V2_strain_timepoints.png')

