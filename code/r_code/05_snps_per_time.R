# How many SNPs are present in each strain group? 
# strain H007_176 je osamelec, ne spada v nobeno group! 

strain_groups <- read.csv('04_ANI/strain_groups.csv') %>% 
  mutate(sample = substr(isolate, 13, 20)) %>% 
  select(sample, group) 
  
snps <- read.table('snippy/gubbins/gubb.per_branch_statistics.csv', sep = '\t', header = TRUE)

strains_snps <- left_join(strain_groups, snps, by = join_by('sample' == 'Node'))

strains_snps %>%  
  group_by(group) %>%  
  reframe(no_snps = sum(Total.SNPs))

strains_snps %>% 
  mutate(timepoint = as.numeric(substr(sample, 2,4))) %>% 
  group_by(timepoint) %>% 
  reframe(no_snps = sum(Total.SNPs))


# Before Gubbins
snps <- read_tsv('~/projects/C_saudiense/snippy/core.tab') %>% 
  pivot_longer(names_to = 'sample', values_to = 'ALT', cols = starts_with('H00')) %>%  
  # if there was a SNP 1, otherwise 0
  mutate(what = ifelse(REF == ALT, 0, 1), 
         timepoint = as.numeric(substr(sample, 2, 4))) %>%  
  left_join(strain_groups, by = 'sample') %>% 
  filter(!is.na(group))

snps %>%  
  group_by(group) %>%  
  reframe(no_snps = sum(REF != ALT))

snps %>% 
  group_by(timepoint) %>% 
  reframe(no_snps = sum(REF != ALT))


# In MEGA I computete the tree and pairwise distances based on SNPs 
dist <- read_csv('snippy/MEGA/DistanceData.csv', col_names = TRUE) %>% 
  rename('isolate1' = '...1') %>% 
  pivot_longer(names_to = 'isolate2', values_to = 'dist', cols = starts_with('H00')) %>% 
  filter(!is.na(dist))

dist 

snippy_gubb_dists <- read_tsv('snippy/snippy_gubb_dists.tab') %>%  
  rename('isolate1' = `snp-dists 0.8.2`) %>% 
  pivot_longer(names_to = 'isolate2', values_to = 'dist', cols = starts_with('H00')) %>% 
  filter(!is.na(dist)) %>%  
  filter(isolate1 != isolate2)

# Define clusters (adjust intervals as needed)
breaks <- seq(0, max(snippy_gubb_dists$dist, na.rm=TRUE) + 50, by=50)

# Cut distances into clusters
snippy_gubb_dists <- snippy_gubb_dists %>%
  mutate(cluster = cut(dist, breaks = breaks, right=FALSE, include.lowest=TRUE))

# Calculate frequency and percent per cluster
freq_table <- snippy_gubb_dists %>%
  group_by(cluster) %>%
  summarise(freq = n()) %>%
  mutate(percent = freq / sum(freq) * 100)

# Plot frequency bar chart
ggplot(freq_table, aes(x=cluster, y=percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("SNP distance cluster") + ylab("Percentage of pairs (%)") +
  ggtitle("SNP distance clustering frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('out/snippy/clusters_snps.png')
