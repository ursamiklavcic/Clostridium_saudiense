library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ape)

theme_set(theme_bw())

#core <- read_tsv('~/projects/C_saudiense/snippy/core.tab')
core <- read_tsv('~/projects/C_saudiense/snippy_H005_98_hybrid/core.tab')
metadata <- read_tsv('~/projects/C_saudiense/snippy/metadata.tsv') 

snp_matrix <- core %>% 
  pivot_longer(names_to = 'sample', values_to = 'ALT', cols = starts_with('H00')) %>%  
  # if there was a SNP 1, otherwise 0
  mutate(what = ifelse(REF == ALT, 0, 1)) %>% 
  pivot_wider(names_from = 'sample', values_from = 'what', values_fill = 0) %>% 
  select(-c(CHR, REF, ALT, POS)) %>%
  as.matrix() %>% 
  t()

snp_dist <- dist(snp_matrix, method = 'euclidian')

# Visualise distance matrix 
snp_dist_plot <- as.data.frame(as.matrix(snp_dist)) %>% 
  rownames_to_column('sample1') %>% 
  pivot_longer(names_to = 'sample2', values_to = 'euclidian_dist', cols = starts_with('H00')) 


ggplot(snp_dist_plot, aes(x = reorder(sample1, euclidian_dist), y = reorder(sample2,euclidian_dist-1), fill = euclidian_dist)) +
  geom_tile() +
  scale_fill_gradient2(low = 'darkgreen',mid = "white",high = 'lightgreen', midpoint = 30) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = '', y = '', fill = 'Euclidian dist based on \nSNPs from SNPpy')
ggsave('~/projects/C_saudiense/out/snippy/clusters.png')

# Visualise dendrogram of SNPs 
annotation <- data.frame(Timepoint = as.factor(metadata$timepoint))
rownames(annotation) <- metadata$name

pheatmap(as.matrix(snp_matrix),
         clustering_distance_rows = snp_dist,
         clustering_method = "ward.D2",
         annotation_row = annotation,
         show_rownames = TRUE, 
         treeheight_row = 100,  # default is 50
         treeheight_col = 100)

# herarcichal clustering 
hc <- hclust(snp_dist, method = 'ward.D2')
# We can observe 4 clusters based on SNPs from pheatmap
clusters <- cutree(hc, k = 4)

strain_df <- data.frame(name = rownames(snp_matrix), 
                        strain = as.factor(clusters)) %>% 
  left_join(metadata, by = 'name')

strain_df %>% 
  ggplot(aes(x = as.factor(timepoint), fill = strain)) +
  geom_bar(position = 'dodge') +
  labs(y = '# strains', x = 'Timepoint', fill = 'Strain')
ggsave('~/projects/C_saudiense/out/snippy/strains_timepoints.png')

strain_df2 <- strain_df %>% 
  group_by(timepoint) %>%  
  reframe(n_strain = n_distinct(strain), 
          n_sample = n_distinct(name)) 

strain_df2 %>% 
  ggplot(aes(x = as.factor(timepoint), y = n_strain)) +
  geom_col() +
  labs(x = 'Timepoint', y = '# strains')
ggsave('~/projects/C_saudiense/out/snippy/strains_timepoints_count.png')

strain_df2 %>%  
  ggplot(aes(x = as.factor(timepoint), y = n_sample)) +
  geom_col() +
  labs(x = 'Timepoint', y = '# isolates')
ggsave('~/projects/C_saudiense/out/snippy/isoaltes_timepoints.png')

# Strain clusters thorugh time 
library(networkD3)

# First, build three levels: person, timepoint, strain
links <- bind_rows(strain_df %>%
    count(person, timepoint, name = "value") %>%
    mutate(source = person, target = paste0(" ", timepoint)) %>%
    select(source, target, value),
  strain_df %>%
    count(timepoint, strain, name = "value") %>%
    mutate(source = paste0(" ", timepoint),
           target = paste0("Strain_", strain)) %>%
    select(source, target, value))

# Make nodes table (all unique names)
nodes <- data.frame(name = unique(c(links$source, links$target)))

# Map source/target to IDs
links <- links %>%
  mutate(source = match(source, nodes$name) - 1,
    target = match(target, nodes$name) - 1 )

# Plot Sankey
sankeyNetwork(Links = links,
  Nodes = nodes,
  Source = "target",
  Target = "source",
  Value = "value",
  NodeID = "name",
  sinksRight = FALSE,
  fontSize = 12,
  nodeWidth = 30)


# Look at overall statistics for variant calling 
core_txt <- read.table('projects/C_saudiense/snippy/core.txt', header =  T) %>%  
  left_join(metadata, by = join_by('ID' == 'name')) 

core_txt %>%  
  filter(!is.na(timepoint)) %>% 
  ggplot(aes(x = (UNALIGNED/LENGTH)*100, y = (ALIGNED/LENGTH)*100)) +
  geom_text(aes(label = ID), show.legend = FALSE, angle = 45) +
  geom_point(size=3, aes(color = as.factor(timepoint))) +
  labs(x = '% unaligned to the reference', y = '% aligned to the reference')
ggsave('~/projects/C_saudiense/out/snippy/aligned_unaligned.png')

core_txt2 <- filter(core_txt, !is.na(timepoint))
min(core_txt2$VARIANT)
max(core_txt2$VARIANT)
mean(core_txt2$VARIANT)
median(core_txt2$VARIANT)

core_txt %>%  
  ggplot(aes(x = ID, y = VARIANT, color = as.factor(timepoint))) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 90) )
ggsave('~/projects/C_saudiense/out/snippy/VARIANT_2.png')


min(core_txt2$LOWCOV)
max(core_txt2$LOWCOV)
mean(core_txt2$LOWCOV)
median(core_txt2$LOWCOV)

core_txt %>%  
  ggplot(aes(x = ID, y = LOWCOV, color = as.factor(timepoint))) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 90) )
ggsave('~/projects/C_saudiense/out/snippy/LOWCOV_2.png')


# Plot the distibution of vriants across the genome 

hist(core$POS)
ggsave('~/projects/C_saudiense/out/snippy/position_SNPs.png')

ggplot(core, aes(x = POS)) +
  geom_histogram(binwidth=100) 
ggsave('~/projects/C_saudiense/out/snippy/position_SNPs_100.png')

# Phylogenetic tree annotation and visualisation https://comparative-genomics.readthedocs.io/en/latest/day3_morning.html

# Phylogenetic tree
tree <- read.tree('~/projects/C_saudiense/snippy/core.tree')

plot(tree)

# Reroot tree in reference (outgroup)
tree_reroted <- root(tree, 'Reference')
save(tree_reroted, file = '~/projects/C_saudiense/snippy/core_reroted.tree')

plot(tree_reroted)


##
# Exclude polymorphic regions of the genome and do the SNP analysis again!
# Phylogenetic tree
tree_gub <- read.tree('snippy/gubbins/gubb.final_tree.tre')

plot(tree_gub)


# Distance matrix after removal of polymorhic sites with Gubbins
clean_aln <- read.dna('snippy/Csaudiense_clean.core.fasta', format = 'fasta')

strain_groups = read_csv('04_ANI/strain_groups.csv') %>% 
  mutate(sample = substr(isolate, 13, 21)) %>% 
  select(sample,group)

# Poisitions character matrix 
msa_char <- as.character(clean_aln) %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  pivot_longer(-sample, names_to = 'POS', values_to = 'base')

reference <- filter(msa_char, sample == 'Reference') %>%  
  select(-sample)

msa_long <- full_join(filter(msa_char, sample != 'Reference'), reference, by = 'POS') %>%  
  mutate(is_change = ifelse(base.y != base.x, TRUE, FALSE)) %>% 
  left_join(strain_groups, by = 'sample') %>% 
  mutate(group = ifelse(is.na(group), 5, group), 
         timepoint = as.numeric(substr(sample, 2, 4)), 
         isolate = substr(sample, 6, 8)) 

# I have to check where all my isolates have a chnage to the reference = and exclude this snps
msa_long %>%  group_by(POS) %>% 
  reframe(n_base = list(unique(base.x)), 
          ref = base.y) %>% 
  unnest(n_base) %>% 
  unique() %>% 
  group_by(POS) %>% 
  mutate(n_b = n_distinct(n_base), 
         n_ref = n_distinct(ref), 
         t = ifelse(n_base == n_ref, TRUE, FALSE )) %>%  
  filter(t == TRUE, n_b == 1)

msa_long %>% 
  filter(is_change = TRUE) %>% 
  group_by(timepoint) %>% 
  reframe(n_isolates = n_distinct(isolate), 
          n_snps = sum(is_change == TRUE), 
          n_isolates_per_timepoint = n_snps/n_isolates)

msa_long %>% 
  filter(is_change = TRUE) %>% 
  group_by(group) %>% 
  reframe(n_isolates = n_distinct(isolate), 
          n_snps = sum(is_change == TRUE), 
          n_isolates_per_strainGroup = n_snps/n_isolates)


# 
base_map <- c('a' = 1, 'c' = 2, 'g' = 3, 't'  = 4)

msa_long %>% 
  filter(is_change == TRUE) %>% 
  ggplot(aes(x = sample, y = POS, fill = base.x)) +
  geom_tile() +
  scale_fill_manual(values = c('a' = "#4daf4a", 'c' = "#F5EB27", 'g' = "#F59827", 
                               't' = "#984ea3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank()) +
  labs(fill = 'Base', x = 'Position on genome', y = '')
ggsave('out/snippy/base_changes.png')

# Trying somethning different 
# https://cran.r-project.org/web/packages/adegraphics/vignettes/adegraphics.html
library(adegenet)
snp_data <- fasta2genlight('snippy/Csaudiense_clean.core.fasta')

# Distribution of SNPs in the genome 
snpposi.plot(position(snp_data), genome.size=3707884, codon=FALSE)
ggsave('out/adegenet/SNP_genome_distribution.png')

snpposi.plot(position(snp_data), genome.size=7500, codon=TRUE)
ggsave('out/adegenet/SNP_genome_distribution_closeup.png')

# SNPs? Where?
glPlot(snp_data)

snp_plot <- as.data.frame(as.matrix(snp_data)) %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(-sample, names_to = "SNP", values_to = "Allele") %>%  
  mutate(REF = substr(SNP, 3, 3), 
         ALT = substr(SNP, 5, 5))

df_clean <- as.data.frame(as.matrix(snp_data)) %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(-sample, names_to = "SNP", values_to = "Allele")  %>% 
  separate(SNP, into = c("SNP_index", "Allele_pair"), sep = "\\.") %>%
  mutate(SNP_index = as.integer(SNP_index)) %>%
  separate(Allele_pair, into = c("REF", "ALT"), sep = "/")

# Rearrange dr_clean so that strain groups are togehter 
df_clean$sample = factor(df_clean$sample, levels = c('H005_98', 'H005_102', 
                                                     'H009_220', 'H009_215', 'H009_221', 'H009_213', 'H009_218', 
                                                     'H007_176', 
                                                     'H005_75', 'H005_73', 'H005_77', 'H005_97', 
                                                     'H009_230', 'H009_234', 'H009_256', 
                                                     'H007_145', 'H009_211', 'H009_236', 'H009_251', 'H009_202', 'H009_259', 'H009_255', 'H009_228', 
                                                     'H009_205', 'H009_219', 'H005_118', 'H009_217', 'H007_194', 'H007_154', 'H005_137', 'H007_157', 
                                                     'H007_187', 'H007_148', 'H007_158', 'H007_181', 'H007_146', 'H007_183', 'H007_156', 'H007_190', 
                                                     'H007_162', 'H007_172', 'H007_188'))


ggplot(df_clean %>% filter(sample != 'Reference'), aes(x = SNP_index, y = sample, fill = factor(Allele))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "2" = "purple")) +
  labs(x = "SNP", y = "", fill = "Number of 2nd allele") +
  theme(axis.text.y = element_text(size = 8))
ggsave('out/adegenet/SNPs_index.png')


df_clean %>%  filter(sample != 'Reference') %>% 
  mutate(change = ifelse(REF != ALT, 'diff', 'same')) %>%  
  filter(change == 'diff') %>% 
  ggplot(aes(x = SNP_index, y = sample, fill = ALT)) +
  geom_tile() +
  scale_fill_manual(values = c('A' = "#4daf4a", 'C' = "#F5EB27", 'G' = "#F59827", 
                                'T' = "#984ea3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = 'Base', x = 'Position on genome', y = '') 


# Plot second allele frequencies
# glMean computes the mean of second alleles (frequencies for each SNP)
myFreq <- glMean(snp_data)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)
ggsave('out/adegenet/sec_allele_freq.png')

# Plot main allele frequencies
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
ggsave('out/adegenet/allele_freq.png')

# PCA 
glPca(snp_data) # we keep the first 3 PCAs 
pca1 <- glPca(snp_data, nf = 2) 

pca_df <- pca1$scores %>%  
  as.data.frame() %>% 
  rownames_to_column('sample') %>%  
  mutate(timepoint = as.numeric(substr(sample, 2, 4)))

library(ggrepel)
ggplot(pca_df %>% filter(sample != 'Reference'), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = as.factor(timepoint)), size = 4) +
  geom_text_repel(aes(label = sample), size = 4, fontface = 'bold',
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  labs(color = 'Timepoint')
ggsave('out/adegenet/PCA.png')

