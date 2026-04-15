# Valles-Colomer et al. (2023)

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggtree)
library(ape)
library(Biostrings)

# StrainPhlan Tree
tree <- read.tree('09_strainPhlan/RAxML_bestTree.t__SGB6178.StrainPhlAn4.tre')

plot(tree)

metadata <- read.table("09_strainPhlan/metadata.tsv", sep = '\t', header = TRUE) %>%  
  rbind(data.frame(sampleID = 'reference', 
                   subjectID = 'reference', 
                   realtion  = 'reference', 
                   timepoint = 0))

# 
alignment <- readDNAStringSet("09_strainPhlan/t__SGB6178.StrainPhlAn4_concatenated.aln")
alignment_mat <- as.matrix(alignment)
rownames(alignment_mat) <- names(alignment)

is_variable <- apply(alignment_mat, 2, function(col) length(unique(col)) > 1)
snp_mat <- alignment_mat[, is_variable, drop = FALSE]

# Convert matrix to tibble for ggplot or other uses
snp_df <- as.data.frame(snp_mat) %>% 
  tibble::rownames_to_column("Sample") %>% 
  pivot_longer(-Sample, names_to = "POS", values_to = "base")

# Now snp_df contains a dataframe with columns Sample, POS, base
# You can use snp_df to plot SNP heatmaps or for further SNP analysis
# SNPs
variable_snps <- snp_df %>% 
  filter(Sample %in% c('MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010', 'MH011', 'MH012', 'MH013'), 
         base != '-') %>%
  group_by(POS) %>% 
  reframe(n = n_distinct(base)) %>% 
  filter(n > 1) %>%  
  pull(POS)

snp_df %>%  left_join(metadata, by = join_by('Sample' == 'sampleID')) %>%  
  filter(subjectID == 'H' & Sample %in% c('MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010', 'MH011', 'MH012', 'MH013'), 
         POS %in% variable_snps) %>%
  #mutate(POS = as.numeric(substr(POS, 2, 10))) %>%
  ggplot(aes(x = POS, y = Sample, fill = base)) +
  geom_tile() +
  scale_fill_manual(values = c("A" = "red", "C" = "yellow", "G" = "green", "T" = "blue", '-' = 'white')) +
  labs(x = 'Position')  
ggsave('out/strainphlan/snps_MH.png')
# This is not ok - tooo many snps 

# All 
ggtree(tree) %<+% metadata + 
  geom_tippoint(aes(color=factor(subjectID)), size=3) +
  geom_tiplab() +
  theme(legend.position = "bottom")

samples_H <- metadata %>%
  filter(subjectID == "H") %>%
  pull(sampleID)

# Prune tree to keep only tips for subject "H"
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, samples_H))

# Plot
ggtree(pruned_tree) %<+% filter(metadata, subjectID == 'H') + 
  geom_tippoint(aes(color=factor(timepoint)), size=3) +
  geom_tiplab() +
  theme(legend.position = "bottom") +
  labs(color = 'Timepoint')
ggsave('out/strainphlan/tree_H.png')

# 
# 
# # Stran tracking 
# H <- c('MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010', 'MH011', 'MH012', 'MH013', 
#        'SH001', 'SH002', 'SH003', 'SH004', 'SH005', 'SH006', 'SH007', 'SH008', 'SH009', 'SH010', 'SH011', 'SH012', 'SH013')
# 
# track <- read.table('09_strainPhlan/transmission_events.info', sep = '\t', header = FALSE) %>%  
#   mutate(sample1 = substr(V1, 1, 5), 
#          sample2 = substr(V1, 11, 16)) %>% 
#   select(-V1) 
# 
# trackH <- track %>%  
#   filter(sample1 %in% H & sample2 %in% H) %>%  
#   mutate(timepoint1 = as.numeric(substr(sample1, 3, 5)),
#          timepoint2 = as.numeric(substr(sample2, 3, 5)))

# Mutattion rates 

mutation_pre <- read.table('09_strainPhlan/t__SGB6178.mutation', sep = '\t', header = T) %>% 
  filter(ids %in% c('MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010', 'MH011', 'MH012', 'MH013', 
                    'SH001', 'SH002', 'SH003', 'SH004', 'SH005', 'SH006', 'SH007', 'SH008', 'SH009', 'SH010', 'SH011', 'SH012', 'SH013')) %>% 
  column_to_rownames('ids') %>%  
  select(c('MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH007', 'MH008', 'MH009', 'MH010', 'MH011',
           'SH001', 'SH002', 'SH003', 'SH005', 'SH008', 'SH009', 'SH010', 'SH011', 'SH012'))

# Function to convert fraction string to decimal numeric
convert_fraction_to_decimal <- function(x) {
  if (grepl("/", x)) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  } else {
    as.numeric(x)
  }
}

# Apply the function to the entire data frame (non-rowname part)
mutation <- mutation_pre %>%
  mutate(across(everything(), ~ sapply(., convert_fraction_to_decimal))) %>% 
  rownames_to_column('sample1') %>%  
  pivot_longer(names_to = 'sample2', values_to = 'mutation_rate', cols = c(-sample1)) %>%  
  filter(sample1 != sample2)


mutation %>%  
  ggplot(aes(x = reorder(sample1, mutation_rate), y = reorder(sample2, mutation_rate), fill = mutation_rate)) +
  geom_tile()


