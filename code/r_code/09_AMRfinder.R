# AMRfinder+ results 

amrs <- read.table('007_AMRfinder+/AMRfinder_results_C_saudiense.tsv', sep = '\t', header = TRUE) %>% 
  filter(Contig.id != 'Contig id') %>% 
  rename('genome' = 'H005_102.fasta')
