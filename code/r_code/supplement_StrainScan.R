# StrainScan 

library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)

theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 11),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 11)))


metadata <- read.table('metadata.tsv', sep = '\t', header = T)

ids  <- sprintf("%02d", 1:13)
paths <- paste0("strainscan/out_MH0", ids, "/final_report.txt")

# keep only files that actually exist
paths_ok <- paths[file.exists(paths)]
ids_ok   <- ids[file.exists(paths)]

strainscan_pre <- map_dfr(ids_ok, ~
                        read.delim(paste0("strainscan/out_MH0", .x, "/final_report.txt")) %>% 
                        mutate(sampleID = paste0("MH0", .x)))

strainscan <- strainscan_pre %>% 
  mutate(strain_id = ifelse(!is.na(Strain_ID), Strain_ID, ID), 
         relabund = ifelse(!is.na(Relative_Abundance_Inside_Cluster), Relative_Abundance_Inside_Cluster, Relative_Abundance), 
         Strain_Name = str_remove_all(Strain_Name, '_hybrid'), 
         strain_id = ifelse(strain_id %in% c(1, 8), '7', 
                            ifelse(strain_id %in% c(2, 4), '6', strain_id))) %>% 
  select(-c(Strain_ID, ID, Relative_Abundance_Inside_Cluster, Relative_Abundance)) %>%  
  group_by(sampleID, strain_id) %>% 
  reframe(relabund = sum(relabund)) %>% 
  left_join(metadata, by = 'sampleID')


strainscan %>% 
  ggplot(aes(x = days_from_first_collection, y = relabund, color = as.factor(strain_id))) +
  geom_line(linewidth= 2) +
  geom_point(size = 3) +
  labs(x = 'Day', y = 'Relative abundance', color = 'StrainScan Strains', subtitle = 'Participant H')

ggsave('out/strainscan/rel_abund_time.png', dpi=600)
