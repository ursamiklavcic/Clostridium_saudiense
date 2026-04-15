# loading packages
library(readr)
library(dplyr)
library(ggplot2)
library(cutpointr)

# Reading the metadata
md <- read_csv2('../longitudinal_shotgun/data/metadata.csv') %>%  
  mutate(days_from_first_collection = day)

# reading the tsv table of pairwise distances
nGD <- read_tsv(file = "t__SGB6178_nGD.tsv", col_names = F, show_col_types = F) %>% 
  filter(X1 != 'reference') %>% 
  filter(X2 != 'reference')

# Adding the metadata to the table
nGD <- left_join(nGD %>% select(sampleID_1 = X1, everything()),
                 md %>% select(sampleID_1 = Group,
                               subjectID_1 = person,
                               days_from_first_collection_1 = days_from_first_collection))
nGD <- left_join(nGD %>% select(sampleID_2 = X2, everything()),
                 md %>% select(sampleID_2 = Group,
                               subjectID_2 = person,
                               days_from_first_collection_2 = days_from_first_collection))
# Computing time difference between sample (important for longitudinal samples)
nGD$time_day_diff <- abs(nGD$days_from_first_collection_1 - nGD$days_from_first_collection_2)

# Annotating pairs of samples. Are they from the same individual?
nGD$same_individual <- ifelse(nGD$subjectID_1 == nGD$subjectID_2, "same_individual", "different_individual")



# Keeping only the training data
nGD_training <- nGD # %>%
  filter(subjectID_1 %in% c("D","H"), 
         subjectID_2 %in% c("D","H"))

  

nGD_training <- rbind(nGD_training %>% 
                        filter(same_individual == "same_individual") %>%
                        filter(time_day_diff <= 180) %>%
                        group_by(subjectID_1) %>% 
                        arrange(subjectID_1, time_day_diff) %>% 
                        slice_head(n = 1) %>% 
                        ungroup(),
                      nGD_training %>% 
                        filter(same_individual == "different_individual") %>%
                        group_by(subjectID_1, subjectID_2) %>% 
                        slice_head(n = 1) %>% 
                        ungroup())


table(nGD_training$same_individual)

# Highly-prevalent SGBs Youden's index
res_youden <- cutpointr(data = nGD_training, 
                        x = X3, class = same_individual,
                        method = maximize_metric, metric = youden)
quantile_3pc <- nGD_training %>% 
  filter(same_individual == "different_individual") %>% 
  pull(X3) %>% quantile(0.03)

# Visualise distributions 

ggplot(data = nGD_training) +
  geom_density(aes(x = X3, fill = same_individual), alpha = 0.5 ) +
  geom_vline(aes(xintercept = res_youden$optimal_cutpoint, color = "youden")) +
  geom_vline(aes(xintercept = quantile_5pc, color = "quantile"), linetype = "dotted", lwd = 1) +
  theme_bw() + 
  labs(x = "StrainPhlAn nGD", y = "") +
  theme(legend.title =  element_blank(),
        legend.position = "bottom") +
  scale_color_manual(name = "statistics", values = c(youden = "blue", quantile = "red"))

prop.table(table(nGD_training %>% 
                   filter(same_individual == "different_individual") %>% 
                   pull(X3) >= res_youden$optimal_cutpoint))

#  TRUE
# 1



# metadata forward 
meta <-  md %>% select(sampleID = Group,
                       subjectID = person,
                       days_from_first_collection = days_from_first_collection)
write_tsv(meta, 'metadata.tsv')
