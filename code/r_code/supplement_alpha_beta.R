# Alpha, beta, normalized abundance participant H 

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggpubr)
library(purrr)
library(stringr)

set.seed(96)
theme_set(theme_bw(base_size = 14))

metadata <- read.table('~/projects/thesis/data/metadata.csv', sep = ';', header = TRUE) %>%  
  filter(person == 'H') %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untretaed sample', biota))
  

alpha <- readRDS('~/projects/thesis/data/longitudinal_shotgun/alpha_diveristy.RDS') %>% select(name, richness, shannon, evenness) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(person == 'H')

beta <- readRDS('~/projects/thesis/data/longitudinal_shotgun/nmds_mpa_positions.rds')

events <- read.table('~/projects/thesis/data/extreme_event_data.csv', sep = ',', header = T) %>% filter(person == 'H')

alpha_plot <- ggplot(alpha, aes(x = day, y = shannon, color = biota)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336', '#3CB371')) +
  labs(x = 'Day', y = 'Shannon diveristy index', color = 'Sample type', fill = 'Event') +
  theme(legend.position = 'bottom')
alpha_plot

beta_plot <- beta %>% 
  filter(!is.na(person)) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size = 5) +
  labs(color = 'Individual') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom") 
beta_plot

normalized <- readRDS('~/projects/longitudinal_amplicons/data/r_data/long_all.RDS') %>% 
  filter(person== 'H') %>% 
  left_join(metadata, by = 'original_sample')

norm_plot <- normalized %>% 
  # group_by(day) %>% 
  # reframe(norm_abund = sum(norm_abund)) %>% 
  ggplot(aes(x = as.factor(day), y = norm_abund)) +
  # geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
  #           alpha = 0.6) +
  # scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = 'Day', y = 'Absolute bacterial abundance', fill = 'Event') +
  theme(legend.position = 'bottom')
norm_plot

ggplot(normalized, aes(x = as.factor(day), y = rel_abund)) +
  geom_boxplot() +
  scale_y_log10()


ggarrange(alpha_plot + labs(tag = 'A'), 
          beta_plot + labs(tag = 'B'), 
          widths = c(1, .8), nrow = 1)
ggsave('out/article/supplement_2.svg')


