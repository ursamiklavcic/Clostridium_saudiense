# MGE detection 

library(readxl)
set.seed(96)
theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 12),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 12)))



metadata <- read.table('data/metadata_genomes.tsv', sep = '\t', header = T)

# MobileElementFinder 
mgf <- read_excel('008_mge/zdruzeni_csv.xlsx') %>% 
  mutate(genomes = str_remove_all(Ime_datoteke, '.csv')) %>% 
  filter(!genomes %in% c('H005_94', 'H005_126')) %>% 
  left_join(metadata, by = join_by('genomes' == 'samples'))

mgf %>%  
  group_by(name) %>%  
  reframe(n = n_distinct(genomes))

mgf %>% 
  group_by(genomes, type) %>% 
  reframe(n = n_distinct(name), name, STRAINS) %>% 
  ggplot(aes(x = genomes, fill = name)) +
  geom_bar() +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  facet_wrap(~STRAINS, scales = 'free') +
  labs(y = 'Number of mobile genetic elements', x = '', fill = 'Name of mobile\ngenetic element') +
  coord_flip() 

ggsave('out/MGE/MobileElementFinder.svg', dpi=400)


# Mob Suite 
files <- list.files(
  path = '008_mge_mob_suite',
  pattern = ".mge_report.txt$",
  recursive = TRUE,
  full.names = TRUE)

# read and combine
mob <- map_dfr(files, function(f) {
  genome_name <- str_extract(f, "H[0-9]+_[0-9]+")
  
  read_tsv(f, col_types = cols(.default = "c"),  # force all columns as character
           progress = FALSE,
           show_col_types = FALSE) %>%
    mutate(genomes = genome_name) }) %>% 
  left_join(metadata, by = join_by('genomes' == 'samples')) 

mob %>% 
  # GenBank Accession of MGE
  group_by(mge_acs) %>%  
  reframe(n = n_distinct(genomes))
# Če imajo MGE imajo po MOB SUITE samo enega in si ga večina deli! 

mob %>%  
  ggplot(aes(x = as.factor(TIMEPOINT), fill = as.factor(STRAINS))) +
  geom_bar() +
  scale_y_continuous(breaks = c(0,2,4,6,8)) +
  labs(x = 'Time point', y = 'Number of genomes', fill = 'Strain')
ggsave('out/MGE/mob_suite.png', dpi=400)


