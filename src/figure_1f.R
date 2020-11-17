library(tidyverse)
library(magrittr)
library(here)
library(ggplot2)

source(here::here('src', 'helper_functions.R'))

## Achilles

sample_info <- read_csv(here::here('data', 'raw', 'public-20q3_v32-sample-info.csv'))
EM_classification <- read_csv(here::here('data', 'processed', 'lineage_EM_classification.csv')) %>%
    mutate(classification= replace(classification, classification=='E', 'Epithelial')) %>%
    mutate(classification= replace(classification, classification=='M', 'Mesenchymal')) 

Achilles.gene.effect.CRISPR <- data.table::fread(here::here('data', 'processed', 'cluster_level_CRISPR_ge.csv')) %>% 
  select(DepMap_ID=Row.name, BIRC6_cluster_mean_CERES) %>%
  left_join(sample_info %>% select(DepMap_ID, lineage))

Achilles.gene.effect.lineage.medians <- Achilles.gene.effect.CRISPR %>% 
  left_join(sample_info %>% select(DepMap_ID, lineage)) %>%
  group_by(lineage) %>%
  summarise(lineage_median=median(BIRC6_cluster_mean_CERES), num_samples = n())

plot_data <- left_join(Achilles.gene.effect.CRISPR, Achilles.gene.effect.lineage.medians, by = 'lineage') %>% 
  mutate(lineage_mod=str_replace_all(lineage, '_', ' ')) %>%
  mutate(labels=paste0(lineage_mod, ' (', num_samples,')') %>% str_to_title()) %>%
  left_join(EM_classification, by=c('lineage'='primary_disease'))

ggplot(plot_data %>% filter(!is.na(classification)), aes(y=reorder(labels, desc(lineage_median)), x=BIRC6_cluster_mean_CERES, fill=classification))+
  scale_fill_manual(values=c('Epithelial'=rgb(240/255, 147/255,78/255), 'Mesenchymal'=rgb(31/255,165/255,78/255)))+
  geom_boxplot(alpha=1)+
  geom_point(alpha=0.2)+
  scale_y_discrete(position = "right")+
  theme_classic()+
  ylab('Primary site')+
  xlab('Gene effect (mean CERES)\n on the BIRC6 module')


ggsave(here::here('output', 'Fig1', '1f_disease_level_BIRC6_cluster_CERES_scores.pdf'), plot, dpi=600)
