library(tidyverse)
library(magrittr)
library(here)
library(ggplot2)

source(here::here('src', 'helper_functions.R'))

## Achilles

sample_info <- read_csv(here::here('data', 'raw', 'public-20q3_v32-sample-info.csv'))
EM_classification <- read_csv(here::here('data', 'feature_engineering', 'EMT_cell_classification_TS.csv'))

Achilles.gene.effect.CRISPR <- data.table::fread(here::here('data', 'processed', 'cluster_level_CRISPR_ge.csv')) %>% 
  select(DepMap_ID=Row.name, BIRC6_cluster_mean_CERES) %>%
  left_join(sample_info %>% select(DepMap_ID, primary_disease))

Achilles.gene.effect.lineage.medians <- Achilles.gene.effect.CRISPR %>% 
  left_join(sample_info %>% select(DepMap_ID, primary_disease)) %>%
  group_by(primary_disease) %>%
  summarise(lineage_median=median(BIRC6_cluster_mean_CERES))

plot_data <- left_join(Achilles.gene.effect.CRISPR, Achilles.gene.effect.lineage.medians, by = 'primary_disease')

plot <- ggplot(plot_data, aes(y=reorder(primary_disease, desc(lineage_median)), x=BIRC6_cluster_mean_CERES))+
  geom_boxplot(alpha=1)+
  geom_point(alpha=0.2)+
  scale_y_discrete(position = "right")+
  theme_classic()+
  ylab('Primary site')+
  xlab('Gene effect (mean CERES)\n on the BIRC6 module')

ggsave(here::here('output', 'Fig1', '1f_disease_level_BIRC6_cluster_CERES_scores.pdf'), plot, dpi=600)