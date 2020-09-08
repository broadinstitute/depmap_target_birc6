library(tidyverse)
library(useful)
library(magrittr)
library(ggplot2)
library(here)

source(here::here('src', 'helper_functions.R'))

Achilles.gene.effect <- data.table::fread(here::here('data', 'raw', 'public-20q3_v32-achilles-gene-effect.csv')) %>% column_to_rownames('V1') %>% extract_hugo_symbol_colnames() 



genes_of_interest <- c('BIRC6', 'UBR4', 'UBA6', 'KCMF1')
gene_effect_subset <- Achilles.gene.effect[,genes_of_interest]

data.frame(BIRC6_cluster_mean_CERES = rowMeans(gene_effect_subset, na.rm = T), BIRC6_cluster_z_score_CERES = rowMeans(scale(gene_effect_subset), na.rm = T), gene_effect_subset) %>% 
  rownames_to_column('Row.name') %>% 
  write_csv(here::here('data', 'processed', 'cluster_level_CRISPR_ge.csv'))




# RNAi

RNAi.gene.effect <- data.table::fread(here::here('data', 'raw', 'demeter2-achilles_v12-gene-effect.csv')) %>% column_to_rownames('V1') %>% extract_hugo_symbol_colnames() 



genes_of_interest <- c('BIRC6', 'UBR4', 'UBA6', 'KCMF1')
gene_effect_subset <- RNAi.gene.effect[,genes_of_interest]

data.frame(BIRC6_cluster_mean_CERES = rowMeans(gene_effect_subset), BIRC6_cluster_z_score_CERES = rowMeans(scale(gene_effect_subset)), gene_effect_subset) %>% 
  rownames_to_column('Row.name') %>% 
  write_csv(here::here('data', 'processed', 'cluster_level_RNAi_ge.csv'))
