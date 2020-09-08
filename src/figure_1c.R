library(tidyverse)
library(magrittr)
library(corrplot)
library(here)
library(ggpubr)

source(here::here('src', 'helper_functions.R'))


## Achilles

Achilles.gene.effect.CRISPR <- data.table::fread(here::here('data', 'public-20q3_v32-achilles-gene-effect.csv')) %>% column_to_rownames('V1') %>% extract_hugo_symbol_colnames() 

genes <- c('UBA6', 'BIRC6', 'KCMF1', 'UBR4')
ge_subset_crispr <- Achilles.gene.effect.CRISPR[, genes]
correlation = cor(ge_subset_crispr[rowSums(is.na(ge_subset_crispr))==0,])

pdf(file =here::here('output', 'Fig1', '1c_CRISPR_correlations.pdf'))
correlation_plot <- corrplot(correlation, method="circle", type='lower', addCoef.col='white')
dev.off()




plot_CRISPR_1 <- ggscatter(ge_subset_crispr, x = "UBA6", y = "BIRC6", add = "reg.line", conf.int = TRUE, color='blue', alpha=0.3,
          add.params = list(color = "black",
                            fill = "lightgray")
          )+
  stat_cor(method = "pearson",  label.x = 0, label.y = -1.5)+
  xlab('UBA6 dependency (CERES)')+
  ylab('BIRC6 dependency (CERES)')+
  theme_classic()

ggsave(here::here('output', 'Fig1', '1c_UBA6_BIRC6_CRISPR_correlation.pdf'), plot_CRISPR_1, dpi=600)

plot_CRISPR_2 <- ggscatter(ge_subset_crispr, x = "KCMF1", y = "UBR4", add = "reg.line", conf.int = TRUE, color='orange', alpha=0.5,
          add.params = list(color = "black",
                            fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 0, label.y = -1.7)+
  xlab('KCMF1 dependency (CERES)')+
  ylab('UBR4 dependency (CERES)')

ggsave(here::here('output', 'Fig1', '1c_KCMF1_BIRC6_CRISPR_correlation.pdf'), plot_CRISPR_2, dpi=600)




## RNAi

gene.effect.RNAi <- data.table::fread(here::here('data', 'demeter2-achilles_v12-gene-effect.csv')) %>% column_to_rownames('V1') %>% extract_hugo_symbol_colnames() 

ge_subset_RNAi <- gene.effect.RNAi[, genes]
correlation = cor(ge_subset_RNAi[rowSums(is.na(ge_subset))==0,])
corrplot(correlation, method="circle", type='lower', addCoef.col='white')

pdf(file =here::here('output', 'Fig1', '1c_RNAi_correlations.pdf'))
correlation_plot <- corrplot(correlation, method="circle", type='lower', addCoef.col='white')
dev.off()