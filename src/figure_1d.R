library(tidyverse)
library(magrittr)
library(here)
library(ggplot2)

source(here::here('src', 'helper_functions.R'))

genes <- c('UBA6', 'BIRC6', 'KCMF1', 'UBR4')

## Achilles

Achilles.gene.effect.CRISPR <- data.table::fread(here::here('data', 'public-20q3_v32-achilles-gene-effect.csv')) %>% 
  column_to_rownames('V1') %>% 
  extract_hugo_symbol_colnames() 



ge_subset_crispr <- Achilles.gene.effect.CRISPR[, genes] %>% rownames_to_column(var='CLs') %>% 
  pivot_longer(-CLs, names_to='Gene', values_to='gene_effect')

corner(ge_subset_crispr)

colors <- c(rgb(232/255,203/255,215/255), 
            rgb(210/255,152/255,176/255),
            rgb(164/255,60/255,103/255),
            rgb(142/255,22/255,70/255))

ge_subset_crispr$Gene <- factor(ge_subset_crispr$Gene, levels = genes)

ggplot(ge_subset_crispr, aes(x=gene_effect))+
  geom_density(aes(fill=Gene), show.legend = F)+
  facet_wrap(~Gene, nrow=4, ncol=1)+
  geom_vline(xintercept = -1, col='red')+
  scale_fill_manual(values=colors)+
  theme_minimal()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab('Gene effect (CERES)')

ggsave(here::here('output', 'Fig1', '1d_CRISPR_gene_effect_density_plots.pdf'), plot_CRISPR_2, dpi=600)

Achilles.gene.dep.CRISPR <- data.table::fread(here::here('data', 'public-20q3_v32-achilles-gene-dependency.csv')) %>% 
  column_to_rownames('V1') %>% 
  extract_hugo_symbol_colnames() 

pr_subset_crispr <- Achilles.gene.dep.CRISPR[, genes]

sink(here::here('output', 'Fig1', '1d_CL_nums.txt'))

print('Strongly dependent CLs')
dep_lines <- (pr_subset_crispr>0.9) %>% colSums()
dep_lines
print(paste('Total Cls', nrow(pr_subset_crispr)))
