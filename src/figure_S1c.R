library(tidyverse)
library(magrittr)
library(here)
library(ggplot2)

source(here::here('src', 'helper_functions.R'))

genes <- c('UBA6', 'BIRC6', 'KCMF1', 'UBR4')

## Achilles

demeter2_RNAi_ge <- data.table::fread(here::here('data', 'demeter2-achilles_v12-gene-effect.csv')) %>% 
  column_to_rownames('V1') %>% 
  extract_hugo_symbol_colnames() 



ge_subset_RNAi <- demeter2_RNAi_ge[, genes] %>% rownames_to_column(var='CLs') %>% 
  pivot_longer(-CLs, names_to='Gene', values_to='gene_effect')

colors <- c(rgb(195/255,223/255,233/255), 
            rgb(137/255,192/255,212/255),
            rgb(80/255,161/255,192/255),
            rgb(18/255,132/255,172/255))

ge_subset_RNAi$Gene <- factor(ge_subset_RNAi$Gene, levels = genes)

ggplot(ge_subset_RNAi, aes(x=gene_effect))+
  geom_density(aes(fill = Gene), show.legend = F)+
  facet_wrap(~Gene, nrow=4, ncol=1)+
  geom_vline(xintercept = -1, col='red')+
  scale_fill_manual(values=colors)+
  theme_minimal()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab('Gene Effect (DEMETER2)')

demeter2_RNAi_pr <- data.table::fread(here::here('data', 'demeter2-achilles_v12-gene-dependency.csv')) %>% 
  column_to_rownames('V1') %>% 
  extract_hugo_symbol_colnames() 

pr_subset_RNAi <- demeter2_RNAi_pr[, genes]

sink(here::here('output', 'FigS1', 'S1c_CL_nums.txt'))

print('Strongly dependent CLs')
dep_lines <- (pr_subset_crispr>0.9) %>% colSums()
dep_lines
print(paste('Total Cls', nrow(pr_subset_crispr)))
