library(tidyverse)
library(cdsr)
library(ggplot2)
library(knitr)
library(taigr)
library(useful)
library(arrow)
library(RcppCNPy)
library(magrittr)

library(ggplot2)


Achilles.gene.effect <- data.table::fread('GLS_analysis/data/CRISPR_gene_effect.csv')
Achilles.gene.effect %<>% column_to_rownames('DepMap_ID')
colnames(Achilles.gene.effect) <- str_split_fixed(colnames(Achilles.gene.effect), ' ',2)[,1]
colnames(Achilles.gene.effect) %<>% str_replace_all('-', '.')

Achilles.gene.dependency <- data.table::fread('GLS_analysis/data/CRISPR_gene_dependency.csv')
Achilles.gene.dependency %<>% column_to_rownames('DepMap_ID')
colnames(Achilles.gene.dependency) <- str_split_fixed(colnames(Achilles.gene.dependency), ' ',2)[,1]
colnames(Achilles.gene.dependency) %<>% str_replace_all('-', '.')

cluster_info <- read_csv('GLS_analysis/data/cluster_2000.csv')



gls_pvalues <-  arrow::read_feather('GLS_analysis/data/GLS_p.ftr') %>% 
  as.matrix()

rownames(gls_pvalues) <- colnames(gls_pvalues)

mask <- arrow::read_feather('GLS_analysis/data/binary_matrix_2000.ftr') %>% 
  column_to_rownames('Gene') %>% 
  as.matrix()

gene_names <- colnames(gls_pvalues)

assertthat::are_equal(rownames(gls_pvalues), rownames(mask))
assertthat::are_equal(colnames(gls_pvalues), colnames(mask))
assertthat::are_equal(colnames(gls_pvalues), colnames(Achilles.gene.effect))


dependent_lines <- Achilles.gene.dependency>0.5

i<-1
row_data <- list()

for (cluster in cluster_info$community_string){
  
  gene_list <- str_split(cluster, ',') %>% unlist()
  
  subset_gene_effect <- Achilles.gene.effect[,gene_list]
  subset_gene_dependency <- dependent_lines[,gene_list]
  
  gene_effect_dependent = subset_gene_effect * subset_gene_dependency
  dependency_line_gene_effects <- list()
  
  for (j in 1:length(gene_list)){
    gene_level <- gene_effect_dependent[,j]
    
    dependency_line_gene_effects[j] <- gene_level[gene_level<0 & !is.na(gene_level)] %>% mean(na.rm=T)
  }
  dependency_line_gene_effects <- setNames(unlist(dependency_line_gene_effects), gene_list)
  
  gls_pvalues_subset <- gls_pvalues[gene_list,gene_list]
  p_values <- sort(gls_pvalues_subset[upper.tri(gls_pvalues_subset)])[1:3]
  count_zeroes <- sum()
  hmp_top_3 <- 3/(sum(1/p_values[1]+1/p_values[2]+1/p_values[3]))
  
  
  mean_of_mean_CERES <- mean(rowMeans(subset_gene_effect, na.rm = T), na.rm = T)
  mean_mean_CERES_dependent_lines <- mean(dependency_line_gene_effects, na.rm = T)
  mean_of_mins <- mean(apply(subset_gene_effect, MARGIN = 2, FUN=min, na.rm = T))
  cluster_variance <- var(rowMeans(subset_gene_effect, na.rm = T), na.rm = T)
  
  row_data[[i]] <- c(mean_of_mean_CERES, mean_mean_CERES_dependent_lines, mean_of_mins, cluster_variance,hmp_top_3)
  i<-i+1
}

row_data <- as.data.frame(do.call(rbind, row_data))

colnames(row_data) <- c('mean_of_mean_CERES', 'mean_of_mean_CERES_dependent_lines', 'mean_of_mins', 'cluster_variance', 'hmp_top_3_interactions')


cluster_data <- cbind(cluster_info, row_data) 

write_csv(cluster_data %>% select(-c(`...1`)), 'GLS_analysis/data/cluster_file_with_metrics.csv')