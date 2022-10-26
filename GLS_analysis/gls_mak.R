library(readr)
library(here)
library(tidyverse)

GLS = function(X){
  require(corpcor)
  require(reshape2)
  require(dplyr)
  require(Rcpp)
  # for fast matrix multiplication
  sourceCpp(here("matmul.cpp"))
  
  D = corpcor::invcov.shrink(t(X), 0, 0)
  X = scale(X, scale =F)
  X.D = eigenMapMatMult(t(X),D)
  
  num = eigenMapMatMult(X.D,X);   den = diag(num)
  colnames(num) = colnames(X); rownames(num) = colnames(X)
  
  beta = num / den; 
  beta.se = sqrt(den) / den
  p = 2*pnorm(-abs(beta / beta.se)); diag(p) = NA
  q = apply(p, 1, function(x) p.adjust(x, method = "BH"))
  
  return(dplyr::filter(cbind(reshape2::melt(q, varnames = c("y", "x"), value.name = "q"),
                             p = c(p), beta = c(beta), beta.se = c(beta.se)), is.finite(q)))
}





gene.effect <- read_csv(here('data', 'CRISPR_gene_effect.csv')) %>% column_to_rownames('DepMap_ID') %>% as.matrix()
gene.effect <- gene.effect[rowSums(is.na(gene.effect))==0,]
X = gene.effect[apply(gene.effect, 1, function(x) all(is.finite(x))),]

res = GLS(X)

save(res, file = 'data/GLS_chronos.rData')
load('data/GLS_chronos.rData')

library(tidyverse)
library(magrittr)

res_wider <- res %>% 
  dplyr::select(c('x','y','p')) %>%
  pivot_wider(names_from=y,values_from=p)

rm(res)

res_wider<- res_wider %>%  column_to_rownames('x')

res_wider %<>%
  as.matrix()

res_wider <- res_wider[,rownames(res_wider)]

### Setting the diag pvalues to 1 as we are not interested in the association of a gene to itself
diag(res_wider) <- 1

rownames(res_wider) <- str_split(rownames(res_wider), pattern=' ', simplify=T)[,1]
colnames(res_wider) <- str_split(colnames(res_wider), pattern=' ', simplify=T)[,1]


library(feather)

res_wider<- data.frame(res_wider)

write_feather(res_wider, 'GLS_p.ftr')
