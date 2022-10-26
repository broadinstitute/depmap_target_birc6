library(igraph) 
library(network) 
library(sna)
library(ggraph)
library(visNetwork)
library(threejs)
library(networkD3)
library(tidyverse)
library(ggplot2)


clusters_to_highlight <- c(10,17,20,25,60)
novel_cluster <- c(17,20)
save <- F


graph_df <- list()
vertex_membership_list <- list()
for (num in 1:178){
  file_name = paste0("GLS_analysis/data/cluster_level/cluster_",num,".csv")
  
  graph_data <- read_csv(file_name) %>% 
    mutate(width = -log10(p_value)) %>%
    mutate(weight=p_value) %>%
    select(-c(`...1`)) %>%
    mutate(cluster_num=num)
  
  vertex_membership_list[[num]] <- data.frame(name=unique(graph_data$gene), cluster=num)
  graph_df[[num]] <- graph_data
}

combined_data <- do.call(rbind.data.frame, graph_df) %>% 
  filter(cluster_num %in% clusters_to_highlight) 

vertex_membership <- do.call(rbind.data.frame, vertex_membership_list) %>% 
  filter(cluster %in% clusters_to_highlight) %>% 
  mutate(node_col=ifelse(cluster %in% novel_cluster,'red','blue'))


net <- igraph::graph_from_data_frame(combined_data, vertex_membership, directed=F)
net <- igraph::simplify(net, remove.multiple = T, remove.loops = T)

V(net)$size <- igraph::degree(net, mode="in")
E(net)$width <- -log10(E(net)$weight)

vertex_info <- data.frame('size'=V(net)$size, 'cluster_num'=V(net)$cluster, 'name'=V(net)$name) %>% 
  group_by(cluster_num) %>%
  mutate(radius=size*size/sum(size*size))


stopifnot(all(vertex_info$name==V(net)$name))

V(net)$radius <- vertex_info$radius


ggraph(net, layout = 'fr') +
  geom_node_point(aes(size=radius, color=node_col),alpha=0.7, show.legend = F)+
  geom_edge_link(aes(width = width), alpha=0.5, color='grey50') +   # add edges to the plot
  geom_node_text(aes(label=name), size=5) +
  facet_nodes(~cluster)+
  scale_size(range = c(1,15))+
  scale_edge_width(range = c(0.25, 2.5))+
  guides(edge_width=guide_legend("-log10(p value)"))+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

if (save){
  ggsave(paste0("output/cluster_plots/cluster_",num,".pdf"))
}
