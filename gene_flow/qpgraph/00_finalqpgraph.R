library(rlang,lib.loc="/data/modules/R/4.5.2/lib64/R/library/")
library(tzdb,lib.loc="/data/modules/R/4.5.2/lib64/R/library/")
library(igraph,lib.loc="/data/modules/R/4.5.2/lib64/R/lib2/")
library(admixtools,lib.loc="/data/modules/R/4.5.2/lib64/R/library/")
library(admixr)
library(tibble)

#Compute F2 first
prefix = 'all'
my_f2_dir = paste0("/home/zajac/SYRI/QPGRAPH/F2")
extract_f2("combined", "F2/")

#Load a scaffold model
yang_graph = read.delim("SYRI/QPGRAPH/tree.scaffold.Yang")
#Check if model is valid structure
is_valid(edges_to_igraph(yang_graph))

#Test graph rearrangements and choose the one with lowest score
setwd("SYRI/QPGRAPH/thinned/")
library(readr)
read_table2 <- read_table
f2_blocks = f2_from_precomp("F2/")
newgraphs = graph_flipadmix(yang_graph)
newgraphs = newgraphs %>%
  rowwise %>%
  mutate(res = list(qpgraph(f2_blocks, graph))) %>%
  unnest_wider(res) %>%
  arrange(score)
winner = newgraphs[newgraphs$score == min(newgraphs$score),]$edges[[1]]
newgraphs = graph_minusone(edges_to_igraph(winner))
newgraphs = newgraphs %>%
  rowwise %>%
  mutate(res = list(qpgraph(f2_blocks, graph))) %>%
  unnest_wider(res) %>%
  arrange(score)
newgraphs = graph_plusone(edges_to_igraph(winner))
newgraphs = newgraphs %>%
  rowwise %>%
  mutate(res = list(qpgraph(f2_blocks, graph))) %>%
  unnest_wider(res) %>%
  arrange(score)
winner2 = newgraphs[newgraphs$score == min(newgraphs$score),]$edges[[1]]
newgraphs = graph_minusone(edges_to_igraph(winner2))
newgraphs = newgraphs %>%
  rowwise %>%
  mutate(res = list(qpgraph(f2_blocks, graph))) %>%
  unnest_wider(res) %>%
  arrange(score)
newgraphs = graph_flipadmix(edges_to_igraph(winner2))
newgraphs = newgraphs %>%
  rowwise %>%
  mutate(res = list(qpgraph(f2_blocks, graph))) %>%
  unnest_wider(res) %>%
  arrange(score)
newgraphs = graph_plusone(edges_to_igraph(winner2))
newgraphs = newgraphs %>%
  rowwise %>%
  mutate(res = list(qpgraph(f2_blocks, graph))) %>%
  unnest_wider(res) %>%
  arrange(score)