#!/usr/bin/R
library(stringr)
library(DT)
library(sys)
library(data.table)
library(dplyr)
library(tidyverse)
library(ezRun)

# Computes dispensable node sharing between all paths in a pangenome variation
# graph (GFA format) in sliding windows of Xbp.
# 
# Dispensable nodes = nodes present in more than one but not all genomes.
# 
# Two modes:
#   Default   : counts shared dispensable bp regardless of node order
#   colinear: finds co-linear blocks globally across the full path,
#               then bins the co-linear bp into windows.
#               Uses sparse DP - only processes focal positions that
#               match each comp node, skipping the vast majority of
#               work for typical pangenome dispensable node distributions.
# 
# Usage:
#   Rscript ~/pggb_allcommunities/06_1Mb_dispnodesharing_perc.diagnostic.R <gfa_file> <colinear> <window_size> <outfile>
# This script does not scale with genome sizes and so serves as a diagnostic/sanity check for the python script pangenome_dispensable_windows.py.


args = commandArgs(trailingOnly=TRUE)
file=args[1]
colinear=args[2]
window=as.numeric(args[3])
out=args[4]

gfa = readLines(file)
paths = length(gfa[grepl("^P", gfa)])
rest = length(gfa) - paths
node_lengths = as.data.frame(do.call(rbind,(str_split(gfa[1:rest], "\t")[sapply(str_split(gfa[1:rest], "\t"), function(z) "S" %in% z)]))) %>% mutate(V3 = nchar(V3))
node_lengths = node_lengths[,c(1:3)]
colnames(node_lengths) = c("S","X3", "Length")

lst = NULL
for (i in 1:paths){
  df = gfa[grepl("^P", gfa)][i]
  ddf = data.frame(values = str_split(df, "\t")[[1]]) 
  ddf = data.frame(t(tibble(ddf)))[,-4] %>% separate_rows(X3, sep = ",", convert = TRUE)
  ddf = ddf %>% rownames_to_column()
  lst[[i]] = ddf
}

lst = lapply(lst, function(x){x$X3 = gsub("[+-]", "", x$X3); return(x)})
lst = lapply(lst, function(x){x %>% left_join(node_lengths, by = "X3")})
lst = lapply(lst, function(x){x$cumulative = cumsum(x$Length); return(x)})
path_names   = sapply(lst, function(x) unique(x$X2))
genome_names = sapply(path_names, function(x) strsplit(x, "#", fixed = TRUE)[[1]][1])
n_genomes = length(unique(genome_names))
genome_to_paths = split(seq_along(lst), genome_names)
genome_node_sets = lapply(genome_to_paths, function(idx){
  unique(unlist(lapply(lst[idx], function(x) unique(x$X3))))
})
all_genome_nodes = unlist(genome_node_sets)

path_lengths = sapply(lst, function(x) max(x$cumulative))
names(path_lengths) = sapply(lst, function(x) unique(x$X2))

counts = data.frame(table(all_genome_nodes), stringsAsFactors = FALSE)
counts$all_genome_nodes = as.character(counts$all_genome_nodes)

core        = counts[counts$Freq == n_genomes, ]$all_genome_nodes
private     = counts[counts$Freq == 1,          ]$all_genome_nodes
dispensable = counts[counts$Freq > 1 & counts$Freq < n_genomes, ]$all_genome_nodes

lst2 = lapply(lst, function(x){x = x[x$X3 %in% dispensable,]; return(x)})
lst2 = lapply(seq(1,length(lst2),1), function(x){
  lst2[[x]]$start = lst2[[x]]$cumulative - lst2[[x]]$Length
  lst2[[x]]$end = lst2[[x]]$cumulative
  return(lst2[[x]])
})

if (colinear == "colinear"){

  lst4 = lapply(lst2, function(x){x %>%
    mutate(
      start = as.integer(start),end   = as.integer(end),X3    = as.integer(X3)
    ) %>%
    arrange(start, end) %>%
    mutate(consecutive = lag(end) == start) %>%
    mutate(consecutive = ifelse(is.na(consecutive), FALSE, consecutive),
           run_id = cumsum(!consecutive)) %>%
    group_by(run_id) %>%
    mutate(run_length = n()) %>%
    ungroup() %>%
    filter(run_length >= 2) %>% arrange(start)})
  
  combs = data.frame(t(combn(seq(1,length(lst),1), 2)))
  combs = rbind(combs, combs[,c(2,1)] %>% dplyr::rename(X2 = X1, X1 = X2))
  combs = combs[genome_names[combs$X1] != genome_names[combs$X2], ]
  
  names(lst4) = lapply(lst4, function(x){unique(x$X2)})
  lst5 = apply(combs, 1, function(x){
    n1 = x[1]
    n1 = names(lst4)[n1]
    df1 = lst4[[n1]]
    genome_length = path_lengths[[n1]]
    n2 = x[2]
    n2 = names(lst4)[n2]
    df2 = lst4[[n2]]
    df1$comp = paste0(n1,"_", n2)
    dfx = merge(df1, df2, by = "X3") %>% 
      arrange(start.x) %>%
      mutate(
        consecutive_x = lag(end.x) == start.x,
        consecutive_y = lag(end.y) == start.y,
        consecutive   = consecutive_x & consecutive_y
      ) %>%
      
      mutate(
        consecutive = ifelse(is.na(consecutive), FALSE, consecutive),
        run_id = cumsum(!consecutive)
      ) %>% 
      group_by(run_id) %>% 
      mutate(run_length = n()) %>%
      filter(run_length >= 2) %>% 
      ungroup() %>% 
      mutate(blocks = ezCut(end.x, seq(1,genome_length,window))) %>% 
      group_by(blocks, comp) %>% 
      summarise(sum(Length.x)) 
    return(dfx)
  })
  
  write_delim(bind_rows(lst5), "outr.csv", delim = "\t")
} else {
  combs = data.frame(t(combn(seq(1,length(lst),1), 2)))
  combs = rbind(combs, combs[,c(2,1)] %>% dplyr::rename(X2 = X1, X1 = X2))
  combs = combs[genome_names[combs$X1] != genome_names[combs$X2], ]
  names(lst2) = lapply(lst2, function(x){unique(x$X2)})
  lst5 = apply(combs, 1, function(x){
    n1 = x[1]
    n1 = names(lst2)[n1]
    df1 = lst2[[n1]]
    genome_length = path_lengths[[n1]]
    n2 = x[2]
    n2 = names(lst2)[n2]
    df2 = lst2[[n2]]
    df1$comp = paste0(n1,"_", n2)
    dfx = merge(df1, df2, by = "X3") %>% 
      arrange(start.x) %>% 
      mutate(blocks = ezCut(end.x, seq(1,genome_length,window))) %>% 
      group_by(blocks, comp) %>% 
      summarise(sum(Length.x)) 
    return(dfx)
  })
  write_delim(bind_rows(lst5), "outr.csv", delim = "\t")
}
