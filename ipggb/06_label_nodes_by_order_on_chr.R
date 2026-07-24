library(dplyr)
library(stringr)
library(readr)
library(data.table)
library(tidyverse)
library(ezRun)

args <- commandArgs(trailingOnly = TRUE)
n <- args[1]

base_dir <- "/home/zajac/pggb_allcommunities"
out_dir  <- "/home/zajac/PhyloTree/GenomicMosaicism"

files <- list.files(base_dir, pattern = n, recursive = TRUE)[grepl("paths", list.files(base_dir, pattern = n, recursive = TRUE))][-1]

chr_map <- fread("/home/zajac/chromosome_sets/community.chr.map")

for (i in files) {
  
  full_path <- file.path(base_dir, i)
  comm <- str_split(i, "/", simplify = TRUE)[1]
  comm_id <- str_remove(str_remove(comm, "community"), ".s5000")
  
  chr <- chr_map$chromosome[chr_map$community == comm_id]
  
  comm_dir <- file.path(base_dir, comm)
  
  # ---- Read and combine category files ----
  cat_df <- bind_rows(
    fread(list.files(comm_dir, pattern = "nodes.dispensable",
                     full.names = TRUE), header = FALSE) %>%
      mutate(cat = "dispensable"),
    fread(list.files(comm_dir, pattern = "nodes.private",
                     full.names = TRUE), header = FALSE) %>%
      mutate(cat = "private"),
    fread(list.files(comm_dir, pattern = "nodes.core",
                     full.names = TRUE), header = FALSE) %>%
      mutate(cat = "core")
  )
  
  colnames(cat_df)[1] <- "V1"
  
  # ---- Read length file ----
  len_df <- fread(list.files(comm_dir,
                             pattern = ".sorted.nodes.fasta.len",
                             full.names = TRUE),
                  header = FALSE)
  
  colnames(len_df)[1:2] <- c("V1", "V2")
  
  # ---- Read main file ----
  df <- fread(full_path, header = FALSE,
              col.names = c("V1", "orientation"))
  
  # ---- LEFT JOINS (order preserved) ----
  df <- df %>%
    left_join(cat_df, by = "V1") %>%
    left_join(len_df, by = "V1")
  
  # ---- Collapse consecutive categories ----
  df <- df %>%
    mutate(block = cumsum(cat != lag(cat, default = dplyr::first(cat)))) %>%
    group_by(block, cat) %>%
    summarise(
      V2 = sum(V2),
      segment_id = paste(V1, collapse = ","),
      .groups = "drop"
    ) %>%
    mutate(
      end   = cumsum(V2),
      start = end - V2 + 1,
      chr   = chr
    )
  
  # ---- Write output ----
  write_delim(
    df,
    file.path(out_dir, paste0(basename(i), ".order.length.txt")),
    delim = "\t",
    quote = "none"
  )
}