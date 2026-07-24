#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(purrr)
  library(tidyverse)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 00_growth_curve.R <community_dir>")
community_dir <- args[1]
if (!dir.exists(community_dir)) stop(paste("Directory does not exist:", community_dir))
cat("Using community directory:", community_dir, "\n")

community_dir <- args[[1]]
out_file      <- paste0(community_dir, "/growth.curve.csv")
ncores        <- max(1, detectCores() - 1)
chunk_size    <- 200  # combos per batch

if (!dir.exists(community_dir)) stop("community_dir does not exist: ", community_dir)
cat("community_dir:", community_dir, "\n")
cat("out_file:", out_file, "\n")
cat("ncores:", ncores, "\n")
cat("chunk_size:", chunk_size, "\n\n")

genome_names_vec <- c(
  "PodMur1","rPodBoc1","rPodCre2.1","rPodErh1","rPodFil1","rPodGai1",
  "rPodLil1.2","rPodLio1","rPodMel1","rPodMur119","rPodPit1","rPodRaf1",
  "rPodSic1","rPodTil1","rPodVau1"
)

## Build all combinations (cap at 200 per combination size)
all_combs <- map(seq_along(genome_names_vec), function(genome_nb) {
  combos_list <- combn(genome_names_vec, genome_nb, simplify = FALSE)
  if (length(combos_list) > 200) combos_list <- sample(combos_list, 200)
  combos_list
})
all_combs_flat <- unlist(all_combs, recursive = FALSE)
total_combos <- length(all_combs_flat)
cat("Total combinations to process:", total_combos, "\n")

if (total_combos == 0) {
  stop("No combinations generated.")
}
## Find files and read genome tables
# gather full paths for each genome file (one per genome in genome_names_vec)

file_map <- setNames(
  vapply(genome_names_vec, function(x) {
    matches <- list.files(community_dir, pattern = paste0("^.*", x, "\\.paths$"), full.names = TRUE)
    if (length(matches) == 0) {
      stop("No .paths file found for genome: ", x)
    } else if (length(matches) > 1) {
      warning("Multiple .paths matches for genome ", x, "; using the first: ", matches[1])
      matches[1]
    } else {
      matches[1]
    }
  }, FUN.VALUE = character(1)),
  genome_names_vec
)

## Load node lengths file (same naming convention)
len_file <- file.path(community_dir, paste0(str_replace(basename(community_dir), ".s5000", ".sorted"), ".nodes.fasta.len"))
if (!file.exists(len_file)) stop("Length file not found: ", len_file)
length_table <- fread(len_file, header = FALSE)
setnames(length_table, c("V1", "V2"))

## Function to process genome combinations
process_combo_internal <- function(genomes, file_map_local, length_table_local) {
  nb_genomes <- length(genomes)
  dfs <- vector("list", length = nb_genomes)
  for (i in seq_len(nb_genomes)) {
    g <- genomes[[i]]
    fp <- file_map_local[[g]]
    dfs[[i]] <- fread(fp, header = FALSE)
  }
  combined <- rbindlist(dfs, use.names = FALSE)[ , .N, by = V1]  # .N = count per V1
  combined <- merge(combined, length_table_local, by = "V1", all.x = TRUE) # adds V2
  combined[ , class := fifelse(N == nb_genomes, "core",
                               fifelse(N == 1, "private", "dispensable"))]
  out <- combined[ , .(value = sum(V2, na.rm = TRUE)), by = class]
  out[ , nb_genomes := nb_genomes]
  setcolorder(out, c("nb_genomes", "class", "value"))
  return(out)
}

cl <- makeCluster(ncores)
on.exit({
  try(stopCluster(cl), silent = TRUE)
}, add = TRUE)
clusterExport(cl, varlist = c("file_map", "length_table", "process_combo_internal"), envir = environment())
clusterEvalQ(cl, {
  library(data.table)
})

## Streaming processing in chunks with parLapplyLB
indices <- seq_len(total_combos)
chunk_starts <- seq(1, total_combos, by = chunk_size)
cat("Processing in", length(chunk_starts), "chunks\n")
pb <- txtProgressBar(min = 0, max = length(chunk_starts), style = 3)
for (ci in seq_along(chunk_starts)) {
  start_idx <- chunk_starts[ci]
  end_idx   <- min(total_combos, start_idx + chunk_size - 1)
  chunk_idx <- indices[start_idx:end_idx]
  
  combos_chunk <- all_combs_flat[chunk_idx]
  results_list <- parLapplyLB(cl, combos_chunk, function(genomes) {
    process_combo_internal(genomes, file_map, length_table)
  })
  chunk_dt <- rbindlist(results_list, fill = TRUE)
  setcolorder(chunk_dt, c("nb_genomes", "class", "value"))
  fwrite(chunk_dt, out_file, append = file.exists(out_file))
  setTxtProgressBar(pb, ci)
}

close(pb)
cat("\nAll done. Results written to:", out_file, "\n")
