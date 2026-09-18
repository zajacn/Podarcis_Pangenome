library(dplyr)
library(tidyverse)
library(ezRun)
library(readr)
library(stringr)
library(DT)
library(sys)
library(data.table)

all_invs_mod = read.delim("/home/zajac/INVERSIONS/inversions.syri.impg.filtered.txt", header = T, sep = "\t")
all_invs_mod = all_invs_mod[,c(1:4)]
mapc_files = list.files("/home/zajac/SYRI/", pattern = "depth$", recursive = TRUE, full.names = TRUE)
mapc_files = mapc_files[!grepl("LacAg", mapc_files)]
args <- commandArgs(trailingOnly = TRUE)
x=as.numeric(args[1])

mapc_all_invs = NULL
df = fread(mapc_files[x], header = F)
df$V2 = as.numeric(df$V2)
species = str_remove(str_remove(sapply(str_split(mapc_files[x], "/"), .subset, 8), ".fasta.sorted.filtered.bam.depth"), ".rev")
invs_to_check = all_invs_mod[all_invs_mod$species == species,]
if (nrow(invs_to_check) != 0){
  invs_qv = apply(invs_to_check, 1, function(x) df[df$V2 > as.numeric(x[["start"]]) & df$V2 < as.numeric(x[["end"]]),] %>% separate_rows(V3, sep = ","))
  invs_qv = lapply(invs_qv, function(x){mean(as.numeric(x[["V3"]]))})
  invs_to_check = cbind(invs_to_check, data.frame("mean_coverage" = t(data.frame(invs_qv))))
  mapc_all_invs = rbind(mapc_all_invs,invs_to_check)
}

write_delim(mapc_all_invs, paste0("/home/zajac/INVERSIONS/temp/impg", x, ".out.txt"), delim = "\t")