#!/usr/bin/env Rscript
library(vcfR)
library(utils)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ezRun)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
n <- args[[1]]
all = NULL
pggbdir=paste0("/home/zajac/pggb_allcommunities/community", n ,".s5000/")
for (i in str_remove(str_remove(list.files(pggbdir, pattern = "SNP.vcf"), paste0("community", n, ".")), ".pggb.waved.norm.SNP.vcf")[!grepl("annot", str_remove(str_remove(list.files(pggbdir, pattern = "SNP.vcf"), paste0("community", n, ".")), ".pggb.waved.norm.SNP.vcf"))]){
  pggb=list.files(pggbdir, pattern = "SNP.vcf")
  pggb = pggb[grepl(i,pggb)][1]
  pggb = read.vcfR(paste0(pggbdir, pggb))
  pggb = cbind(data.frame(pggb@fix), data.frame(pggb@gt))
  all[[i]] = pggb
}
all = lapply(all, function(x){x[,c(1,2,4,5,10)]})
merged_df <- Reduce(function(x, y) merge(x, y, by = c("CHROM","POS","REF","ALT"), all = TRUE), all)
merged_df$N = rowSums(!is.na(merged_df[,c(5:18)]))
pos = read.delim(paste0("/home/zajac/SYRI/community", n , "/community", n, ".rPodCre2.1_to_all/positions.txt"), header = F)
print(paste("percent of supported singletons: ", length(intersect(merged_df[merged_df$N == 1,]$POS,  pos$V1))/length(merged_df[merged_df$N == 1,]$POS)))
print(paste("percent of supported nonsingletons: ",length(intersect(merged_df[merged_df$N > 1,]$POS,  pos$V1))/length(merged_df[merged_df$N > 1,]$POS)))
merged_df = merged_df[merged_df$POS %in% pos$V1,]
print(paste0("Number of SNPS: ", nrow(merged_df)))
write_delim(merged_df, paste0(pggbdir, "merged.supported.pggb.waved.norm.SNP.txt"), delim = "\t")