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
bcfdir=paste0("/home/zajac/SYRI/community", n ,"/community", n, ".rPodCre2.1_to_all/")
for (i in str_remove(str_remove(list.files(pggbdir, pattern = "SNP.vcf"), paste0("community", n, ".")), ".pggb.waved.norm.SNP.vcf")[!grepl("annot", str_remove(str_remove(list.files(pggbdir, pattern = "SNP.vcf"), paste0("community", n, ".")), ".pggb.waved.norm.SNP.vcf"))]){
  pggb=list.files(pggbdir, pattern = "SNP.vcf")
  pggb = pggb[grepl(i,pggb)][1]
  pggb = read.vcfR(paste0(pggbdir, pggb))
  pggb = cbind(data.frame(pggb@fix), data.frame(pggb@gt))
  bcftools=list.files(bcfdir, pattern = "snps_only.vcf$")
  bcftools = bcftools[grepl(i,bcftools)][1]
  bcftools = read.vcfR(paste0(bcfdir, bcftools))
  bcftools = cbind(data.frame(bcftools@fix), data.frame(bcftools@gt))
  print(paste0("Number of confirmed positions for ", i, " ", length(intersect(pggb$POS, bcftools$POS))*100/ nrow(pggb)))
  pggb =  pggb[pggb$POS %in% bcftools$POS,]
  all[[i]] = pggb
}
all = lapply(all, function(x){x[,c(1,2,4,5,10)]})
merged_df <- Reduce(function(x, y) merge(x, y, by = c("CHROM","POS","REF","ALT"), all = TRUE), all)
merged_df$N = rowSums(!is.na(merged_df[,c(5:18)]))
print(paste0("Number of SNPS: ", nrow(merged_df)))
write_delim(merged_df, paste0(pggbdir, "merged.supported.pggb.waved.norm.SNP.txt"), delim = "\t")
