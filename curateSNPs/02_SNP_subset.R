library(data.table)
library(dplyr)
library(tidyverse)
library(vcfR)

args <- commandArgs(trailingOnly = TRUE)
n=args[1]

#Run as Rsctipt 02_SNP_subset.R <community_number>
#Subset to SNPs present in both methods
confirmed_subset=NULL
pggbdir=paste0("/home/zajac/pggb_allcommunities/community", n ,".s5000/")
bcfdir=paste0("/home/zajac/SYRI/community", n, "/community", n, ".rPodCre2.1_to_all/")
for (i in str_remove(str_remove(list.files(pggbdir, pattern = "SNP.vcf"), paste0("community", n, ".")), ".pggb.waved.norm.SNP.vcf")[!grepl("annot", str_remove(str_remove(list.files(pggbdir, pattern = "SNP.vcf"), paste0("community", n, ".")), ".pggb.waved.norm.SNP.vcf"))]){
  pggb=list.files(pggbdir, pattern = "SNP.vcf")
  pggb = pggb[grepl(i,pggb)][1]
  bcf = list.files(bcfdir, pattern = "snps_only.vcf")
  bcf = bcf[grepl(i, bcf)][1]
  pggb = read.vcfR(paste0(pggbdir, pggb))
  bcf = read.vcfR(paste0(bcfdir, bcf))
  pggb = cbind(data.frame(pggb@fix), data.frame(pggb@gt))
  bcf = cbind(data.frame(bcf@fix), data.frame(bcf@gt))
  res = intersect(pggb$POS, bcf$POS)
  pggb = pggb[pggb$POS %in% res,]
  confirmed_subset[[i]] = pggb
}

confirmed_subset_positions = lapply(confirmed_subset, function(x){x$POS})
values_to_keep = data.frame(values = unlist(unname(confirmed_subset_positions))) %>% group_by(values) %>% summarise(count = n()) %>% filter(count >= 2) %>% pull(values)
values_to_keep = sort(as.numeric(values_to_keep))
confirmed_subset_temp = lapply(confirmed_subset, function(x) {x %>% filter(POS %in% values_to_keep)})

combined = NULL
for (i in names(confirmed_subset_temp)){
  if (is.null(combined)){
    combined = confirmed_subset_temp[[i]][,c(-3,-8)] 
  } else {
    combined = merge(combined, confirmed_subset_temp[[i]][,c(-3,-8)], by = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER" ,"FORMAT"), all = TRUE)
  }
}

#Keep SNPS separated by at least 100bp
thin_chr <- function(df, min_dist = 50) {
  keep <- logical(nrow(df))
  last_pos <- -Inf
  
  for (i in seq_len(nrow(df))) {
    if (df$POS[i] - last_pos >= min_dist) {
      keep[i] <- TRUE
      last_pos <- df$POS[i]
    }
  }
  
  df[keep, ]
}


thinned = combined %>% mutate(POS = as.numeric(POS)) %>% arrange(CHROM, POS) %>%
  group_by(CHROM) %>% 
  group_modify(~ thin_chr(.x, 50)) %>%
  ungroup() %>% data.frame()

thinned$INFO = "LEN=1;TYPE=snp"
thinned$ID = seq_along(thinned$POS)
thinned = thinned[,c(1,2,23,3,4,5,6,22,7,8:21)]
thinned = thinned %>% mutate(across(10:23, ~ ifelse(is.na(.x), 0, .x)))

#Convert to VCF
meta <- c('##fileformat=VCFv4.2',
          '##source=R_thinned_data',
          paste0('##fileDate=', Sys.Date()),
          '##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">',
          '##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">',
          '##INFO=<ID=ORIGIN,Number=1,Type=String,Description="Decomposed from a complex record using vcflib vcfwave and alignment with WFA2-lib.">',
          '##INFO=<ID=INV,Number=0,Type=Flag,Description="Inversion detected">',
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

vcf <- new("vcfR",
             meta = meta,
             fix = as.matrix(thinned[, 1:8]),
             gt = as.matrix(thinned[, 9:23, drop = FALSE]))
#Save VCF
filename=paste0("/home/zajac/INVERSIONS/community", n, "/combined.thinned50.vcf.gz")
write.vcf(vcf, file = filename)
