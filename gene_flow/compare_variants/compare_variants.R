library(Matrix)          
library(vegan)            
library(ape)               
library(dendextend)        
library(ggplot2)
library(patchwork)
library(viridis)
library(reshape2)
library(ComplexHeatmap)   
library(data.table)
library(vcfR)
library(stringr)
library(dplyr)
library(tidyverse)

sample_data = read.delim("/groups/mpistaff/Zajac/genomes/sample_data.tsv")
sample_data$clade=factor(sample_data$clade,levels=c("Balkan","Sicilian-Maltese","Western","Siculus","Muralis","Iberian"))

##Load thinned SNP dataset - SNPs thinned to 1 SNP every 100bp confirmed with bcftools
snps = list.files("INVERSIONS/", pattern = "combined.thinned.vcf.gz$", recursive = TRUE, full.names = TRUE)
all_snps = NULL
for (i in snps){
  snps = read.vcfR(i)
  snps = data.frame(snps@gt)
  snps = data.frame(snps)[2:15] 
  snps = snps %>% dplyr::mutate(across(colnames(snps)[1:14], ~ case_when(. == "0|0" ~ 0, . == "0" ~ 0, . == "1|1" ~ 1, . == "1" ~ 1,  . == "0|1" ~ 1, . == "1|0" ~ 1, . == ".|1" ~ 1, . == "0|."~ 0, . == ".|0" ~ 0, . == "1|." ~ 1 )))
  all_snps = rbind(all_snps, snps)
}
all_snps_sparse <- Matrix(t(all_snps), sparse = TRUE)     
rm(all_snps); gc()
sampled_rows <- sample(ncol(all_snps_sparse), 5000000)
all_snps_sparse_sampled = all_snps_sparse[,sampled_rows]
##Load SVs
svs=list.files("/home/zajac/SV_intersects/largeSVs/", pattern = "merged.DEFAULT", full.names = TRUE)[!grepl("annot",list.files("/home/zajac/SV_intersects/largeSVs/", pattern = "merged.DEFAULT", full.names = TRUE))]
all_svs = NULL
for (i in svs){
  comb = read.vcfR(i)
  comb = cbind(data.frame(comb@fix), data.frame(comb@gt))
  comb = comb %>% dplyr::mutate(across(10:23, ~ str_remove(., ":.*"))) 
  comb = comb %>% dplyr::mutate(across(colnames(comb)[10:23], ~ case_when(. == "0|0" ~ 0, . == "0" ~ 0, . == "1|1" ~ 1, . == "1" ~ 1,  . == "0|1" ~ 1, . == "1|0" ~ 1, . == ".|1" ~ 1, . == "0|."~ 0, . == ".|0" ~ 0, . == "1|." ~ 1 )))
  comb = comb[abs(nchar(comb$REF) - nchar(comb$ALT)) > 50 & abs(nchar(comb$REF) - nchar(comb$ALT)) < 10000,]
  all_svs= rbind(all_svs,comb[10:23])
}
all_svs_sparse <- Matrix(t(all_svs), sparse = TRUE)  
rownames(all_svs_sparse) = sapply(str_split(rownames(all_svs_sparse), "_"), .subset, 3)
rm(all_svs); gc()


##Load Inversions
invs = read.delim("INVERSIONS/inversions.clustered.0.9proc.mapq.mapc.filtered.txt")
shared_haplotypes = unique(invs[invs$is_singleton == "FALSE",] %>% 
                             dplyr::select(3,15,16,17)) %>% 
  separate(species_list, into = c("species", "haplotype"), sep = "#") %>% 
  pivot_wider(names_from = haplotype, values_from = mean_length)
shared_species = invs[invs$is_singleton == "FALSE",] %>% 
  dplyr::select(3,15,16,17) %>% 
  separate(species_list, into = c("species", "haplotype"), sep = "#") %>%
  dplyr::select(1,3,4,5) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = mean_length) %>% 
  replace(is.na(.),0) 
all_invs = shared_species[rowSums(shared_species[3:16] > 0) > 1,c(3:16)] %>% replace(. > 0, 1)
all_invs <- t(all_invs) 

cat("SNP matrix:", dim(all_snps_sparse), "\n")   
cat("SV matrix:",  dim(all_svs_sparse),  "\n")       
cat("INV matrix:", dim(all_invs), "\n") 

##Distance matrices (Jaccard, binary)
jaccard_sparse <- function(mat_sparse) {
  co_presence <- as.matrix(mat_sparse %*% t(mat_sparse))   
  totals <- diag(co_presence)                                
  n <- nrow(co_presence)
  d <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      union_ij <- totals[i] + totals[j] - co_presence[i, j]
      d[i, j] <- if (union_ij == 0) 0 else 1 - co_presence[i, j] / union_ij
    }
  }
  rownames(d) <- colnames(d) <- rownames(mat_sparse)
  as.dist(d)
}

d_snp <- jaccard_sparse(all_snps_sparse)
d_sv <- jaccard_sparse(all_svs_sparse)
d_inv <- vegdist(all_invs, method = "jaccard", binary = TRUE)

pcoa_snp <- cmdscale(d_snp, k = 2)
pcoa_sv  <- cmdscale(d_sv,  k = 2)
pcoa_inv <- cmdscale(d_inv, k = 2)

proc_sv  <- procrustes(pcoa_snp, pcoa_sv)
proc_inv <- procrustes(pcoa_snp, pcoa_inv)

proc_df <- rbind(
  data.frame(Species=rownames(pcoa_snp), x=pcoa_snp[,1], y=pcoa_snp[,2], Type="SNP"),
  data.frame(Species=rownames(proc_sv$Yrot), x=proc_sv$Yrot[,1], y=proc_sv$Yrot[,2], Type="SV"),
  data.frame(Species=rownames(proc_inv$Yrot), x=proc_inv$Yrot[,1], y=proc_inv$Yrot[,2], Type="INV")
)

pD <- ggplot(proc_df, aes(x, y, color=Type, group=Species)) +
  geom_line(aes(group=Species), color="grey70", alpha=0.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#440154FF","#21908CFF","#FDE725FF")) +
  theme_minimal(base_size=10) +
  labs(title="Procrustes-aligned ordination", x="PCo1", y="PCo2") + 
  geom_text(aes(label = Species))

pD$data$Species = sample_data[match(pD$data$Species, sample_data$MyName),]$Organism.Name


##Mantel test 
mantel_snp_sv  <- mantel(d_snp, d_sv,  permutations = 999)
mantel_snp_inv <- mantel(d_snp, d_inv, permutations = 999)
mantel_sv_inv  <- mantel(d_sv,  d_inv, permutations = 999)
mantel_results <- data.frame(
  Comparison = c("SNP vs SV", "SNP vs INV", "SV vs INV"),
  r = c(mantel_snp_sv$statistic, mantel_snp_inv$statistic, mantel_sv_inv$statistic),
  p = c(mantel_snp_sv$signif,    mantel_snp_inv$signif,    mantel_sv_inv$signif)
)
mantel_results$p_adj <- p.adjust(mantel_results$p, method = "fdr")
print(mantel_results)

##Tanglegrams
d_snp_euc = dist(all_snps_sparse_sampled, method = "binary")
d_inv_euc = dist(all_invs, method = "binary")
d_svs_euc = dist(all_svs_sparse, method = "binary")
dend_snp <- as.dendrogram(hclust(d_snp_euc, "ward.D"))
dend_sv  <- as.dendrogram(hclust(d_svs_euc,  "ward.D"))
dend_inv <- as.dendrogram(hclust(d_inv_euc, "ward.D"))

dend_snp <- as.dendrogram(as.dist(d_snp) %>% hclust())
dend_sv  <-as.dendrogram(as.dist(d_sv) %>% hclust())
dend_inv <-as.dendrogram(as.dist(d_inv) %>% hclust())

tanglegram(dend_snp, dend_sv, main_left = "SNP", main_right = "SV",
           highlight_distinct_edges = TRUE, common_subtrees_color_lines = TRUE,
           lwd = 1.5, columns_width = c(4, 2, 4))
tanglegram(dend_snp, dend_inv, main_left = "SNP", main_right = "Inversion",
           highlight_distinct_edges = TRUE, common_subtrees_color_lines = TRUE,
           lwd = 1.5, columns_width = c(4, 2, 4))
tanglegram(dend_sv, dend_inv, main_left = "SNP", main_right = "Inversion",
           highlight_distinct_edges = TRUE, common_subtrees_color_lines = TRUE,
           lwd = 1.5, columns_width = c(4, 2, 4))
