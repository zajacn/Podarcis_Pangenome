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
library(ape)
library(phytools)
library(ggrepel)

sample_data = read.delim("/groups/mpistaff/Zajac/genomes/sample_data.tsv")
sample_data$clade=factor(sample_data$clade,levels=c("Balkan","Sicilian-Maltese","Western","Siculus","Muralis","Iberian"))

##Load thinned SNP dataset - SNPs thinned to 1 SNP every 50bp confirmed with bcftools
snps = list.files("INVERSIONS/", pattern = "combined.thinned50.vcf.gz$", recursive = TRUE, full.names = TRUE)
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
d_snp <- vegdist(all_snps_sparse, method = "jaccard", binary = TRUE)
d_sv <- vegdist(all_svs_sparse, method = "jaccard", binary = TRUE)
d_inv <- vegdist(all_invs, method = "jaccard", binary = TRUE)

##Build trees consistently (NJ is preferable to hclust for cross-marker comparison)
tr_snp <- nj(d_snp)
tr_sv  <- nj(d_sv)
tr_inv <- nj(d_inv)

##Root NJ trees
tr_snp$node.label <- as.character((Ntip(tr_snp)+1):(Ntip(tr_snp)+Nnode(tr_snp)))
tr_snp_r = midpoint.root(tr_snp)
tr_snp_r <- root(tr_snp_r, node = 17, resolve.root = TRUE)
tr_sv$node.label <- as.character((Ntip(tr_sv)+1):(Ntip(tr_sv)+Nnode(tr_sv)))
tr_sv_r  <- midpoint.root(tr_sv)
tr_sv_r <- root(tr_sv_r, node = 17, resolve.root = TRUE)
tr_inv$node.label <- as.character((Ntip(tr_inv)+1):(Ntip(tr_inv)+Nnode(tr_inv)))
tr_inv_r <- midpoint.root(tr_inv)
tr_inv_r <- root(tr_inv_r, node = 17, resolve.root = TRUE)

##Check tip labels
tr_snp_r$tip.label = str_replace(sample_data[match(tr_snp_r$tip.label,sample_data$MyName),]$Organism.Name, "Podarcis ", "P. ")
tr_sv_r$tip.label = str_replace(sample_data[match(tr_sv_r$tip.label,sample_data$MyName),]$Organism.Name, "Podarcis ", "P. ")
tr_inv_r$tip.label = str_replace(sample_data[match(tr_inv_r$tip.label,sample_data$MyName),]$Organism.Name, "Podarcis ", "P. ")

##Matching the trees
## SNP vs SV
pdf("../Figures/SNP_SV_tree.pdf", width = 5, height = 5)
cophy_snp_sv <- cophylo(tr_snp_r, tr_sv_r, rotate = TRUE)  # rotate=TRUE minimizes crossings
plot(cophy_snp_sv,
     link.lwd = 2,
     link.type = "curved",
     link.col = "#21908CFF",
     fsize = 0.8)
title("SNP tree vs SV tree")
dev.off()

## SNP vs Inversions
pdf("../Figures/SNP_INV_tree.pdf", width = 5, height = 5)
cophy_snp_inv <- cophylo(tr_snp_r, tr_inv_r, rotate = TRUE)
plot(cophy_snp_inv,
     link.lwd = 2,
     link.type = "curved",
     link.col = "darkorange",
     fsize = 0.8)
title("SNP tree vs Inversion tree")
dev.off()

## SV vs Inversions
pdf("../Figures/SV_INV_tree.pdf", width = 5, height = 5)
cophy_sv_inv <- cophylo(tr_sv_r, tr_inv_r, rotate = TRUE)
plot(cophy_sv_inv,
     link.lwd = 2,
     link.type = "curved",
     link.col = "#440154FF",
     fsize = 0.8)
title("SV tree vs Inversion tree")
dev.off()

##Procrustes on ordinations
stopifnot(all(rownames(pcoa_snp) %in% rownames(pcoa_sv)))
stopifnot(all(rownames(pcoa_snp) %in% rownames(pcoa_inv)))

common_order <- rownames(pcoa_snp)          # use SNP tree order as the reference
pcoa_snp <- pcoa_snp[common_order, ]
pcoa_sv  <- pcoa_sv[common_order, ]
pcoa_inv <- pcoa_inv[common_order, ]

proc_snp_sv  <- protest(pcoa_snp, pcoa_sv,  permutations = 999)
proc_snp_inv <- protest(pcoa_snp, pcoa_inv, permutations = 999)
rownames(proc_snp_sv$Yrot)  <- common_order
rownames(proc_snp_inv$Yrot) <- common_order

proc_df <- rbind(
  data.frame(Species = common_order,
             x = pcoa_snp[, 1], y = pcoa_snp[, 2], Type = "SNP"),
  data.frame(Species = common_order,
             x = proc_snp_sv$Yrot[, 1], y = proc_snp_sv$Yrot[, 2], Type = "SV"),
  data.frame(Species = common_order,
             x = proc_snp_inv$Yrot[, 1], y = proc_snp_inv$Yrot[, 2], Type = "INV")
)
proc_df$Species <- sample_data$Organism.Name[match(proc_df$Species, sample_data$MyName)]

pD <- ggplot(proc_df, aes(x, y, color = Type, group = Species)) +
  geom_line(color = "grey70", alpha = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "darkorange")) +
  theme_bw(base_size = 10) +
  labs(title = "Procrustes-aligned ordination", x = "PCo1", y = "PCo2") +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = 20, show.legend = FALSE)

pdf("../Figures/SV_INV_SV_procrustes.pdf", height = 7, width = 11)
pD
dev.off()

##Mantel test visualised in a heatmap
mantel_snp_sv  <- mantel(d_snp, d_sv,  permutations = 999)
mantel_snp_inv <- mantel(d_snp, d_inv, permutations = 999)
mantel_sv_inv  <- mantel(d_sv,  d_inv, permutations = 999)
mantel_results <- data.frame(
  Comparison = c("SNP vs SV", "SNP vs INV", "SV vs INV"),
  r = c(mantel_snp_sv$statistic, mantel_snp_inv$statistic, mantel_sv_inv$statistic),
  p = c(mantel_snp_sv$signif,    mantel_snp_inv$signif,    mantel_sv_inv$signif)
)
mantel_results$p_adj <- p.adjust(mantel_results$p, method = "fdr")

congruence_mat <- matrix(
  c(NA, mantel_snp_sv$statistic, mantel_snp_inv$statistic,
    mantel_snp_sv$statistic, NA, mantel_sv_inv$statistic,
    mantel_snp_inv$statistic, mantel_sv_inv$statistic, NA),
  nrow = 3, dimnames = list(c("SNP","SV","INV"), c("SNP","SV","INV"))
)

pdf("../Figures/SV_INV_SV_mantel.pdf", height = 5, width = 7)
pheatmap(congruence_mat, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, color = viridis(50), na_col = "white",
         main = "Mantel r between marker-type distance matrices")
dev.off()

##Objectively choose linkage via cophenetic correlation
methods <- c("single", "complete", "average", "ward.D2")

pick_method <- function(d) {
  cors <- sapply(methods, function(m) cor(d, cophenetic(hclust(d, method = m))))
  names(sort(cors, decreasing = TRUE))[1]
}

best_snp <- pick_method(d_snp)
best_sv  <- pick_method(d_sv)
best_inv <- pick_method(d_inv)
best_snp; best_sv; best_inv   #Chek

##Build tanglegrams
dend_snp <- as.dendrogram(hclust(d_snp, method = best_snp))
dend_sv  <- as.dendrogram(hclust(d_sv,  method = best_sv))
dend_inv <- as.dendrogram(hclust(d_inv, method = best_inv))

tanglegram(dend_snp, dend_sv,  main_left = "SNP", main_right = "SV")
tanglegram(dend_snp, dend_inv, main_left = "SNP", main_right = "INV")
