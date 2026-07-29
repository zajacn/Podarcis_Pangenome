## This script analyses colinear disp node sharing in Podarcis raffonei

library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)

##Load inversions
invs = read.delim("../INVERSIONS/inversions.clustered.0.9proc.mapq.mapc.filtered.node_annotation.txt")
shared = invs[invs$is_singleton == "FALSE",] %>% 
  dplyr::select("species_list","max_start","min_end","mean_length") %>% 
  separate(species_list, into = c("species", "haplotype"), sep = "#") %>%
  dplyr::select(1,3,4,5) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = mean_length) %>% 
  replace(is.na(.),0) %>% mutate(across(c(3:16), ~ if_else(. > 0, 1,0))) %>%
  mutate(N = rowSums(.[3:16])) 

invs$shared_between_species = ""
invs[paste0(invs$max_start, "-", invs$min_end) %in% paste0(shared[shared$N >1,]$max_start,"-",shared[shared$N >1,]$min_end),]$shared_between_species = "Yes"

##Load bed files for all species
beds = list.files("/groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_unmasked/", pattern = "fai", full.names = TRUE)
beds = lapply(beds, function(x){ 
  bed = read.delim(x, header = F)
  bed$Start = 1
  bed = bed[,c(1,6,2)]
  colnames(bed) = c("Chr", "Start", "End")
  bed  = bed %>% separate(Chr, into = c("Focal_genome",  "Focal_hap", "Focal_chr"), sep = "#")
  return(bed)
})
beds = bind_rows(beds)

##Load in the colinear disp node sharing across chromosomes
all_noncolinear = NULL
for (i in seq(0,18,1)){
  df = read.delim(paste0("../PhyloTree/GenomicMosaicism/Disp.node.sharing.combined.50kb.community", i, ".colinear.csv"), sep = ",")
  df  = df %>% 
    separate(Compared_Path, into = c("Compared_genome", "Compared_hap", "Compared_chr"), sep = "#") %>% 
    separate(Focal_Path, into = c("Focal_genome", "Focal_hap", "Focal_chr"), sep = "#")
  df = df %>% 
    group_by(Focal_genome, Focal_hap, Focal_chr,Chromosome_Window_Start, Chromosome_Window_End) %>% 
    slice_max(as.numeric(Shared_pct), n = 1, with_ties = TRUE) %>% 
    unique() %>% 
    left_join(beds, by = c("Focal_genome", "Focal_hap", "Focal_chr"))
  df = merge(df, sample_data[,c(3,16,17)], by.x = "Focal_genome", by.y = "MyName", all.x = TRUE) %>% unique()
  df = merge(df, sample_data[,c(16,17)] %>% dplyr::rename("Compared_clade" = clade), by.x = "Compared_genome", by.y = "MyName", all.x = TRUE) %>% unique()
  print(df[df$Shared_pct > 5,] %>% group_by(Focal_genome, Focal_hap, Focal_chr, Chromosome_Window_Start) %>% dplyr::summarise(n_distinct(Compared_genome)) %>% filter(`n_distinct(Compared_genome)` > 1) %>% nrow())
  all_noncolinear = rbind(all_noncolinear,df)
}
all_noncolinear$clade = factor(all_noncolinear$clade, levels = rev(c("Iberian", "Muralis", "Siculus", "Sicilian-Maltese", "Western", "Balkan")))
all_noncolinear$Compared_clade = factor(all_noncolinear$Compared_clade, levels = rev(c("Iberian", "Muralis", "Siculus", "Sicilian-Maltese", "Western", "Balkan")))
all_noncolinear = all_noncolinear %>% mutate(queries = paste0(Focal_genome, "#", Focal_hap, "#", Focal_chr, ":", Chromosome_Window_Start, "-", Chromosome_Window_End)) 
map = read.delim("/home/zajac/chromosome_sets/community.chr.map")
all_noncolinear$community = map[match(all_noncolinear$Focal_chr, map$chromosome),]$community
all_noncolinear$Focal_chr = factor(all_noncolinear$Focal_chr, levels = c(as.character(seq(1,18,1)), "Z"))

#Filter is by how many species turned out to be similar to the focal species in each region:
xx = all_noncolinear %>% 
  group_by(Focal_genome, Focal_hap, Focal_chr,Chromosome_Window_Start, Chromosome_Window_End) %>% 
  summarise(n_distinct(Compared_genome), n_distinct(Compared_clade))
#Similarity to single genome
single_entries = xx[xx$`n_distinct(Compared_genome)` == 1 & xx$`n_distinct(Compared_clade)` == 1,]
single_entries_og = all_noncolinear[paste0(all_noncolinear$Focal_genome, "-", all_noncolinear$Focal_hap, "-", all_noncolinear$Focal_chr, "-", all_noncolinear$Chromosome_Window_Start, "-", all_noncolinear$Chromosome_Window_End) %in%  paste0(single_entries$Focal_genome, "-", single_entries$Focal_hap, "-", single_entries$Focal_chr, "-", single_entries$Chromosome_Window_Start, "-", single_entries$Chromosome_Window_End),]

#Similarity to multiple genomes but from same clade
single_clade = xx[xx$`n_distinct(Compared_genome)` > 1 & xx$`n_distinct(Compared_clade)` == 1,]
single_clade_og = all_noncolinear[paste0(all_noncolinear$Focal_genome, "-", all_noncolinear$Focal_hap, "-", all_noncolinear$Focal_chr, "-", all_noncolinear$Chromosome_Window_Start, "-", all_noncolinear$Chromosome_Window_End) %in%  paste0(single_clade$Focal_genome, "-", single_clade$Focal_hap, "-", single_clade$Focal_chr, "-", single_clade$Chromosome_Window_Start, "-", single_clade$Chromosome_Window_End),]
single_clade_og = single_clade_og %>% group_by(Focal_genome, Focal_hap, Focal_chr, Chromosome_Window_Start,Chromosome_Window_End) %>% slice_head(n=1) %>% ungroup()

#Similarity to multiple genomes but from different clades
multiple_entries = xx[xx$`n_distinct(Compared_genome)` > 1 & xx$`n_distinct(Compared_clade)` > 1,]
multiple_entries_og = all_noncolinear[paste0(all_noncolinear$Focal_genome, "-", all_noncolinear$Focal_hap, "-", all_noncolinear$Focal_chr, "-", all_noncolinear$Chromosome_Window_Start, "-", all_noncolinear$Chromosome_Window_End) %in%  paste0(multiple_entries$Focal_genome, "-", multiple_entries$Focal_hap, "-", multiple_entries$Focal_chr, "-", multiple_entries$Chromosome_Window_Start, "-", multiple_entries$Chromosome_Window_End),]

choose_mcr_first = multiple_entries_og[multiple_entries_og$Compared_genome == multiple_entries_og$Most_closely_related,]

multiple_entries_og = multiple_entries_og[!paste0(multiple_entries_og$Focal_genome, "-", multiple_entries_og$Focal_hap, "-", multiple_entries_og$Focal_chr, "-", multiple_entries_og$Chromosome_Window_Start, "-", multiple_entries_og$Chromosome_Window_End) %in%  paste0(choose_mcr_first$Focal_genome, "-", choose_mcr_first$Focal_hap, "-", choose_mcr_first$Focal_chr, "-", choose_mcr_first$Chromosome_Window_Start, "-", choose_mcr_first$Chromosome_Window_End),]

choose_same_clade = multiple_entries_og[multiple_entries_og$Compared_clade == multiple_entries_og$clade,]

multiple_entries_og = multiple_entries_og[!paste0(multiple_entries_og$Focal_genome, "-", multiple_entries_og$Focal_hap, "-", multiple_entries_og$Focal_chr, "-", multiple_entries_og$Chromosome_Window_Start, "-", multiple_entries_og$Chromosome_Window_End) %in%  paste0(choose_same_clade$Focal_genome, "-", choose_same_clade$Focal_hap, "-", choose_same_clade$Focal_chr, "-", choose_same_clade$Chromosome_Window_Start, "-", choose_same_clade$Chromosome_Window_End),]

multiple_entries_og = multiple_entries_og %>% group_by(Focal_genome, Focal_hap, Focal_chr, Chromosome_Window_Start,Chromosome_Window_End) %>% slice_head(n=1) %>% ungroup()

##Bind all of this selection back into a single data frame
all_noncolinear_modified = rbind(single_entries_og,single_clade_og,choose_mcr_first,choose_same_clade,multiple_entries_og)

##ABB BAB test
merge(all_noncolinear_modified[all_noncolinear_modified$Focal_genome == "rPodRaf1",] %>% 
        group_by(Compared_genome) %>% 
        summarise(windows = n_distinct(Focal_chr,Chromosome_Window_Start, Chromosome_Window_End)) %>% 
        mutate(sum(windows)), 
      all_noncolinear_modified[all_noncolinear_modified$Focal_genome == "rPodFil1",] %>% 
        group_by(Compared_genome) %>% 
        summarise(windows = n_distinct(Focal_chr,Chromosome_Window_Start, Chromosome_Window_End)) %>% 
        mutate(sum(windows)), by = "Compared_genome") %>% 
  mutate((`windows.x`/`sum(windows).x`)/(`windows.y`/`sum(windows).y`))  %>% 
  mutate((`windows.y`/`sum(windows).y`)/(`windows.x`/`sum(windows).x`))


#Check the values of Dinvestigate for the different regions coming from P.muralis and P.siculus
dsuite_files = list.files("../SYRI/Dsuite/Dinvestigate_50kb/", pattern = "rPodRaf1")[grepl("_100_100.txt", list.files("../SYRI/Dsuite/Dinvestigate_50kb/", pattern = "rPodRaf1"))]
trios=unique(gsub("_localFstats_community[0-9]+.Dinvestigate_100_100.txt", "", dsuite_files))
dsuite_results = NULL
for (i in trios){
  files = dsuite_files[grepl(i, dsuite_files)]
  tmp=NULL
  for (x in files){
    df2 = read.delim(paste0("../SYRI/Dsuite/Dinvestigate_50kb/",x))
    tmp = rbind(tmp,df2)
  }
  dsuite_results[[i]] = tmp
}
dsuite_results = lapply(dsuite_results, function(x){ x$incretensis = paste0(x$chr,":", x$windowStart, "-", x$windowEnd); return(x)})


all_disp = NULL
for (x in all_noncolinear_modified[all_noncolinear_modified$Focal_genome == "rPodRaf1",]$queries){
  chr=sapply(str_split(gsub(":", "#", x), "#"), .subset, 3)
  community = chr_map[chr_map$chromosome == chr,]$community
  cmd <- paste0("impg query -r ", x,
                " -p ", paste0("../pggb_allcommunities/community",community, ".s5000.paf"),
                "| grep rPodCre2.1")
  result <- system(cmd, intern = TRUE)
  if (length(result) == 0) {
    next
  } else if (length(result) == 1) {
    result <- paste0(sapply(str_split(result, "\t"), .subset, 2), "-", sapply(str_split(result, "\t"), .subset, 3))
    xm = data.frame(queries = x, incretensis= paste0("rPodCre2.1#1#",chr,":", result))
  } else {
    result <- paste0(sapply(str_split(result[1], "\t"), .subset, 2), "-",
                     sapply(str_split(result[max(length(result))], "\t"), .subset, 3))
    xm = data.frame(queries = x, incretensis= paste0("rPodCre2.1#1#",chr,":", result))}
  all_disp = rbind(all_disp,xm)
}
dsuite_raffonei_check = merge(all_noncolinear_modified[all_noncolinear_modified$Focal_genome == "rPodRaf1",], all_disp, by = "queries", all.x = TRUE)


dsuite_results_merged = lapply(dsuite_results, function(x){
  x = makeGRangesFromDataFrame(x, seqnames.field = "chr", start.field = "windowStart", end.field = "windowEnd", keep.extra.columns = T)
  df2gr = makeGRangesFromDataFrame(all_disp %>% separate(incretensis, into = c("seqnames", "q"), sep = ":") %>% separate(q, into = c("start", "end"), sep = "-"), keep.extra.columns = T)
  hits = findOverlaps(x, df2gr)
  query_hits   <- x[queryHits(hits)]
  subject_hits <- df2gr[subjectHits(hits)]
  overlaps <- pintersect(query_hits, subject_hits)
  final = data.frame(data.frame(query_hits), data.frame(subject_hits)) 
})

#Plot
dsuite_raffonei_check_fdm = lapply(dsuite_results_merged, function(x){ 
  df = merge(dsuite_raffonei_check, x[,c(8,16)], by = "queries", all.x = TRUE)
  return(df)})
pdf("../Figures/Dsuite_f_dM_Raffonei.pdf", width = 6, height = 6)
ggplot(bind_rows(dsuite_raffonei_check_fdm[c(1:3)], .id = "Comparison") %>% 
         mutate(Comparison = if_else(Comparison %in% c("rPodFil1_rPodRaf1_PodMur1", "rPodFil1_rPodRaf1_rPodMur119"), 
                                     "rPodFil1_rPodRaf1_Muralis", 
                                     "rPodFil1_rPodRaf1_rPodSic1")), 
       aes(Compared_clade, f_dM)) + 
  geom_boxplot(outliers = F) + 
  theme_bw() +
  facet_grid2(Comparison~., scales = "free", strip = strip_themed(background_y = elem_list_rect(fill = c("green","yellow")))) 
dev.off()

#Because the all_noncolinear and colinear dataframe is windows of 50kb while the inversions are less than that I will calculate the similarity of raffonei to other species within inversion separately (recompute for precision)
files_dispnodes_within_inv=list.files("../INVERSIONS/", pattern = "paf.dispnodes.csv", recursive = T, full.names = T)
final = NULL
for (x in files_dispnodes_within_inv){
  df = read.delim(x, sep = ",")
  df = df[grepl("rPodRaf1", df$Focal_Path),] %>% slice_max(Shared_pct, with_ties = T)
  final = rbind(final,df)
}
final = final %>% 
  separate(Compared_Path, into = c("Compared_Path", "coordinates"), sep = ":") %>%  
  separate(Compared_Path, into = c("MyName", "Compared_Hap", "Compared_CHr"), sep = "#") %>% 
  left_join(sample_data[,c(16,17)]) %>% 
  mutate(inv_ids = str_replace(str_replace(Focal_Path, ":", "_"), "-", "_")) %>% dplyr::select(13,8,12) 
final = final %>% separate(inv_ids, into = c("species", "subject_start", "subject_end"), sep = "_")

##Load trees built from whole sequence alignments
snp_trees = read.delim("../INVERSIONS/snp.data.trees.summary.txt", header = F)
snp_trees$V1 = str_remove(sapply(str_split(str_remove(snp_trees$V1, ".vcf.vcf.gz.min4.phy.iqtree"), "/"), .subset, 2), "inv_") 
snp_trees = snp_trees %>% separate(V1, into = c("chr", "coord"), ":") %>% separate(coord, into = c("start", "end"), "-")
snp_trees = merge(snp_trees, invs[grepl("rPodRaf1", invs$species) & invs$width > 1000, c("chr","start","end", "species", "subject_start", "subject_end")], by = c("chr","start","end"))
colnames(snp_trees) = c("chr","start","end", "clade", "bootstrap", "species", "subject_start", "subject_end")

#Check how many snps have been used per tree
nb_snps = NULL
for (i in list.files("../INVERSIONS/", pattern = "vcf.gz$", recursive = T, full.names = T)[grepl("inv", list.files("../INVERSIONS/", pattern = "vcf.gz$", recursive = T, full.names = T))]){
  df = read.vcfR(i)
  df = length(df@fix)
  nb_snps = rbind(nb_snps, data.frame(i, df))
}
mean(nb_snps$df)

##Combine all of this information into a single plot
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

all_invs_raffonei = invs[invs$species_list == "rPodRaf1#1",]
shared_raffonei = merge(shared[shared$rPodRaf1 == 1,-13] %>% 
                          pivot_longer(cols = c(3:15)) %>% 
                          filter(value > 0) %>% 
                          mutate(species_list = "rPodRaf1#1") %>% 
                          filter(name != "rPodRaf1"), 
                        sample_data[c(16,17)], by.x = "name", by.y = "MyName") %>% 
  dplyr::rename(presenceinotherclades = clade)
all_invs_raffonei = merge(all_invs_raffonei, shared_raffonei, by = c("species_list", "max_start", "min_end"), all.x = TRUE) %>% unique()
all_invs_raffonei = merge(all_invs_raffonei, final, by = c("species","subject_start","subject_end"), all = TRUE) %>% unique()

all_invs_raffonei$inv_ids = paste0(all_invs_raffonei$species,"_", all_invs_raffonei$subject_start, "_", all_invs_raffonei$subject_end)
all_invs_raffonei = all_invs_raffonei %>% mutate(presenceinotherclades = if_else(is.na(presenceinotherclades), "unique", presenceinotherclades))
all_invs_raffonei = all_invs_raffonei[all_invs_raffonei$width > 1000,]

presence_df = all_invs_raffonei[,c("inv_ids", "presenceinotherclades", "N.y")] 
presence_df = presence_df[presence_df$presenceinotherclades != "unique",]
presence_df$present <- TRUE
colnames(presence_df) = c("inv_ids", "clade", "NumberofSpecies", "present")
presence_df = unique(presence_df)
seq_df = all_invs_raffonei[,c("inv_ids", "clade", "Shared_pct")] 
seq_df = seq_df[!is.na(seq_df$clade),]
seq_df$dispsimilarity = TRUE
colnames(seq_df) = c("inv_ids", "clade", "DispNodePerc", "dispsimilarity")
seq_df = unique(seq_df)
trees_df = snp_trees
trees_df$inv_ids = paste0(snp_trees$species,"_", snp_trees$subject_start, "_", snp_trees$subject_end)
trees_df$snp_tree = TRUE
trees_df = trees_df[,c("inv_ids","clade", "snp_tree", "bootstrap")]
trees_df = trees_df[!is.na(trees_df$clade),]
colnames(trees_df) = c("inv_ids","clade", "snp_tree", "bootstrap")
inversion_order <- unique(all_invs_raffonei$inv_ids)
species_order <- c("Iberian", "Muralis", "Siculus", "Sicilian-Maltese", "Western", "Balkan")

hits_long = merge(merge(presence_df, seq_df, by = c("inv_ids", "clade"), all = TRUE), trees_df, by = c("inv_ids", "clade"), all = TRUE)
order = hits_long %>% mutate(Chr = sapply(str_split(str_replace(inv_ids, "_", "#"), "#"), .subset, 3)) %>% arrange(as.numeric(Chr)) %>% pull(inv_ids)
order = unique(order)
text = hits_long[,c(1,2,3,5,8)] %>% dplyr::rename(present = NumberofSpecies, dispsimilarity = DispNodePerc, snp_tree = bootstrap) %>% reshape2::melt(id.vars = c("inv_ids", "clade")) %>% dplyr::rename(category = variable) %>% mutate( n = if_else(is.na(value), 0, value), inv_ids = factor(inv_ids, levels = rev(order)), clade = factor(clade, levels = species_order),category = factor( category, levels = c("present", "dispsimilarity", "snp_tree") ))
hits_long = hits_long[,c(1,2,4,6,7)] %>% reshape2::melt(id.vars = c("inv_ids", "clade")) %>% dplyr::rename(category = variable) %>% mutate( n = if_else(is.na(value), 0, 1), inv_ids = factor(inv_ids, levels = rev(order)), clade = factor(clade, levels = species_order),category = factor( category, levels = c("present", "dispsimilarity", "snp_tree") )) 
hits_long = hits_long %>% mutate(category2 = if_else(inv_ids %in% setdiff(seq_df$inv_ids, presence_df$inv_ids), "Unique to Raffonei", "Shared with other species"))
text = text %>% mutate(category2 = if_else(inv_ids %in% setdiff(seq_df$inv_ids, presence_df$inv_ids), "Unique to Raffonei", "Shared with other species"))


p2 = ggplot( hits_long, aes(x = category, y = inv_ids)) +
  geom_tile(
    aes(fill = ifelse( n == 1, as.character(category), "empty")),
    color = "grey85",
    linewidth = 0.3
  ) + geom_text(data = text[text$n != 0,], aes(x = category, y = inv_ids, label = round(n,0))) +
  facet_grid(
    category2 ~ clade,
    scales = "free",
    space = "free"
  ) +
  scale_fill_manual(
    values = c(
      present = "#2C7FB8",
      dispsimilarity = "#F2B701",
      snp_tree = "purple",
      empty = "white"
    ),
    breaks = c("present", "dispsimilarity", "snp_tree"),
    labels = c("Matching breakpoints", "Dispensable nodes sharing", "SNP tree relatedness"),
    name = "Match type"
  ) +
  labs(
    x = NULL,
    y = "Inversion"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "grey90", color = "grey60"),
    strip.text.x = element_text(size = 8, face = "bold"),
    legend.position = "left",
    text = element_text(size = 16)
  ) 

p3 = hits_long[match(order, hits_long$inv_ids),c(1,6), drop = FALSE] %>% unique() %>% 
  mutate(ll = inv_ids) %>% 
  separate(ll, into = c("Seqname", "start", "end"), sep = "_") %>% 
  mutate(length = as.numeric(end)- as.numeric(start) ) %>% 
  ggplot(aes(length,inv_ids)) + 
  geom_col() + theme_bw() + 
  facet_grid(category2 ~ . ,space = "free", scales = "free") + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), text = element_text(size = 14)) + 
  scale_x_log10()

pdf("../Figures/Raffonei_full_analysis.pdf", width = 15, height = 9)
p2 + p3
dev.off()

##Check how the results compare if you overlap the inversions with the genome wide estimation in 50kb windows
check = NULL
for (i in unique(all_noncolinear_modified$Chr)) {
  
  chr_data = all_noncolinear_modified[all_noncolinear_modified$Chr == i, ]
  
  cld = unique(chr_data$clade)
  if (length(cld) != 1) {
    warning(sprintf("Chr '%s' has %d distinct 'clade' values (expected 1) - using the first one. Please check your data.", i, length(cld)))
    cld = cld[1]
  }
  
  windows_df = chr_data %>%
    dplyr::select(Compared_clade, queries, Window_length_bp) %>%
    separate(queries, into = c("seqnames", "ranges"), sep = ":") %>%
    separate(ranges, into = c("start", "end"), sep = "-", convert = TRUE)  # convert = TRUE -> numeric start/end
  
  windows = makeGRangesFromDataFrame(windows_df, keep.extra.columns = TRUE)
  
  inv_sub = invs[invs$species == i, ] %>%
    dplyr::select(species, subject_start, subject_end, shared_between_species)
  
  final = windows_df[0, ]  # empty template; filled in below if there are overlaps
  
  if (nrow(inv_sub) > 0) {
    inversions = makeGRangesFromDataFrame(inv_sub,
                                          seqnames.field = "species",
                                          start.field    = "subject_start",
                                          end.field      = "subject_end",
                                          keep.extra.columns = TRUE)
    
    hits = findOverlaps(windows, inversions)
    
    if (length(hits) > 0) {
      query_hits   = windows[queryHits(hits)]
      subject_hits = inversions[subjectHits(hits)]
      overlaps     = pintersect(query_hits, subject_hits)
      
      final = unique(data.frame(data.frame(query_hits), data.frame(subject_hits),
                                overlap_length = width(overlaps)))
    }
  }
  
  ## ---- record which (Chr, window, Compared_clade) combos overlap an inversion ----
  if (nrow(final) > 0) {
    check = rbind(check, final)
  }
  
}
check[grepl("rPodRaf1#1", check$seqnames.1),] %>% 
  group_by(Compared_clade,seqnames.1,start.1,end.1, width.1) %>% 
  summarise(sum(overlap_length)) %>% 
  group_by(perc = `sum(overlap_length)`*100/width.1) %>% 
  group_by(seqnames.1,start.1,end.1, width.1) %>% 
  slice_max(perc, with_ties = TRUE) %>% data.frame()

## Plot karyotype
library(RIdeogram)
bed_raf = beds[beds$Focal_genome == "rPodRaf1",]
bed_raf = bed_raf %>% mutate("Chr" = Focal_chr)
bed_raf = bed_raf[,c("Chr", "Start", "End")]
bed_raf = bed_raf[bed_raf$Chr != "W",]
mosaic = all_noncolinear_modified[all_noncolinear_modified$Focal_genome == "rPodRaf1",c("Focal_genome" ,"Focal_hap" ,"Focal_chr" ,"Chromosome_Window_Start", "Chromosome_Window_End", "Compared_clade")] %>% mutate("Chr" = Focal_chr, Value = case_when(Compared_clade == "Sicilian-Maltese" ~ 0, Compared_clade == "Muralis" ~ 50, Compared_clade == "Siculus" ~ 100, .default = 0)) %>% dplyr::select(Chr, Chromosome_Window_Start, Chromosome_Window_End, Value) 
colnames(mosaic) = c("Chr", "Start", "End", "Value")
mosaic$Chr <- unname(as.character(mosaic$Chr))
ideogram(karyotype = bed_raf, overlaid = mosaic, output = "../Figures/Raffonei.mosaic.karyotype.svg", colorset1 = c("cyan4", "green", "yellow"))
