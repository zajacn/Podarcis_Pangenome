##This script is for analysis of gene flow across all species

library(ggplot2)
library(ggh4x)


chr_map = read.delim("../chromosome_sets/community.chr.map")

#Load bed files for all species
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

#Now load in the disp node sharing across chromosomes
all_noncolinear = NULL
for (i in seq(0,18,1)){
  df = read.delim(paste0("../PhyloTree/GenomicMosaicism/Disp.node.sharing.combined.50kb.community", i, ".csv"), sep = ",")
  df  = df %>% 
    separate(Compared_Path, into = c("Compared_genome", "Compared_hap", "Compared_chr"), sep = "#") %>% 
    separate(Focal_Path, into = c("Focal_genome", "Focal_hap", "Focal_chr"), sep = "#")
  df = df %>% 
    group_by(Focal_genome, Focal_hap, Focal_chr,Chromosome_Window_Start, Chromosome_Window_End) %>% 
    slice_max(as.numeric(Shared_pct), n = 1, with_ties = TRUE) %>% 
    unique() %>% 
    left_join(beds, by = c("Focal_genome", "Focal_hap", "Focal_chr")) %>% 
    mutate(blocks = round(End/8,0))
  df = df %>%
    group_by(Focal_genome, Focal_hap, Focal_chr) %>%
    group_modify(~ {
      breaks <- seq(unique(.x$blocks), unique(.x$End) - unique(.x$blocks), length.out = 15)
      .x %>%
        mutate(block = ezCut(Chromosome_Window_End, 
                             breaks, labels = as.character(seq(1,16,1))))
    }) %>%
    ungroup()
  df = merge(df, 
             df %>% group_by(Focal_genome, Focal_hap, Focal_chr, Compared_genome) %>% 
               dplyr::summarise(n_distinct(Chromosome_Window_Start)) %>% 
               group_by(Focal_genome, Focal_hap, Focal_chr) %>% 
               slice_max(as.numeric(`n_distinct(Chromosome_Window_Start)`), n = 1, with_ties = TRUE) %>% 
               select(1:4) %>% 
               dplyr::rename(Most_closely_related = Compared_genome ), 
             by = c("Focal_genome", "Focal_hap", "Focal_chr"), 
             all.x = TRUE)
  df = merge(df, sample_data[,c(3,16,17)], by.x = "Focal_genome", by.y = "MyName", all.x = TRUE) %>% unique()
  df = merge(df, sample_data[,c(16,17)] %>% dplyr::rename("Compared_clade" = clade), by.x = "Compared_genome", by.y = "MyName", all.x = TRUE) %>% unique()
  df = df %>% 
    mutate(window_description = case_when(Compared_genome == Most_closely_related ~ "phylogenetically_concordant",
                                          Compared_genome != Most_closely_related & clade == Compared_clade ~"phylogenetically_concordant", 
                                          .default = "phylogenetically_discordant"))
  print(df[df$Shared_pct > 5,] %>% group_by(Focal_genome, Focal_hap, Focal_chr, Chromosome_Window_Start) %>% dplyr::summarise(n_distinct(Compared_genome)) %>% filter(`n_distinct(Compared_genome)` > 1) %>% nrow())
  all_noncolinear = rbind(all_noncolinear,df)
}
all_noncolinear$clade = factor(all_noncolinear$clade, levels = c("Iberian", "Muralis", "Siculus", "Sicilian-Maltese", "Western", "Balkan"))
all_noncolinear$Compared_clade = factor(all_noncolinear$Compared_clade, levels = c("Iberian", "Muralis", "Siculus", "Sicilian-Maltese", "Western", "Balkan"))

#Filter the results
#Check if mutliple species equally related:
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
all_noncolinear_modified$queries = paste0(all_noncolinear_modified$Focal_genome, "#", 
                                          all_noncolinear_modified$Focal_hap, "#", 
                                          all_noncolinear_modified$Focal_chr, ":", 
                                          all_noncolinear_modified$Chromosome_Window_Start, "-", 
                                          all_noncolinear_modified$Chromosome_Window_End)

##Plot matrix of shared ancestry

p = all_noncolinear_modified %>% 
  group_by(Focal_genome,clade,Compared_genome, Compared_clade) %>% 
  summarise(sum(Shared_bp)) %>% group_by(Focal_genome,clade) %>%  mutate(total = sum(`sum(Shared_bp)`)) %>% mutate(perc = `sum(Shared_bp)` *100/total) %>% filter(perc >= 1) %>% ggplot(aes(Focal_genome,Compared_genome, fill = perc)) + 
  geom_tile() + 
  scale_fill_viridis_c() + 
  geom_text(aes(label = round(perc,0)), colour = "white", size = 9) + 
  facet_grid2(Compared_clade~clade, space = "free", scales = "free", strip = strip_themed(
    background_x = elem_list_rect(fill = c("magenta2", "green","yellow", "cyan4","blue","brown")),
    background_y = elem_list_rect(fill = c("magenta2", "green","yellow", "cyan4","blue","brown")))
  ) + theme_bw()

p$data$Focal_genome = sample_data[match(p$data$Focal_genome, sample_data$MyName),]$Organism.Name
p$data$Compared_genome = sample_data[match(p$data$Compared_genome, sample_data$MyName),]$Organism.Name

p = p + 
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1, face = "italic"),
        axis.text.y = element_text(size = 13, face = "italic"), 
        text = element_text(size = 13)) + 
  labs(x = "Focal genome", y = "Compared genome", fill = "Percentage")

pdf("../Figures/Matrix_shared_ancestry.pdf", width = 14, height = 10)
p
dev.off()


## Check for distribution of phylogenetically discordant windows across chromosomes of the most admixed species
p1 = all_noncolinear_modified %>% filter(Focal_genome %in% c("rPodRaf1", "rPodFil1", "rPodTil1", "rPodSic1")) %>% 
  group_by(Focal_genome, Focal_hap, Focal_chr, block, window_description, Organism.Name) %>% 
  dplyr::summarise(windows = n_distinct(Chromosome_Window_Start)) %>% group_by(Focal_genome, Focal_hap, Focal_chr, block) %>% mutate(total_windows = sum(windows)) %>% 
  mutate(windows = windows * 100/total_windows) %>% 
  group_by(Focal_genome, Focal_chr, block, window_description, Organism.Name)  %>% 
  dplyr::summarise(windows = mean(windows)) %>% filter(window_description == "phylogenetically_discordant") %>%
  ggplot(aes(block, windows)) + 
  geom_violin(fill = "darkblue") + 
  facet_grid2(Organism.Name~., scales = "free", strip = strip_themed(,
                                                                     background_y = elem_list_rect(fill = c("cyan4", "cyan4","yellow","blue")))
  ) + 
  theme_bw() + 
  labs(x = "Chromosomal block", y = "Percentage of windows") +
  theme(text = element_text(size = 18),
        legend.position = "top")

pdf("../Figures/Distribution_discordant_concordant_regions.genome.pdf", width = 9, height = 9)
p1
dev.off()

#Load SVs and block
all_comb = NULL
for (i in list.files("../SV_intersects/largeSVs/", pattern = "merged.DEFAULT")[!grepl("annot", list.files("../SV_intersects/largeSVs/", pattern = "merged.DEFAULT"))]){
  comb = read.vcfR(paste0("../SV_intersects/largeSVs/", i))
  comb = cbind(data.frame(comb@fix), data.frame(comb@gt))
  comb = comb %>% dplyr::mutate(across(10:23, ~ str_remove(., ":.*"))) %>% 
    mutate(REFLEN = nchar(REF), ALTLEN = nchar(ALT)) %>% 
    mutate(SVTYPE = sub(".*SVTYPE=([^;]+).*", "\\1", INFO),
           TYPE=sub(".*;TYPE=([^;]+).*", "\\1", INFO),
           AVG_LEN = sub(".*AVG_LEN=([^;]+).*", "\\1", INFO), 
           AVG_START = sub(".*AVG_START=([^;]+).*", "\\1", INFO), 
           AVG_END = sub(".*AVG_END=([^;]+).*", "\\1", INFO)) %>% 
    mutate(SVTYPE = case_when(SVTYPE == "INS" ~ "INS", 
                              SVTYPE == "DEL" ~ "DEL", 
                              (!SVTYPE %in% c("INS", "DEL")) &  !is.na(TYPE) ~ TYPE, 
                              (!SVTYPE %in% c("INS", "DEL")) &  is.na(TYPE) ~ "NA"))
  comb = comb[,c(1,2,3,10:27)]
  comb = comb %>% dplyr::mutate(across(colnames(comb)[4:17], ~ case_when(. == "0|0" ~ 0, . == "0" ~ 0, . == "1|1" ~ 1, . == "1" ~ 1,  . == "0|1" ~ 1, . == "1|0" ~ 1, . == ".|1" ~ 1, . == "0|."~ 0, . == ".|0" ~ 0, . == "1|." ~ 1 )))
  all_comb[[i]] = comb
}
all_comb = bind_rows(all_comb)
all_comb$nb_species = rowSums(all_comb[,c(4:17)] > 0)
for_karyotype = all_comb[all_comb$nb_species >=1,]
kary = NULL
for (i in unique(for_karyotype$CHROM)){
  tmp = for_karyotype[for_karyotype$CHROM ==i,]
  bed1 = bed[bed$Chr ==i,]$End
  tmp$blocks = ezCut(x = as.numeric(tmp$POS), breaks = seq(1,as.numeric(bed1), 500000))
  kary = rbind(kary,tmp)
}
kary$POS = as.numeric(kary$POS)
kary = kary %>%
  group_by(CHROM) %>%
  mutate(
    max_pos = max(POS, na.rm = TRUE),
    blocks = if_else(
      str_detect(blocks, "^>\\s*\\d+"),
      str_replace(blocks, "^>\\s*(\\d+)", "(\\1 - ") %>%
        paste0(max_pos, "]"),
      blocks
    )
  ) %>%
  select(-max_pos)
kary = kary[,c(1,2,23)] %>% group_by(blocks, CHROM) %>% summarise(n_distinct(POS)) %>% mutate(Start =as.numeric(sapply(str_split(str_remove(as.character(blocks), "\\("), " - "), .subset, 1)), End = as.numeric(sapply(str_split(str_remove(as.character(blocks), "\\]"), " - "), .subset, 2))) %>% ungroup() %>% dplyr::select(c(2,4,5,3))
colnames(kary) = c("Chr", "Start", "End", "Value")

## And load inversion data and mark if shared between species
invs = read.delim("../INVERSIONS/inversions.clustered.0.9proc.mapq.mapc.filtered.txt")
shared = invs[invs$is_singleton == "FALSE",] %>% 
  dplyr::select(3,15,16,17) %>% 
  separate(species_list, into = c("species", "haplotype"), sep = "#") %>%
  dplyr::select(1,3,4,5) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = mean_length) %>% 
  replace(is.na(.),0) %>% mutate(across(c(3:16), ~ if_else(. > 0, 1,0))) %>%
  mutate(N = rowSums(.[3:16])) 
invs$shared_between_species = ""
invs[paste0(invs$max_start, "-", invs$min_end) %in% paste0(shared[shared$N >1,]$max_start,"-",shared[shared$N >1,]$min_end),]$shared_between_species = "Yes"
  
##Compute proportion of phylogenetically discordant nodes in different regions in the genome (compared to phylogenetically concordant nodes)
###high_SV_regions

all_noncolinear_modified$Chr = sapply(str_split(all_noncolinear_modified$queries, ":"), .subset,1)
high_SV_regions = all_noncolinear_modified %>%
  group_by(Chr) %>%
  mutate(.cutoff = sort(Focal_disp_bp, decreasing = TRUE)[ceiling(n() * 0.20)]) %>%
  filter(Focal_disp_bp >= .cutoff) %>%
  ungroup() %>%
  dplyr::select(-.cutoff)

#Check where these regions are 
ggplot(high_SV_regions[high_SV_regions$Focal_genome == "PodMur1",], aes(x=Chromosome_Window_Start, xend = Chromosome_Window_End, y = 1, yend = 1)) + geom_segment(linewidth = 2) + facet_grid(Chr~.)


###Overlap phylogenetically discordant nodes with inversions and with high SV regions
phylogenetic_discordance_two_categories = NULL
inversion_overlap_keys = NULL 
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
    inversion_overlap_keys = rbind(
      inversion_overlap_keys,
      unique(data.frame(Chr = i, start = final$start, end = final$end,
                        Compared_clade = final$Compared_clade))
    )
  }
  
  ## ---- category 1: within inversions ----
  ## NB: numerator uses full Window_length_bp (not the clipped overlap width),
  ## so it's on the same scale as the other two categories below.
  if (nrow(final) > 0) {
    within_inversions = final %>%
      distinct(seqnames, start, end, Compared_clade, shared_between_species, Window_length_bp) %>%
      mutate(desc = if_else(Compared_clade == cld, "Phylogenetically concordant", "Phylogenetically discordant")) %>%
      group_by(Compared_clade, shared_between_species, desc) %>%
      summarise(windows = sum(Window_length_bp), .groups = "drop") %>%
      group_by(shared_between_species) %>%
      mutate(total = sum(windows)) %>%   # denominator computed *within* the shared/unique subset
      ungroup() %>%
      mutate(Chr = i, clade = cld, category2 = "within inversions")
  } else {
    within_inversions = NULL
  }
  
  ## ---- category 2: highly variable regions ----
  ## Exclude any window already counted in "within inversions" for this Chr,
  ## matching on window coordinates AND Compared_clade (since desc depends on both).
  hv = high_SV_regions[high_SV_regions$Chr == i, ]
  
  if (nrow(final) > 0) {
    hv = anti_join(hv, final,
                   by = c("Chromosome_Window_Start" = "start",
                          "Chromosome_Window_End"   = "end",
                          "Compared_clade"          = "Compared_clade"))
  }
  
  highly_variable_regions = hv %>%
    mutate(desc = if_else(Compared_clade == clade, "Phylogenetically concordant", "Phylogenetically discordant")) %>%
    group_by(Chr, Compared_clade, desc, clade) %>%
    summarise(windows = sum(Window_length_bp), .groups = "drop") %>%
    group_by(Chr) %>%
    mutate(total = sum(windows)) %>%
    ungroup() %>%
    mutate(category2 = "highly variable regions")
  
  out = dplyr::bind_rows(within_inversions, highly_variable_regions)
  phylogenetic_discordance_two_categories = dplyr::bind_rows(phylogenetic_discordance_two_categories, out)
}

#"Rest of chromosome" = everything NOT in high_SV_regions AND NOT
rest_of_chr = all_noncolinear_modified %>%
  anti_join(high_SV_regions,
            by = c("Chr", "Chromosome_Window_Start", "Chromosome_Window_End", "Compared_clade"))

if (!is.null(inversion_overlap_keys)) {
  rest_of_chr = rest_of_chr %>%
    anti_join(inversion_overlap_keys,
              by = c("Chr" = "Chr",
                     "Chromosome_Window_Start" = "start",
                     "Chromosome_Window_End"   = "end",
                     "Compared_clade"          = "Compared_clade"))
}

rest_of_chr = rest_of_chr %>%
  mutate(desc = if_else(Compared_clade == clade, "Phylogenetically concordant", "Phylogenetically discordant")) %>%
  group_by(Chr, Compared_clade, desc, clade) %>%
  summarise(windows = sum(Window_length_bp), .groups = "drop") %>%
  group_by(Chr) %>%
  mutate(total = sum(windows)) %>%
  ungroup() %>%
  mutate(category2 = "rest of chromosome") 

#Combine all three categories
phylogenetic_discordance_three_categories = dplyr::bind_rows(phylogenetic_discordance_two_categories, rest_of_chr)

# Manual correction: rPodSic1 vs Muralis comparisons are treated as concordant
phylogenetic_discordance_three_categories = phylogenetic_discordance_three_categories %>%
  mutate(desc = if_else(grepl("rPodSic1", Chr) & Compared_clade == "Muralis",
                        "Phylogenetically concordant", desc))

phylogenetic_discordance_three_categories = phylogenetic_discordance_three_categories %>%
  mutate(species = sapply(str_split(Chr, "#"), .subset, 1)) %>%
  mutate(category2 = case_when(
    category2 == "within inversions" & shared_between_species == "Yes" ~ "within inversions shared between species",
    category2 == "within inversions" & shared_between_species == ""   ~ "within inversions unique to each species",
    TRUE ~ category2
  ))

#Group into the desc categories and plot
my_plot = phylogenetic_discordance_three_categories %>%
  group_by(Chr, shared_between_species, desc, total, category2, clade, species) %>%
  summarise(windows = sum(windows), .groups = "drop") %>%
  mutate(perc = windows * 100 / total) %>%
  filter(desc == "Phylogenetically discordant")

my_plot = merge(my_plot, sample_data[, c(3, 17)], by.x = "species", by.y = "MyName") %>%
  unique() %>%
  filter(species != "rPodCre2.1")

my_plot_data = my_plot
comparisons_list = combn(as.character(unique(my_plot_data$category2)), 2, simplify = FALSE)
p = ggplot(my_plot_data, aes(category2, perc, fill = category2)) +
  geom_boxplot() + 
  facet_wrap2(~clade, strip = strip_themed(
    background_x = elem_list_rect(fill = c("magenta2", "green","yellow", "cyan4","blue","brown")))) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = "% phylogenetically discordant windows")

my_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


pdf("../Figures/Proportion_of_phylogenetically_disc_per_feature.pdf", width = 6, height = 10)
p + my_theme + guides(fill=guide_legend(ncol=2)) 
dev.off()

##Network plots
library(igraph)
## keep one row per (Chr, Compared_clade, clade, category2, desc) - no info lost yet
disc_conc = phylogenetic_discordance_three_categories %>%
  distinct(Chr, Compared_clade, clade, category2, desc, windows, total)

## Denominator: total windows for a given (clade, category2), pooled across its chromosomes.
## `total` is per-Chr (identical across Compared_clade rows for that Chr), so dedupe by Chr first
## to avoid summing the same chromosome's total multiple times.
denom = disc_conc %>%
  distinct(Chr, clade, category2, total) %>%
  group_by(clade, category2) %>%
  summarise(total = sum(total), .groups = "drop")

## Self-loops: concordant windows where clade == Compared_clade (own-clade agreement)
self_loops = disc_conc %>%
  filter(clade == Compared_clade, desc == "Phylogenetically concordant") %>%
  group_by(clade, category2) %>%
  summarise(windows = sum(windows), .groups = "drop") %>%
  left_join(denom, by = c("clade", "category2")) %>%
  mutate(perc = windows * 100 / total,
         Compared_clade = clade)   # from == to -> self-loop

## Cross-clade edges: discordant windows where clade != Compared_clade
## (this naturally EXCLUDES the rPodSic1/Muralis override rows, since those were
## manually recoded to "concordant" - so they reduce this edge's weight without
## being misattributed to a self-loop or elsewhere)
cross_edges = disc_conc %>%
  filter(clade != Compared_clade, desc == "Phylogenetically discordant") %>%
  group_by(clade, Compared_clade, category2) %>%
  summarise(windows = sum(windows), .groups = "drop") %>%
  left_join(denom, by = c("clade", "category2")) %>%
  mutate(perc = windows * 100 / total)

for_igraph = dplyr::bind_rows(self_loops, cross_edges)

edges_rest        = for_igraph %>% filter(category2 == "rest of chromosome") %>%
  select(from = clade, to = Compared_clade, weight = perc)
edges_within_sp    = for_igraph %>% filter(category2 == "within inversions shared between species") %>%
  select(from = clade, to = Compared_clade, weight = perc)
edges_within_uniq  = for_igraph %>% filter(category2 == "within inversions unique to each species") %>%
  select(from = clade, to = Compared_clade, weight = perc)
edges_outside      = for_igraph %>% filter(category2 == "highly variable regions") %>%
  select(from = clade, to = Compared_clade, weight = perc)

filter_edges = function(df, thresh = 5) {
  df %>% filter(from == to | weight > thresh)   # always keep self-loops; threshold only cross-edges
}

edges_rest        = filter_edges(edges_rest)
edges_within_sp    = filter_edges(edges_within_sp)
edges_within_uniq  = filter_edges(edges_within_uniq)
edges_outside      = filter_edges(edges_outside)
g_within_uniq <- graph_from_data_frame(edges_within_uniq, directed = FALSE)
g_within_sp   <- graph_from_data_frame(edges_within_sp,   directed = FALSE)
g_outside     <- graph_from_data_frame(edges_outside,     directed = FALSE)
g_rest        <- graph_from_data_frame(edges_rest,        directed = FALSE)


par(mfrow = c(2, 2),
    mar = c(3, 3, 3, 3),
    xpd = TRUE)

plot(g_within_uniq,
     edge.width = E(g_within_uniq)$weight / max(E(g_within_uniq)$weight) * 5,
     edge.label = round(E(g_within_uniq)$weight, 1),
     vertex.size = 30,
     vertex.color = "lightblue",
     vertex.label.cex = 0.8,
     main = "Within Inversions unique to each species")

plot(g_within_sp,
     edge.width = E(g_within_sp)$weight / max(E(g_within_sp)$weight) * 5,
     edge.label = round(E(g_within_sp)$weight, 1),
     vertex.size = 30,
     vertex.color = "lightblue",
     vertex.label.cex = 0.8,
     main = "Within Inversions shared between species")

plot(g_outside,
     edge.width = E(g_outside)$weight / max(E(g_outside)$weight) * 5,
     edge.label = round(E(g_outside)$weight, 1),
     vertex.size = 30,
     vertex.color = "lightgreen",
     vertex.label.cex = 0.8,
     main = "Within high SV regions")

plot(g_rest,
     edge.width = E(g_rest)$weight / max(E(g_rest)$weight) * 5,
     edge.label = round(E(g_rest)$weight, 1),
     vertex.size = 30,
     vertex.color = "lightpink",
     vertex.label.cex = 0.8,
     main = "Rest of Chromosome length")

## Check if the inversions unique to each species don't consist of mostly private sequence
invs_annot = read.delim("../INVERSIONS/inversions.clustered.0.9proc.mapq.mapc.filtered.node_annotation.txt")
for (i in c("private", "core")){
  print(invs_annot[invs_annot$shared_between_species == "Yes",] %>% replace(is.na(.), 0) %>% pull(i) %>% mean())
  print(invs_annot[invs_annot$shared_between_species == "Yes",] %>% replace(is.na(.), 0) %>% pull(i) %>% median())
  print(invs_annot[invs_annot$shared_between_species == "",] %>% replace(is.na(.), 0) %>% pull(i) %>% mean())
  print(invs_annot[invs_annot$shared_between_species == "",] %>% replace(is.na(.), 0) %>% pull(i) %>% median())
}

#First check if density of inversions and SVs and proportion of phylogenetically discordant regions correlates
library(ggrepel)
inversion_density = invs %>% 
  separate(species_list, into = c("Focal_genome", "Focal_hap"), sep = "#") %>% 
  group_by(Chrom, clusters.membership, Focal_genome, Focal_hap) %>% 
  summarise(n_distinct(start, end)) %>% 
  ungroup() %>% 
  group_by(Chrom, clusters.membership) %>% 
  summarise(n_distinct(Focal_genome)) %>% 
  ungroup() %>% 
  group_by(Chrom) %>% 
  summarise(total_inversions = sum(`n_distinct(Focal_genome)`)) %>% 
  mutate(Focal_chr = str_remove(Chrom, "rPodCre2.1#1#")) %>% 
  left_join(beds[beds$Focal_genome == "rPodCre2.1",], by = "Focal_chr") %>% 
  mutate(density = total_inversions * 1000000/End) %>% dplyr::select(c(3,8)) %>% mutate(cat = "inversions")

sv_density = kary %>% group_by(Chr) %>% summarise(sum(Value)) %>% 
  mutate(Focal_chr = str_remove(Chr, "rPodCre2.1#1#")) %>% 
  left_join(beds[beds$Focal_genome == "rPodCre2.1",], by = "Focal_chr") %>% 
  mutate(density = `sum(Value)` * 1000000/End) %>% dplyr::select(c(3,8)) %>% mutate(cat = "SVs")

hvr = rbind(inversion_density, sv_density)

pdf("../Figures/Phylogenetically_disc_disp_vs_inv_sv_density.pdf", width = 9, height = 8)
all_noncolinear_modified %>% 
  group_by(Focal_genome, Focal_hap, Focal_chr, window_description) %>% 
  summarise(wind = n_distinct(Chromosome_Window_Start, Chromosome_Window_End)) %>% 
  group_by(Focal_genome, Focal_hap, Focal_chr) %>%
  mutate(total = sum(wind)) %>% 
  mutate(perc = wind*100/total) %>% 
  group_by(Focal_chr, window_description) %>% 
  summarise(perc = mean(perc)) %>% 
  left_join(hvr, by = "Focal_chr") %>% 
  filter(window_description == "phylogenetically_discordant") %>% 
  ggplot(aes(perc, density)) + facet_grid(cat~., scales = "free_y") +
  geom_point() + 
  geom_text_repel(aes(label = Focal_chr)) + 
  geom_smooth(method = "lm") + stat_cor() + theme_bw()
dev.off()



