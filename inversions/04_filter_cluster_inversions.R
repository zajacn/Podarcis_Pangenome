#This script filters inversions for analysis
library(GenomicRanges)
library(utils)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ezRun)
library(data.table)
library(DT)
library(igraph)
library(RIdeogram)

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

##Load inversion data
all_invs = read.delim("/home/zajac/INVERSIONS/inversions.syri.impg.filtered.txt")

##Load mapping QC
mapc = read.delim("/home/zajac/INVERSIONS/temp/mapc.scan.of.syri.inversions.txt")
all_invs = merge(all_invs, mapc, by = c("chr", "start", "end", "species"), all.x = TRUE)

##Filter inversions by mapping cov 
all_invs_filtered = all_invs %>% filter(mean_coverage > 0.25 & mean_coverage < 20)
write_delim(all_invs_filtered, file = "INVERSIONS/inversions.syri.impg.mapc.filtered.txt")

##Cluster the inversions to find common inversions across species
find_overlapping_inversions <- function(df_list, min_overlap) {
  
  # Convert to GRanges and extract data
  gr_list <- lapply(names(df_list), function(sp) {
    df <- df_list[[sp]]
    GRanges(seqnames = df$chr, 
            ranges = IRanges(start = df$start, end = df$end),
            species = sp,                 
            subject_start = df$start_species_syri,     
            subject_end = df$end_species_syri)
  })
  
  all_inversions <- do.call("c", gr_list)
  
  # Create data.table
  dt <- data.table(
    idx = 1:length(all_inversions),
    chr = as.character(seqnames(all_inversions)),
    start = GenomicRanges::start(all_inversions),
    end = GenomicRanges::end(all_inversions),
    width = GenomicRanges::width(all_inversions),
    species = all_inversions$species,
    subject_start = all_inversions$subject_start,     
    subject_end = all_inversions$subject_end
  )
  
  setkey(dt, chr, start, end)
  setnames(dt, c("start", "end"), c("start", "end"))
  
  cat("Total inversions:", nrow(dt), "\n")
  cat("Finding overlaps...\n")
  
  # Self-join to find overlaps
  overlaps <- foverlaps(
    dt, dt,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "start", "end"),
    type = "any",
    nomatch = NULL
  )
  
  # Remove self-matches and duplicates
  overlaps <- overlaps[idx != i.idx]
  overlaps <- overlaps[idx < i.idx]  # Keep only unique pairs
  
  # Calculate reciprocal overlap
  overlaps[, int_start := pmax(start, i.start)]
  overlaps[, int_end := pmin(end, i.end)]
  overlaps[, int_width := int_end - int_start + 1]
  overlaps[, overlap_1 := int_width / width]
  overlaps[, overlap_2 := int_width / i.width]
  
  # Filter by threshold
  overlaps <- overlaps[overlap_2 >= min_overlap & overlap_1 >= min_overlap,]
  
  cat("Found", nrow(overlaps), "overlapping pairs\n")
  
  # Build clusters
  cluster_id <- 1:nrow(dt)
  
  if (nrow(overlaps) > 0) {
    edges <- data.frame(from = overlaps$idx, to = overlaps$i.idx)
    g <- graph_from_data_frame(
      edges,
      directed = FALSE,
      vertices = data.frame(name = dt$idx)
    )
    clusters <- igraph::components(g)
    
    node_ids <- as.numeric(V(g)$name)
    cluster_id[node_ids] <- clusters$membership
    
    # Only merge if clusters exists
    dt <- merge(dt, 
                data.frame(clusters.membership = clusters$membership) %>% 
                  rownames_to_column() %>% 
                  mutate(rowname = as.integer(rowname)), 
                by.x = "idx", by.y = "rowname", all = TRUE)
  } else {
    # No overlaps, add NA column
    dt$clusters.membership <- NA_integer_
  }
  dt[, cluster := ifelse(is.na(clusters.membership),
                         idx,
                         clusters.membership)]
  cluster_sizes <- dt[, .N, by = cluster]
  dt <- merge(dt, cluster_sizes, by = "cluster")
  dt[, is_singleton := (N == 1)]
  
  # Summarize
  summary_dt <- dt[, .(
    n_species = uniqueN(species),
    species_list = paste(unique(species), collapse = ", "),
    chr = dplyr::first(chr),
    mean_start = as.integer(mean(start)),
    mean_end = as.integer(mean(end)),
    min_start = min(start),
    max_start = max(start),
    min_end = min(end),
    max_end = max(end),
    mean_length = as.integer(mean(width)),
    min_length = min(width),
    max_length = max(width),
    n_inversions = .N
  ), by = clusters.membership][order(-n_species)]
  
  cat("Done!\n")
  
  return(list(
    detailed = as.data.frame(dt),
    shared_only = as.data.frame(summary_dt[!is.na(clusters.membership) & n_species > 1])
  ))
}

all_clustering=NULL
all_detailed = NULL
singletons = NULL
for (pid in c(0.5, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8,0.85, 0.90, 0.95, 1)){ #Check for 50-100% overlap
  clustering = NULL
  detailed = NULL
  for (chr in unique(all_invs_filtered$chr)){
    print(chr)
    subset = all_invs_filtered[all_invs_filtered$chr == chr,]
    df_list=split(subset, subset$species)
    df_list = lapply(df_list, function(x) unique(x[,c(1:6)]))
    print(names(df_list))
    results <- find_overlapping_inversions(df_list, min_overlap=pid)
    stopifnot(
      nrow(results$detailed) ==
        nrow(unique(results$detailed[, c("chr", "start", "end", "species")]))
    )
    singletons = rbind(singletons, results$detailed %>% filter(is_singleton) %>% dplyr::select(3,4,5,6,7)  %>% mutate(pid = pid))
    if ("shared_only" %in% names(results) && !is.null(results$shared_only)) {
      clustering[[chr]] = results$shared_only %>% separate_rows(species_list, sep = ", ") %>% dplyr::select(1,3,8,9,11) %>% 
        mutate(species_list = paste0(sapply(str_split(species_list, "#"), .subset ,1), "#", sapply(str_split(species_list, "#"), .subset ,2)))
      detailed[[chr]] = results$detailed
    }
  }
  name=paste0("min_overlap", pid)
  all_clustering[[name]] = clustering
  all_detailed[[name]] = detailed
}

##Select 85%,90% and 95% sequence overlap to save. Primary choice 90%
comb = merge(all_invs_filtered, bind_rows(all_detailed$min_overlap0.9, .id = "Chrom")[,c(3,4:6,8:13)], by.x = c("chr","start", "end", "species", "start_species_syri", "end_species_syri"), by.y = c("chr","start", "end", "species" ,"subject_start", "subject_end")) 
comb$species_list = paste0(sapply(str_split(comb$species, "#"), .subset, 1), "#", sapply(str_split(comb$species, "#"), .subset, 2))
comb = merge(comb, bind_rows(all_clustering$min_overlap0.9, .id = "chr"), by=c("clusters.membership", "chr", "species_list"), all.x = TRUE)
write_delim(comb, "/home/zajac/INVERSIONS/inversions.syri.impg.mapc.filtered.clustered.0.9.txt", delim = "\t", quote = "none")

##Add info on core, private, dispensable node density within them
list_of_chr = unique(comb$species)
list_of_files = list.files("/home/zajac/PhyloTree/GenomicMosaicism/", pattern = ".paths.order.length.txt", full.names = T) #files created with ippg/06_label_nodes_by_order_on_chr.R

divergence_in_inversions = NULL
for (xm in list_of_chr[i]){
  inversions = comb[comb$species == xm,c(6:8)]
  file_of_nodes = list_of_files[grepl(paste0(xm,".paths"),list_of_files)]
  if (nrow(inversions) > 0){
    file_of_nodes = fread(file_of_nodes)
  }
  
  summary_of_nodes = lapply(seq_len(nrow(inversions)), function(i) {
    x <- inversions[i, ]
    
    file_of_nodes[
      file_of_nodes$start >= x$start_species_syri &
        file_of_nodes$end <= x$end_species_syri,
    ] %>%
      group_by(cat) %>%
      summarise(length = sum(V2), .groups = "drop") %>%
      mutate(
        species = x$species,
        subject_start = x$start_species_syri,
        subject_end = x$end_species_syri
      )
  })
  
  summary_of_nodes <- dplyr::bind_rows(summary_of_nodes)
  divergence_in_inversions = rbind(divergence_in_inversions,bind_rows(summary_of_nodes)) 
}

comb = merge(comb, 
             divergence_in_inversions %>% 
               unique() %>% 
               group_by(species, start_species_syri, end_species_syri) %>% 
               mutate(sum(length)) %>% 
               ungroup() %>% 
               mutate(perc = length*100/`sum(length)`) %>% 
               dplyr::select(-length) %>% 
               unique() %>% 
               pivot_wider(names_from = cat, values_from = perc), 
             by = c("species" ,"start_species_syri", "end_species_syri"), all.x = TRUE) 
#comb[comb$species == "rPodTil1#2#8" & comb$subject_start == 78122500 & comb$subject_end == 78142455,]$private = 100
write_delim(comb, "/home/zajac/INVERSIONS/inversions.syri.impg.mapc.filtered.clustered.0.9.node_annotation.txt", delim = "\t")

##Plot
library(RIdeogram)
bed = read.delim("/groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_unmasked/rPodCre2.1#1.all.fasta.fai", header = F)
bed$Start = 1
bed = bed[,c(1,6,2)]
colnames(bed) = c("Chr", "Start", "End")
bed = bed %>% mutate(CH = str_remove(Chr, "rPodCre2.1#1#")) %>% arrange(as.numeric(CH)) %>% dplyr::select(1,2,3)

gr_merged <- GRanges(
  seqnames = comb$chr,
  ranges   = IRanges(start = comb$start, 
                     end = comb$end), 
)

whole_chr <- GRanges(
  seqnames = bed$Chr,
  ranges = IRanges(start = bed$Start, end = bed$End)
)
non_inv <- setdiff(whole_chr, gr_merged)
mcols(gr_merged)$value <- 100
mcols(non_inv)$value <- 0
final_gr <- c(gr_merged, non_inv)
final_gr <- sort(final_gr)
final_table <- data.frame(
  chr   = as.character(seqnames(final_gr)),
  start = start(final_gr),
  end   = end(final_gr),
  value = mcols(final_gr)$value
)
final_table = final_table %>% 
  left_join(comb[,c(2,4,5,14)]) %>% 
  mutate(value = case_when(
    is.na(is_singleton) ~ 0, 
    is_singleton == "FALSE" ~ 100, 
    is_singleton == "TRUE" ~ 50)) %>% 
  dplyr::select(1,2,3,4,5) 
colnames(final_table) = c("Chr", "Start", "End", "Value")
final_table = unique(final_table)
final_table <- final_table %>%
  ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE)
final_table$Chr = factor(final_table$Chr, levels = paste0("rPodCre2.1#1#", c(seq(1,18,1), "Z")))
final_table$Chr <- unname(as.character(final_table$Chr))
ideogram(karyotype = bed, overlaid = final_table, colorset1 = c("white","lightpink1", "plum4"), output = "Figures/Karyotype.invs.svg")

##With centromeres
centrs = read.delim("/groups/mpistaff/Zajac/analyses/Centromers/manually_curated_centromeres_quartet.txt")
centrs = merge(centrs, data.frame(Species = unique(centrs$Species), color = c("7FB800", "E62B9B", "7F4031", "7F4031", "299578", "7F4031","1642E1", "E62B9B", "7F4031", "7FB800", "1642E1", "299578", "ECE74E", "1642E1","E62B9B"), Shape = rep("box", 15)), by = "Species") %>% dplyr::select(1,12,2,9,10,11) %>% unique() %>% dplyr::rename("Type" = "Species", "Start" = "centromere_start_inCretensis", "End" = "centromere_end_inCretensis")  %>% mutate(Chr = paste0("rPodCre2.1#1#", Chr)) %>% arrange(Chr) %>% filter(!is.na(Start))
ideogram(karyotype = bed, overlaid = final_table, label = centrs, label_type = "marker", colorset1 = c("white","lightpink1", "plum4"), output = "Figures/Karyotype.invs2.svg")