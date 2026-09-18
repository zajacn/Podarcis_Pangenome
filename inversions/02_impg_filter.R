library(GenomicRanges)
library(data.table)
library(igraph)

invs = read.delim("/home/zajac/INVERSIONS/temp/impg.scan.of.syri.inversions.txt", header = F)
invs$V13 = str_remove(str_remove(sapply(str_split(invs$V13, "_"), .subset, 2), "syri.vcf.forimpg.txt"), ".rev") #fix names
invs = invs[invs$V1 == invs$V13,] # select only target species from the alignment
invs$leninsp = invs$V3 - invs$V2 # measure total bp in alignment in target species
invs$lenincre = invs$V6 - invs$V5 # measure total bp in alignment in P. cretensis
invs = invs %>% mutate(V8 = sapply(str_split(V7, ":"), .subset, 2), chr = sapply(str_split(V7, ":"), .subset, 1)) %>% 
  separate(V8, into = c("start", "end"), sep = "-") %>% 
  mutate(originallen = as.numeric(end) - as.numeric(start)) 
filteredinvs = invs %>% 
  group_by(V1, V7) %>% 
  dplyr::summarise(sum(leninsp), sum(lenincre), unique(originallen)) %>% 
  filter(`unique(originallen)` > 1000) %>% #filter inversions longer than 1000bp
  filter(`sum(lenincre)`/`unique(originallen)` > 0.5) %>% #filter inversions where total bp in wfmash alignment of P. cretensis sequence match the syri inversion length in at least 70%
  filter(`sum(leninsp)`/`sum(lenincre)` > 0.7) %>% #filter inversions where total bp in wfmash alignment of target species sequence match the length of P.cretensis alignment within each inversion at least 70%
  pull(V7)  %>% unique() # pull list of inversions
invs = invs[invs$V7 %in% filteredinvs,] %>% # filter inversions
  group_by(V1,V7) %>% 
  dplyr::summarise(min(V2), max(V3), unique(start), unique(end), unique(chr)) %>% #get unique coordinates
  dplyr::select(7,5,6,1,3,4)
colnames(invs) = c("chr", "start", "end", "species", "start_species_impg", "end_species_impg")
invs$start = as.numeric(invs$start)
invs$end = as.numeric(invs$end)
syri=read.delim("/home/zajac/INVERSIONS/inversions.syri.txt")
syri = merge(syri, invs, by = c("chr", "start", "end", "species"))
write_delim(syri, "/home/zajac/INVERSIONS/inversions.syri.impg.filtered.txt", delim = "\t")