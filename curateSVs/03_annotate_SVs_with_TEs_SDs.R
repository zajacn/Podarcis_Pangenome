##Check which SVs are TEs/ SDs/ other repeats
library(GenomicRanges)
library(vcfR)
library(utils)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ezRun)

##Load Variants
all_comb = NULL
for (i in list.files("/home/zajac/SV_intersects/largeSVs/", pattern = "merged.DEFAULT")[!grepl("annot",list.files("/home/zajac/SV_intersects/largeSVs/", pattern = "merged.DEFAULT"))]){
  comb = read.vcfR(paste0("/home/zajac/SV_intersects/largeSVs/", i))
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
all_comb$LEN = abs(all_comb$ALTLEN - all_comb$REFLEN)


##DELETIONS vs TEs/Repeats - deletions are all present in reference
repeats = read.delim("/home/zajac/Species_Specific_Masking/repeats.allhaplotypes.tsv", sep = "\t") #For TE/Repeat annotation see genome_annotation/RepeatMasking

cretensis_repeats = makeGRangesFromDataFrame(
  repeats[repeats$species == "rPodCre2.1",c(1,3,4,5,6,7)] %>% 
    filter(repeat_class != "Unknown"), 
  seqnames.field = "chromosome", 
  start.field = "positioninquery_begin", 
  end.field = "positioninquery_end", 
  keep.extra.columns = TRUE)

deletions = makeGRangesFromDataFrame(
  all_comb %>% 
    filter(SVTYPE == "DEL") %>% 
    dplyr::select(1,2,4:17,22) %>% 
    mutate(END = as.numeric(POS)+as.numeric(LEN)) %>% 
    separate(CHROM, into = c("species", "haplotype", "chrom"), sep = "#") %>% 
    dplyr::select(1,3,4:18,20), 
  seqnames.field = "chrom", 
  start.field = "POS", 
  end.field = "END", 
  keep.extra.columns = T)

hits = findOverlaps(deletions, cretensis_repeats)
query_hits   <- deletions[queryHits(hits)]
subject_hits <- cretensis_repeats[subjectHits(hits)]
overlaps <- pintersect(query_hits, subject_hits)
deletions_tes <- data.frame(
    query_seqname   = as.character(seqnames(query_hits)),
    query_start     = start(query_hits),
    query_end       = end(query_hits),
    subject_seqname = as.character(seqnames(subject_hits)),
    subject_start   = start(subject_hits),
    subject_end     = end(subject_hits),
    overlap_start   = start(overlaps),
    overlap_end     = end(overlaps),
    repeats = subject_hits$repeat_class,
    subrepeats = subject_hits$repeat_subclass, 
    overlap_length  = width(overlaps)
) %>% unique()

deletions_tes$original_len = deletions_tes$query_end - deletions_tes$query_start
deletions_tes = deletions_tes %>% 
  group_by(query_seqname, query_start, query_end, repeats,original_len) %>% 
  summarise(overlap_length = sum(overlap_length)) %>% 
  ungroup() %>% 
  group_by(query_seqname, query_start, query_end) %>% slice_max(overlap_length, n = 1, with_ties = FALSE) %>% 
  ungroup()
deletions_tes = merge(deletions_tes, data.frame(deletions)[,c(1,2,3)] %>% 
                        mutate(cat = "deletion"), by.x = c( "query_seqname" ,"query_start" ,"query_end"), 
                      by.y = c( "seqnames" ,"start" ,"end"), all = TRUE)
deletions_tes = deletions_tes %>% mutate(repeats = if_else(is.na(repeats), "Not repeats", repeats))
deletions_tes = deletions_tes %>% mutate(original_len = if_else(is.na(original_len), query_end - query_start, original_len))

##INSERITONS vs TEs/REPEATS - insertions are not all present in reference - need annotation first
##first run 03_annotate_INS.sh - then load output
files = list.files("/home/zajac/Species_Specific_Masking/transposable_elements_fasta", pattern = "INS.", full.names = TRUE)
insertions_tes=NULL
for (i in files){
  df = read.delim(i, header = F)
  df = df %>% separate(V2, into = c("X", "positioninquery_end"), sep = "-") %>% separate(X, into = c("seqid", "positioninquery_begin"), sep = ":") %>% separate(seqid, into = c("species", "haplotype", "chromosome"), sep = "#")
  df = merge(df, repeats[,c(1:9)], by = c("species" ,"haplotype", "chromosome", "positioninquery_begin", "positioninquery_end"), all.x = TRUE)
  df = df %>% separate(V1, into = c("X", "rest"), sep = ":") %>% separate(rest, into = c("CHROM", "POS"), sep = "_")
  insertions_tes = rbind(insertions_tes, df)
}
insertions_tes = insertions_tes[insertions_tes$V3 > 85,] %>% 
  group_by(CHROM, POS, repeat_class) %>% 
  slice_max(V4, n = 1, with_ties = FALSE) %>% ungroup()
insertions_tes = merge(insertions_tes, all_comb[all_comb$SVTYPE == "INS",c(1,2,22)], by = c("CHROM","POS"), all = TRUE)
insertions_tes = insertions_tes %>% 
  mutate(repeat_class = if_else(is.na(repeat_class), "Not repeats", repeat_class)) 

##Maybe those non-TE insertions are deletions are larger segmental duplications
sds = read.delim("/home/zajac/biser/rPodCre2.1#1.all.fasta.SD.out.elem.txt", header = F, sep = "\t") # Run with biser, see genome_annotation/SegmentalDuplications
sds = makeGRangesFromDataFrame(
  sds[,c(2,3,4,5)], 
  seqnames.field = "V2", 
  start.field = "V3", 
  end.field = "V4", 
  keep.extra.columns = TRUE)

dels = makeGRangesFromDataFrame(
  deletions_tes[deletions_tes$repeats == "Not repeats",c(1,2,3)] %>% 
    mutate(CHROM = paste0("rPodCre2.1#1#",query_seqname)),
  seqnames.field = "CHROM", 
  start.field = "query_start", 
  end.field = "query_end", 
  keep.extra.columns = T)

hits = findOverlaps(dels, sds)
query_hits   <- dels[queryHits(hits)]
subject_hits <- sds[subjectHits(hits)]
overlaps <- pintersect(query_hits, subject_hits)
dels_sds <- data.frame(
  query_seqname   = as.character(seqnames(query_hits)),
  query_start     = start(query_hits),
  query_end       = end(query_hits),
  subject_seqname = as.character(seqnames(subject_hits)),
  subject_start   = start(subject_hits),
  subject_end     = end(subject_hits),
  overlap_start   = start(overlaps),
  overlap_end     = end(overlaps),
  repeats = subject_hits$V5,
  overlap_length  = width(overlaps)
) %>% unique()

dels_sds$original_len = dels_sds$query_end - dels_sds$query_start
length(unique(paste0(dels_sds$query_seqname, "-", dels_sds$query_start, "-", dels_sds$query_end)))/length(dels)

ins = makeGRangesFromDataFrame(
  insertions_tes[insertions_tes$repeat_class == "Not repeats",c(1,2)],
  seqnames.field = "CHROM", 
  start.field = "POS", 
  end.field = "POS", 
  keep.extra.columns = T)

hits = findOverlaps(ins, sds)
query_hits   <- ins[queryHits(hits)]
subject_hits <- sds[subjectHits(hits)]
overlaps <- pintersect(query_hits, subject_hits)
ins_sds <- data.frame(
  query_seqname   = as.character(seqnames(query_hits)),
  query_start     = start(query_hits),
  query_end       = end(query_hits),
  subject_seqname = as.character(seqnames(subject_hits)),
  subject_start   = start(subject_hits),
  subject_end     = end(subject_hits),
  overlap_start   = start(overlaps),
  overlap_end     = end(overlaps),
  repeats = subject_hits$V5,
  overlap_length  = width(overlaps)
) %>% unique()

ins_sds$original_len = ins_sds$query_end - ins_sds$query_start
length(unique(paste0(ins_sds$query_seqname, "-", ins_sds$query_start, "-", ins_sds$query_end)))/length(ins)
