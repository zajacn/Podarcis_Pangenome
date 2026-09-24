#Check overlap between SVs and TEs
library(vcfR)
library(utils)
library(dplyr)
library(tidyverse)
library(reshape2)
library(cowplot)
library(patchwork)
library(ezRun)
library(ggbeeswarm)
library(data.table)
library(DT)

# Metadata
sample_data = read.delim("/groups/mpistaff/Zajac/genomes/sample_data.tsv")
sample_data$clade=factor(sample_data$clade,levels=c("Balkan","Sicilian-Maltese","Western","Siculus","Muralis","Iberian"))


# Bed coordinates
bed = read.delim("/groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_unmasked/rPodCre2.1#1.all.fasta.fai", header = F)
bed$Start = 1
bed = bed[,c(1,6,2)]
colnames(bed) = c("Chr", "Start", "End")

# Repeats and TEs
repeats = fread("../Species_Specific_Masking/repeats.allhaplotypes.tsv", sep = "\t")

# SVs (INDELs)
all_comb = NULL
for (i in list.files("../SV_intersects/largeSVs/", pattern = "merged.DEFAULT")[!grepl("annot",list.files("../SV_intersects/largeSVs/", pattern = "merged.DEFAULT"))]){
    comb = read.vcfR(paste0("../SV_intersects/largeSVs/", i))
    comb = cbind(data.frame(comb@fix), data.frame(comb@gt))
    comb = comb %>% 
      mutate(countALT = nchar(ALT), countN = str_count(ALT, "N")) %>% 
      mutate(perc = countN/countALT) %>% filter(perc < 0.5) 
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
    comb = comb[,c(1,2,3,10:23,27:33)]
    comb = comb %>% dplyr::mutate(across(colnames(comb)[4:17], ~ case_when(. == "0|0" ~ 0, . == "0" ~ 0, . == "1|1" ~ 1, . == "1" ~ 1,  . == "0|1" ~ 1, . == "1|0" ~ 1, . == ".|1" ~ 1, . == "0|."~ 0, . == ".|0" ~ 0, . == "1|." ~ 1 )))
    all_comb[[i]] = comb
}
all_comb = bind_rows(all_comb)
all_comb$LEN = abs(all_comb$ALTLEN - all_comb$REFLEN)

# Annotation of DELETIONs for TEs and Repeats 

cretensis_repeats = makeGRangesFromDataFrame(repeats[repeats$species == "rPodCre2.1",c(1,3,4,5,6,7)] %>% 
                                             seqnames.field = "chromosome", 
                                             start.field = "positioninquery_begin", 
                                             end.field = "positioninquery_end", 
                                             keep.extra.columns = TRUE)

deletion = makeGRangesFromDataFrame(all_comb %>% 
                                       filter(SVTYPE == "DEL") %>% 
                                       dplyr::select(1,2,4:17,25) %>% 
                                       mutate(END = as.numeric(POS)+as.numeric(LEN)) %>% separate(CHROM, into = c("species", "haplotype", "chrom"), sep = "#") %>% dplyr::select(1,3,4:18,20), seqnames.field = "chrom", start.field = "POS", end.field = "END", keep.extra.columns = T)

hits = findOverlaps(deletion, cretensis_repeats)
query_hits   <- deletion[queryHits(hits)]
subject_hits <- cretensis_repeats[subjectHits(hits)]
overlaps <- pintersect(query_hits, subject_hits)
result <- data.frame(
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
result$original_len = result$query_end - result$query_start
result$length_of_TE = result$subject_end - result$subject_start
result = setdiff(result, result[result$overlap_length/result$length_of_TE < 0.1 & result$overlap_length/result$original_len < 0.1 & result$overlap_length < 200,]) #filtering criterium
deletions = merge(data.frame(deletion)[,c(1:4)], 
                 result, 
                 by.x = c("seqnames", "start", "end"), 
                 by.y = c("query_seqname", "query_start", "query_end"), all.x = TRUE) %>% 
  mutate(overlap_length = if_else(is.na(overlap_length), 0, as.numeric(overlap_length))) %>% 
  group_by(seqnames, start, end, width) %>% 
  slice_max(overlap_length)
deletions = deletions %>% mutate(Category = 
                                 case_when(repeats %in% c("LTR","DNA","SINE","LINE", "PLE") ~ "Transposable Elements", 
                                           repeats %in% c("Satellite","Simple_repeat","Low_complexity","snRNA","tRNA","rRNA", "Unknown") ~ "Repeats", 
                                           is.na(repeats) ~ "Not TE or Repeat")) 
write_delim(deletions, "../Species_Specific_Masking/transposable_elements_fasta/temp.DEL.filt.txt", delim  = "\t")
# Annotation of INSERTIONs for TEs and Repeats - performed with blast
TESINS = list.files("/home/zajac/Species_Specific_Masking/transposable_elements_fasta", pattern = "INS.", full.names = TRUE)[-1]
insertions=NULL
for (i in TESINS[-4]){
  df = fread(i, header = F)
  df$subject_len = abs(df$V10 - df$V9) # calculate legnth of the transposable element overlap
  df$query_len = abs(df$V8 - df$V7) # calculate length of the insertion overlap
  df = df[,c(1:4,13:14)]
  setDT(df)
  df[, c("X", "positioninquery_end") := tstrsplit(V2, "-", fixed = TRUE)] # get TE coordinates
  df[, c("seqid", "positioninquery_begin") := tstrsplit(X, ":", fixed = TRUE)] # get TE coordinates
  df[, c("species", "haplotype", "chromosome") := tstrsplit(seqid, "#", fixed = TRUE)] 
  df[, c("file", "POS") := tstrsplit(V1, "_", fixed = TRUE)] # get insertion coordinates
  df[, c("file", "CHROM") := tstrsplit(file, ":", fixed = TRUE)] # get insertion coordinates
  df[, X := NULL]
  df = df %>% filter(as.character(chromosome) == as.character(sapply(str_split(CHROM, "#"), .subset, 3))) # filter if mapped within the same chromosome
  coln = colnames(all_comb)[grepl(unique(df$species), colnames(all_comb))][1]

  df = merge(
    merge(df, 
          all_comb[all_comb[[coln]] == 1 & all_comb$SVTYPE == "INS",c(1,2,25)], 
          by = c("CHROM", "POS"), 
          all.y = TRUE) %>% 
      filter(query_len/LEN < 3), # filter if any overlap between INS and TE is 3x greater than the length of the TE (allows for filtering matches for different INS at same positions)
    all_comb[all_comb[[coln]] == 1 & all_comb$SVTYPE == "INS",c(1,2,25)], 
    by = c("CHROM", "POS", "LEN"), 
    all.y = TRUE)
  
  df$species = if_else(is.na(df$species), sapply(str_split(coln, "_"), .subset, 3), df$species)
  df$chromosome = if_else(is.na(df$chromosome), sapply(str_split(df$CHROM, "#"), .subset, 3), df$chromosome)
  
  insertions = rbind(insertions, df)
}

repeats = repeats %>% mutate(haplotype = as.character(haplotype), chromosome = as.character(chromosome), positioninquery_begin = as.character(positioninquery_begin), positioninquery_end = as.character(positioninquery_end))
setDT(insertions)
setDT(repeats)

insertions[
  repeats,
  repeat_class := i.repeat_class,
  on = .(species, haplotype, chromosome,
         positioninquery_begin, positioninquery_end)
]
insertions = insertions[,-c("V1", "file", "V2")]
insertions = insertions %>% 
  mutate(V4 = if_else(is.na(V4),"0",as.character(V4) )) %>% 
  mutate(V3 = if_else(is.na(V3),"0",as.character(V3) ))
insertions$leninTE = as.numeric(insertions$V4)
insertions$lennotinTE = as.numeric(insertions$LEN) - as.numeric(insertions$V4)

insertions = unique(insertions)
insertions = insertions  %>% group_by(CHROM, POS, LEN, species, chromosome) %>% slice_max(as.numeric(V4))  %>% data.frame()
insertions = insertions %>% filter(lennotinTE > -5000) %>% mutate(repeat_class = if_else(repeat_class == "SINE?", "SINE", repeat_class)) 
insertions$original_TE_LEN = if_else(!is.na(insertions$positioninquery_end), as.numeric(insertions$positioninquery_end) - as.numeric(insertions$positioninquery_begin), NA)
insertions = insertions %>% mutate(Category = 
                                     case_when(repeat_class %in% c("LTR","DNA","SINE","LINE", "PLE", "RC", "Retroposon") ~ "Transposable Elements", 
                                               repeat_class %in% c("Unknown","Satellite","Simple_repeat","Low_complexity","snRNA","tRNA","rRNA") ~ "Repeats", 
                                               is.na(repeat_class) ~ "Not TE or Repeat")) 
write_delim(insertions, "../Species_Specific_Masking/transposable_elements_fasta/temp.INS.filt.txt", delim  = "\t")
