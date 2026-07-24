library(ezRun)
library(dplyr)
library(tidyverse)
library(sys)
library(vcfR)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

i=args[1]

#Define files for running IMPG
fai_file <- paste0("/home/zajac/chromosome_sets/community", i, "_unmasked/rPodCre2.1*.fai")
tmp_bed  <- paste0("/home/zajac/pggb_allcommunities/community", i, ".s5000/chroms_", i, ".bed")
out_file <- paste0("/home/zajac/pggb_allcommunities/community", i, ".s5000/chuncks.txt")
paf_file <- paste0("/home/zajac/pggb_allcommunities/community", i, ".s5000.paf")
invs_bed <- paste0("/home/zajac/pggb_allcommunities/community", i, ".s5000/invs.txt")

#Run IMPG
cmd1 <- paste0("awk '{print $1\"\\t\"1\"\\t\"$2}' ", fai_file, " > ", tmp_bed)
ezSystem(cmd1)
cmd2 <- paste0("/data/opt/bedtools2/bin/bedtools makewindows -b ", tmp_bed, " -w 10000000 > ", out_file)
ezSystem(cmd2)
cmd3 <- paste0("rm ", tmp_bed)
ezSystem(cmd3)
cmd4 = paste0("impg query -p ", paf_file, " -b ", out_file, " --output-format bedpe --merge-distance 1000 > ", invs_bed)
ezSystem(cmd4)

for (x in list.files(paste0("/home/zajac/SYRI/community",i, "/community", i, ".rPodCre2.1_to_all/"), pattern = "sorted.bam$")){
  #Define files for running SYRI
  BAM=x
  A=list.files(paste0("/home/zajac/chromosome_sets/community",i, "_unmasked/"), pattern = "rPodCre2.1")[1]
  sp1=str_remove(A, ".fasta")
  A=paste0(paste0("/home/zajac/chromosome_sets/community",i, "_unmasked/"),A)
  sp2=str_remove(str_remove(BAM, ".fasta.sorted.bam"), ".rev")
  B=paste0(paste0("/home/zajac/chromosome_sets/community",i, "_unmasked/"),list.files(paste0("/home/zajac/chromosome_sets/community",i, "_unmasked/"), pattern = sp2)[1])
  BAM=paste0(paste0("/home/zajac/SYRI/community",i, "/community", i, ".rPodCre2.1_to_all/"), BAM)
  #Run SYRI
  cmd5=paste("source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh; conda activate syri; syri -c", BAM, "-r", A, "-q", B, "-F B --nc 10 --prefix", paste0(sp1,"_",sp2, "."), "--samplename", sp2, "--dir", paste0("/home/zajac/SYRI/community",i, "/community", i, ".rPodCre2.1_to_all/"))
  ezSystem(cmd5)
}

#Merge with Jasmine
cmd6=paste0("source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh; conda activate jasmine; ls ", dirname(BAM), "/rPodCre2.1*syri.vcf > ", dirname(BAM), "/list_of_files")
ezSystem(cmd6)
cmd7=paste0("jasmine file_list=", dirname(BAM), "/list_of_files --output_genotypes --default_zero_genotype max_dist=1000 --keep_var_ids out_file=", dirname(BAM), "/community", i, ".combined.syri.vcf")
ezSystem(cmd7)

#Load IMPG files
invs = read.delim(invs_bed, header = F)
invs = invs %>% arrange(V5, V6, V1) %>% mutate(run_id = rep(seq_along(rle(V9)$lengths), times = rle(V9)$lengths) ) %>% # Identify run groups - each continuous stretch of same orientation gets a unique ID
  filter(V9 == "-") %>%
  # Group by run_id, P.cretensis chromosome and focal genome chromosome
  group_by(V1,V4, run_id) %>%
  # Get min start and max end for each continuous section
  dplyr::summarise(
    Cretensis_start = min(as.numeric(V5)),
    Cretensis_end = max(as.numeric(V6)),
    X_start = min(as.numeric(V2)),
    X_end = max(as.numeric(V3)),
    .groups = "drop"
  ) 
colnames(invs) = c("X", "Cretensis", "inv_id", "Cretensis_start", "Cretensis_end", "X_start", "X_end") 

#Load SYRI files
options(scipen = 999)
syriinv = read.vcfR(paste0("/home/zajac/SYRI/community",i, "/community", i, ".rPodCre2.1_to_all/community", i, ".combined.syri.vcf"))
syriinv = cbind(data.frame(syriinv@fix), data.frame(syriinv@gt))
syriinv$ID2 = gsub('[[:digit:]]+', '', syriinv$ID)
syriinv = syriinv[syriinv$ID2 == "INV",]
syriinv = syriinv %>% separate(INFO, into = paste0("V", seq(1,7,1)), sep = ";") %>% mutate(V3 = as.numeric(str_remove(V3,"StartB=")))  %>% mutate(V4 = as.numeric(str_remove(V4,"EndB="))) %>% mutate(V1 = as.numeric(str_remove(V1,"END=")))  
cols=colnames(syriinv)[16:(ncol(syriinv)-1)]
syriinv = syriinv %>% mutate(dplyr::across(dplyr::all_of(cols) , ~ sapply(str_split(., ":"), .subset, 1))) %>% dplyr::select(1,2,3,8,cols) %>% reshape2::melt(id.vars = c("CHROM", "POS", "V1", "ID")) %>% dplyr::filter(value != "0|0") %>% mutate(variable = sapply(str_split(variable, "_"), .subset, 2)) %>% mutate(variable = str_replace(str_replace(variable,  "\\.", "#"),  "\\.", "#")) 
chr=sapply(str_split(unique(syriinv$CHROM), "#"), .subset, 3)
syriinv = syriinv %>% mutate(variable = if_else(grepl("rPodLil1#2#1.", variable), paste0("rPodLil1.2#1#",chr), variable))
syriinv = syriinv[,c(1:5)]
colnames(syriinv) = c("Cretensis", "Cretensis_start", "Cretensis_end","ID", "X")

#Combine into a final set of inversions identified by both tools
syriinv = makeGRangesFromDataFrame(syriinv %>% unique() %>% dplyr::rename(start = Cretensis_start, end = Cretensis_end), seqnames.field = "X", keep.extra.columns=TRUE)
invs = makeGRangesFromDataFrame(invs %>% unique() %>% dplyr::rename(start = Cretensis_start, end = Cretensis_end), seqnames.field = "X", keep.extra.columns=TRUE)
hits = findOverlaps(syriinv,invs)
query_hits   <- syriinv[queryHits(hits)]
subject_hits <- invs[subjectHits(hits)]
overlaps <- pintersect(query_hits, subject_hits)
result <- data.frame(
  query_genome   = as.character(seqnames(query_hits)),
  query_start    = start(query_hits),
  query_end      = end(query_hits),
  query_ID       = query_hits$ID,
  Length_SYRI     = end(query_hits) - start(query_hits),
  subject_seqname     = subject_hits$Cretensis,
  subject_start   = start(subject_hits),
  subject_end     = end(subject_hits),
  subject_X  = as.character(seqnames(subject_hits)),
  subject_X_start = subject_hits$X_start,
  subject_X_end = subject_hits$X_end) %>% 
  unique()
write_delim(result, paste0("/home/zajac/SV_intersects/community",i, ".inversions.combined.txt"))


