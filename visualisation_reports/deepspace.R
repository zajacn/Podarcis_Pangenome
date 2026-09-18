
library(dplyr)
library(DEEPSPACE)
library(stringr)

options(timeout = 1000)
workHere <- "/home/zajac/DeepSpace" 
dir.create(workHere)

fastaFiles <- list.files("/home/zajac/chromosome_sets/ALL/", pattern = "all.fasta")
fastaFiles <- paste0("/home/zajac/chromosome_sets/ALL/", fastaFiles)
names(fastaFiles) <- str_remove_all(str_remove_all(fastaFiles, ".all.fasta"),"/home/zajac/chromosome_sets/ALL/")
fastaFiles = fastaFiles[c(2,15)] #fastaFiles[c(2,15,8,1,10,13,3,4,6,9,5,12,14,7,11)]
test <- clean_windows(
  faFiles = fastaFiles,
  genomeIDs = names(fastaFiles)[c(2,15)],
  wd = workHere,
  
  preset = "fast",
  
  stripChrname = ".*#1#",
  minChrLen = 2e6, 
  
  nCores = 8,
  MCScanX_hCall = "/home/zajac/scripts/Software/MCScanX-1.0.0/MCScanX_h",
  minimap2call = "/data/biosoftware/minimap2/minimap2-2.28_x64-linux/minimap2")

#sbatch --ntasks=4 --nodes=1 --time=24:00:00 --mem=50G --error=dpspc.err --output=dpspc.out --mail-user=zajac@evolbio.mpg.de --partition=standard --wrap="Rscript --vanilla /home/zajac/scripts/deepspace.R"
