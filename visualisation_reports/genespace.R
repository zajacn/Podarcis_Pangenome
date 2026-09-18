library(GENESPACE)
genomeRepo <- "../GeneSpace2/results_one_hap/"
wd <- "../GeneSpace2/results_one_hap/"
path2mcscanx <- "/home/zajac/scripts/Software/MCScanX-1.0.0/"
gp_par = init_genespace(
  wd,
  genomeIDs = str_remove(list.files("../GeneSpace2/results_one_hap/bed"), ".bed"),
  ploidy = 1,
  ignoreTheseGenomes = NULL,
  path2orthofinder = "/home/zajac/.conda/envs/of3_env/bin/",
  path2diamond = "/data/biosoftware/diamond/diamond",
  path2mcscanx = path2mcscanx,
  rawOrthofinderDir = "../GeneSpace2/results_one_hap/orthofinder/",
  nCores=10
)
gpar <- run_genespace(gsParam = gp_par) 