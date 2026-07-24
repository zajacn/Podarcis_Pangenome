#!/bin/bash
 
#SBATCH --job-name=snpeff
#SBATCH --ntasks=20
#SBATCH --array=0-18
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=15G
#SBATCH --error=/home/zajac/snpeff.%J.err
#SBATCH --output=/home/zajac/snpeff.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global

#Define params
k=$SLURM_ARRAY_TASK_ID
DIR=/home/zajac/pggb_allcommunities/community${k}.s5000
module load java/x64/24u2

#large SVs
for i in $(ls /home/zajac/SV_intersects/largeSVs/merged.DEFAULT.*.vcf);
do
  java -jar /home/zajac/scripts/Software/SnpEff/snpEff/snpEff.jar ann -c /data/biosoftware/SnpEff/snpEff/snpEff.config rPodCre2.1#1 $i > ${i}.annot.vcf;
done