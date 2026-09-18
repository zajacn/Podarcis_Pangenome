#!/bin/bash
#SBATCH --job-name=inversions
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=1-464
#SBATCH --time=1-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/mapc.%J.err
#SBATCH --output=/home/zajac/mapc.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=standard

module load R/4.4.2
module load gcc/6.2.1
k=$SLURM_ARRAY_TASK_ID 

Rscript --vanilla /home/zajac/INVERSIONS/temp.script.R ${k} 
sed -i '1d' /home/zajac/INVERSIONS/temp/*out.txt > /home/zajac/INVERSIONS/temp/mapc.scan.of.syri.inversions.txt 

