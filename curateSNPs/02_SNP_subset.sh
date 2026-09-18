#!/bin/bash
 
#SBATCH --job-name=snp_subset
#SBATCH --ntasks=4
#SBATCH --array=0-18
#SBATCH --time=3-00:00:00
#SBATCH --mem=75G
#SBATCH --error=/home/zajac/snp_subset.%J.err
#SBATCH --output=/home/zajac/snp_subset.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global

DIR=/home/zajac
k=$SLURM_ARRAY_TASK_ID
module load R/4.4.2

Rscript 02_SNP_subset.R $k
