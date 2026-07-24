#!/bin/bash
#SBATCH --job-name=inversions
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=5,7
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/invs.%J.err
#SBATCH --output=/home/zajac/invs.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global

source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
module load R/4.4.2
module load gcc/6.2.1
k=$SLURM_ARRAY_TASK_ID 

Rscript --vanilla /home/zajac/pggb_allcommunities/01_inversions.R ${k} 
