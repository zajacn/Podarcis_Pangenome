#!/bin/bash
#SBATCH --job-name=dispnodesharing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-18
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/dispnodesharing.%J.err
#SBATCH --output=/home/zajac/dispnodesharing.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new

DIR=/home/zajac
k=$SLURM_ARRAY_TASK_ID
module load python/3.9.13

python 00_pipeline/pangenome_dispensable_windows.py -i community${k}.s5000/community${k}.sorted.gfa -w 50000 -o GenomicMosaicism/Disp.node.sharing.combined.50kb.community${k}.colinear.csv -t 8 --colinear
