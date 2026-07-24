#!/bin/bash

#SBATCH --job-name=biser
#SBATCH --ntasks=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/biser.%J.err
#SBATCH --output=/home/zajac/biser.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new

source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate python3.10

# biser -o /home/zajac/biser/all_genomes.SD.out -t 32 /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_masked/*#1.all.fasta

biser -o /home/zajac/biser/rPodMel1#1.all.fasta.SD.out -t 4 /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_masked/rPodMel1#1.all.fasta --gc-heap 32G