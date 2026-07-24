#!/bin/bash
 
#SBATCH --job-name=pankmer
#SBATCH --ntasks=20
#SBATCH --time=1-00:00:00
#SBATCH --mem=200G
#SBATCH --error=/home/zajac/pankmer.%J.err
#SBATCH --output=/home/zajac/pankmer.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new
 
DIR=/home/zajac/CompareHaplotypes

source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate pankmer

pankmer index -g $DIR/ANI_PanKmer/genomes/ -o $DIR/ANI_PanKmer/genomes_index.tar --threads 20
pankmer adj-matrix -i $DIR/ANI_PanKmer/genomes_index.tar -o $DIR/ANI_PanKmer/genomes_adj_matrix.csv