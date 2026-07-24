#!/bin/bash
 
#SBATCH --job-name=repeatmodeler
#SBATCH --ntasks=15
#SBATCH --array=1-15
#SBATCH --nodes=3
#SBATCH --time=6-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/repeatmodeler.%J.err
#SBATCH --output=/home/zajac/repeatmodeler.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

DIR=/home/zajac/Species_Specific_Masking
module load python/3.9.13
module load perl/5.30.1
export PATH=$PATH:/data/biosoftware/ncbi-blast/ncbi-blast-2.16.0+/bin/

k=$SLURM_ARRAY_TASK_ID
name=`sed -n ${k}p < /home/zajac/busco_list.txt`

#Run Repeatmodeler
singularity exec --bind /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ --cleanenv /home/zajac/scripts/Software/dfam-tetools-latest.sif BuildDatabase -name ${name}.db /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_unmasked/${name}
singularity exec --bind /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ --cleanenv /home/zajac/scripts/Software/dfam-tetools-latest.sif RepeatModeler -threads 15 -database ${name}.db -LTRStruct >& $DIR/${name}.db.masked.out
