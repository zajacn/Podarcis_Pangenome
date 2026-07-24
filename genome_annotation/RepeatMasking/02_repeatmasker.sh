#!/bin/bash
 
#SBATCH --job-name=repeatmasker
#SBATCH --ntasks=15
#SBATCH --array=1-15
#SBATCH --nodes=3
#SBATCH --time=6-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/repeatmasker.%J.err
#SBATCH --output=/home/zajac/repeatmasker.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=standard

DIR=/home/zajac/Species_Specific_Masking
module load python/3.9.13
module load perl/5.30.1
export PATH=$PATH:/data/biosoftware/ncbi-blast/ncbi-blast-2.16.0+/bin/

k=$SLURM_ARRAY_TASK_ID
name=`sed -n ${k}p < /home/zajac/busco_list.txt`

singularity exec --bind /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ --cleanenv /home/zajac/scripts/Software/dfam-tetools-latest.sif RepeatMasker -pa 15 -lib $DIR/${name}.db-families.fanoprot -xsmall -gff -s -no_is -cutoff 255 -frag 20000 /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_unmasked/${name}
singularity exec --bind /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ --cleanenv /home/zajac/scripts/Software/dfam-tetools-latest.sif rmToUCSCTables.pl -out /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_unmasked/${name}.out