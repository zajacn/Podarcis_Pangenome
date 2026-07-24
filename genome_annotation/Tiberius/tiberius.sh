#!/bin/bash
#SBATCH --job-name=tiberius-gpu
#SBATCH -t 10:00:00
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --array=24
#SBATCH --mail-type=all
#SBATCH -c 40
#SBATCH -p scc-gpu 	             # the partition
#SBATCH -G A100:2                    # For requesting 2 GPUs.

module load miniforge3
module load apptainer
module load gcc/14.2.0
module load cuda

# Print out some info.
echo "Submitting job with sbatch from directory: ${SLURM_SUBMIT_DIR}"
echo "Home directory: ${HOME}"
echo "Working directory: $PWD"
echo "Current node: ${SLURM_NODELIST}"

k=$SLURM_ARRAY_TASK_ID
name=`sed -n ${k}p < /user/natalia.zajac01/u17422/Tiberius/list_of_genomes`


apptainer run --nv /user/natalia.zajac01/u17422/Software/tiberius/tiberius.sif tiberius.py \
	--genome /user/natalia.zajac01/u17422/Tiberius/${name}.fna \
	--out /user/natalia.zajac01/u17422/Tiberius/ALL_corrected_masked_annotated/${name}.v.gtf \
	--codingseq /user/natalia.zajac01/u17422/Tiberius/ALL_corrected_masked_annotated/${name}.v.all.codingseq \
	--protseq /user/natalia.zajac01/u17422/Tiberius/ALL_corrected_masked_annotated/${name}.v.aa \
	--batch_size 16 \
	--model /user/natalia.zajac01/u17422/Software/tiberius/vertebrates_weights \
	--no_softmasking