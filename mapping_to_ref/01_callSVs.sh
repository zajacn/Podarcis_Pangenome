#!/bin/bash
#SBATCH --job-name=svimasm
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-18
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/svimasm.%J.err
#SBATCH --output=/home/zajac/svimasm.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new

source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate svimasm_env ##v1.0.3
k=$SLURM_ARRAY_TASK_ID 

##Run svimasm using bam files
dir=/home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all
fastadir=/home/zajac/chromosome_sets/community${k}_unmasked/

##Run in haploid mode for haploid genomes
for i in PodMur1 rPodLil1.2 rPodRaf1;
do
  svim-asm haploid $dir $dir/${i}#1#*sorted.bam $fastadir/rPodCre2.1#1*.fasta --sample ${i} --types INS,DEL
  mv $dir/variants.vcf $dir/community${k}.${i}.svimasm.norm.SV.vcf;
done

##Run in diploid mode for diploid genomes
for i in rPodBoc1 rPodErh1 rPodFil1 rPodGai1 rPodLio1 rPodMel1 rPodMur119 rPodPit1 rPodSic1 rPodTil1 rPodVau1;
do
  svim-asm diploid $dir $dir/${i}#1#*sorted.bam $dir/${i}#2#*sorted.bam $fastadir/rPodCre2.1#1*.fasta --sample ${i} --types INS,DEL
  mv $dir/variants.vcf $dir/community${k}.${i}.svimasm.norm.SV.vcf;
done
