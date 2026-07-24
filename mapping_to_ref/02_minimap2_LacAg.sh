#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-8,10-18
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --error=/home/zajac/minimap2.%J.err
#SBATCH --output=/home/zajac/minimap2.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

k=$SLURM_ARRAY_TASK_ID 

## Map LacAg to rPodCre2.1
/data/biosoftware/minimap2/minimap2/minimap2 -ax asm10 -t 12 $DIR/chromosome_sets/community${k}_unmasked/rPodCre2.1#1*.fasta $DIR/chromosome_sets/community${k}_unmasked/LacAg#1*.fasta > $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sam
/data/biosoftware/samtools/samtools-1.20/samtools sort -m4G -@12 -o $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sorted.bam $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sam
/data/biosoftware/samtools/samtools-1.20/samtools index $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sorted.bam
