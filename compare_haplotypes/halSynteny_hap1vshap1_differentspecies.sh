#!/bin/bash
#SBATCH --job-name=hal
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-18
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/hal.%J.err
#SBATCH --output=/home/zajac/hal.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new


x=$SLURM_ARRAY_TASK_ID 
DIR=/home/zajac/CompareHaplotypes/HalSynteny_Hap1vsHap1_DifferentSpecies

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodTil1_1" --targetGenome "rPodLil1.2_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodTil1_rPodLil1.2.psl
#cat community*.rPodTil1_rPodLil1.2.psl > rPodTil1_rPodLil1.2.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodPit1_1" --targetGenome "rPodLil1.2_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodPit1_rPodLil1.2.psl
#cat community*.rPodPit1_rPodLil1.2.psl > rPodPit1_rPodLil1.2.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodFil1_1" --targetGenome "rPodRaf1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodFil1_rPodRaf1.psl
#cat community*.rPodFil1_rPodRaf1.psl > rPodFil1_rPodRaf1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodMur119_1" --targetGenome "PodMur1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodMur119_PodMur1.psl
#cat community*.rPodMur119_PodMur1.psl > rPodMur119_PodMur1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodSic1_1" --targetGenome "PodMur1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodSic1_PodMur1.psl
#cat community*.rPodSic1_PodMur1.psl > rPodSic1_PodMur1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodMel1_1" --targetGenome "rPodCre2.1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodMel1_rPodCre2.1.psl
#cat community*.rPodMel1_rPodCre2.1.psl > rPodMel1_rPodCre2.1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodErh1_1" --targetGenome "rPodCre2.1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodErh1_rPodCre2.1.psl
#cat community*.rPodErh1_rPodCre2.1.psl > rPodErh1_rPodCre2.1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodGai1_1" --targetGenome "rPodCre2.1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodGai1_rPodCre2.1.psl
#cat community*.rPodGai1_rPodCre2.1.psl > rPodGai1_rPodCre2.1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodLio1_1" --targetGenome "rPodVau1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodLio1_rPodVau1.psl
#cat community*.rPodLio1_rPodVau1.psl > rPodLio1_rPodVau1.all.psl

singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "rPodBoc1_1" --targetGenome "rPodVau1_1" /home/zajac/cactus/community${x}/community${x}.hal $DIR/community${x}.rPodBoc1_rPodVau1.psl
#cat community*.rPodBoc1_rPodVau1.psll > rPodBoc1_rPodVau1.all.psl
