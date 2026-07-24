#!/bin/bash

#SBATCH --job-name=minigraph
#SBATCH --ntasks=10
#SBATCH --array=0-18
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/minigraph.%J.err
#SBATCH --output=/home/zajac/minigraph.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new

#see https://github.com/AnimalGenomicsETH/superpangenome_construction/blob/main/snakepit/construct_pangenomes.smk
#Needs special version of mgutils - commit 38f04593f9c9ef8b1085481d0b50040bec83de89

DIR=/home/zajac/minigraph
k=$SLURM_ARRAY_TASK_ID
RefGraph=$DIR/community${k}/community${k}.graph.out
OUT_DIR=$(dirname $RefGraph)
REF=/home/zajac/chromosome_sets/community${k}_unmasked/rPodCre2.1*.fasta
BCFTOOLS_DIR=/data/biosoftware/bcftools/bcftools-1.21
VCF=$OUT_DIR/community${k}.graph.waved.vcf
SAMPLES=$(grep -v "##" $VCF | head -n1 | sed 's/\t/\n/g' | tail -n14)

##Normalize and filter for large SVs/ splity by sample
for SAMPLE in $SAMPLES;
do
##Select one sample
  BASE="$OUT_DIR/community${k}.${SAMPLE}.graph.waved.norm"
  $BCFTOOLS_DIR/bcftools view -s ${SAMPLE} $VCF \
  | $BCFTOOLS_DIR/bcftools norm -f ${REF} -c s -m -any --threads 12 \
  | $BCFTOOLS_DIR/bcftools view -f 'PASS,.' -e 'GT="." | GT="0" | GT="0|0" | GT=".|." | GT= "0|." | GT= ".|0"' --threads 12 \
  | $BCFTOOLS_DIR/bcftools sort -m 100G -T $OUT_DIR/bcftools-sort.XXXXXX \
  | $BCFTOOLS_DIR/bcftools annotate -x ^INFO/LEN,^INFO/TYPE,^FORMAT/GT  \
  | $BCFTOOLS_DIR/bcftools norm -d exact -Oz --threads 12 -o $BASE.vcf.gz \
  && $BCFTOOLS_DIR/bcftools index -t $BASE.vcf.gz \
  && $BCFTOOLS_DIR/bcftools view -i "STRLEN(REF)>49 | STRLEN(ALT)>49" --threads 12 -Ov -o $BASE.SV.vcf $BASE.vcf.gz;
done