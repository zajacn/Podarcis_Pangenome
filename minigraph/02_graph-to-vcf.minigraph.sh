#!/bin/bash

#SBATCH --job-name=minigraph
#SBATCH --ntasks=10
#SBATCH --array=1
#SBATCH --time=24:00:00
#SBATCH --mem=60G
#SBATCH --error=/home/zajac/minigraph.%J.err
#SBATCH --output=/home/zajac/minigraph.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global

#see https://github.com/AnimalGenomicsETH/superpangenome_construction/blob/main/snakepit/construct_pangenomes.smk
#Needs special version of mgutils - commit 38f04593f9c9ef8b1085481d0b50040bec83de89

DIR=/home/zajac/minigraph
k=$SLURM_ARRAY_TASK_ID

RefGraph=$DIR/community${k}/community${k}.graph.out
OUT_DIR=$(dirname $RefGraph)
reference="rPodCre2.1#1"
chrom=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/* | grep -v "combined" | head -n1 | cut -d"/" -f6 | sed 's/.fasta//g' | sed 's/.rev.fasta//g' | cut -d"#" -f3)
BEDS=$(ls $OUT_DIR/*.bed | grep -v $reference)
echo $BEDS | sed 's"'$OUT_DIR/'""g' | sed 's/.fasta.bed//g' | sed 's/ /\n/g' | sed 's/.rev//g' > $OUT_DIR/samples
BCFTOOLS_DIR=/data/biosoftware/bcftools/bcftools-1.21/

##Call variants
singularity run /home/zajac/scripts/Software/pggb/pggb_latest.sif vg convert -r 0 -g $RefGraph -f > $OUT_DIR/part1.gfa
paste $BEDS | /home/zajac/scripts/Software/k8-1.2/k8-x86_64-Linux /home/zajac/scripts/Software/minigraph-0.21/misc/mgutils_AL.js path $OUT_DIR/samples - > $OUT_DIR/part2.gfa
sed -i 's/s//g' $OUT_DIR/part2.gfa
cat $OUT_DIR/part1.gfa $OUT_DIR/part2.gfa > $OUT_DIR/community${k}.graph.gfa
singularity run /home/zajac/scripts/Software/pggb/pggb_latest.sif vg deconstruct -P ${reference}#${chrom} -a $OUT_DIR/community${k}.graph.gfa > $OUT_DIR/community${k}.graph.vcf

##Pop bubbles
singularity run /home/zajac/scripts/Software/pggb/pggb_latest.sif vcfbub -l 0 -a 100000 --input $OUT_DIR/community${k}.graph.vcf > $OUT_DIR/community${k}.graph.sorted.gfa.bub.vcf #allele length no greater than 100k
bgzip $OUT_DIR/community${k}.graph.sorted.gfa.bub.vcf
tabix $OUT_DIR/community${k}.graph.sorted.gfa.bub.vcf.gz

##Split the vcf into chuncks
/data/biosoftware/bedtools2/bedtools2/bin/bedtools makewindows -b <(awk '{print $1"\t"1"\t"$2}' $DIR/../chromosome_sets/community${k}_unmasked/rPodCre2.1*.fai) -w 1000000 | awk '{print $1":"$2"-"$3"\t""block_"$2}' > $OUT_DIR/chuncks.txt
cat $OUT_DIR/chuncks.txt | while read -r line;
do
  region=$(echo $line | cut -d" " -f1)
  prefix=$(echo $line | cut -d" " -f2)
  $BCFTOOLS_DIR/bcftools view -r $region --output-file $OUT_DIR/${prefix}.vcf  $OUT_DIR/community${k}.graph.sorted.gfa.bub.vcf.gz;
done

##Vcfwave each chnuck separately
for file in  $OUT_DIR/block_*.vcf;
do
  sbatch --job-name=vcfwaveinMING --partition=standard --time=2-00:00:00 --mem=40G --ntasks=15 --error=/home/zajac/vcfwave.%J.err --output=/home/zajac/vcfwave.%J.out --wrap="echo $file; singularity run /home/zajac/scripts/Software/pggb/pggb_latest.sif vcfwave -I 1000 -t 15  $file  >  $file.waved.vcf";
done

# #Concatenate
VCFWAVE=$(awk '{print $2".vcf.waved.vcf"}'  $OUT_DIR/chuncks.txt)
for i in $VCFWAVE;
do
  $BCFTOOLS_DIR/bcftools sort --output   $OUT_DIR/${i}.sorted.vcf.gz --write-index --output-type z $OUT_DIR/${i};
done
VCFWAVE=$(awk -v dir=$OUT_DIR/ '{print dir$2".vcf.waved.vcf.sorted.vcf.gz"}' $OUT_DIR/chuncks.txt)
$BCFTOOLS_DIR/bcftools concat $VCFWAVE --allow-overlaps --output $OUT_DIR/community${k}.graph.waved.vcf

#Compare with output of gfatools
cat $OUT_DIR/community${k}.graph.waved.vcf | grep -v "##" | awk '{if(length($4) > length($5)) print $1"\t"$2"\t"($2+length($4)-1); else print $1"\t"$2"\t"($2+length($5)-1)}' > $OUT_DIR/community${k}.vgdeconstruct.bed
cat $OUT_DIR/community${k}.vgdeconstruct.bed > $OUT_DIR/community${k}.overlap_with_gfatools.txt
/home/zajac/scripts/Software/gfatools/gfatools bubble $RefGraph > $OUT_DIR/community${k}.gfatools.bed
cat $OUT_DIR/community${k}.gfatools.bed | awk '{ print $1"\t"$2"\t"$2+length($14)}' | sort -k 1,1 -k2,2n > $OUT_DIR/community${k}.gfatools.SVs.bed
/data/biosoftware/bedtools2/bedtools-2.30.0/bin/bedtools intersect -a $OUT_DIR/community${k}.vgdeconstruct.bed -b $OUT_DIR/community${k}.gfatools.SVs.bed -u | wc -l  >> $OUT_DIR/community${k}.overlap_with_gfatools.txt
