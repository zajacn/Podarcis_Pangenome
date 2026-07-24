#!/bin/bash
 
#SBATCH --job-name=vg
#SBATCH --ntasks=20
#SBATCH --array=0-18
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=70G
#SBATCH --error=/home/zajac/vg.%J.err
#SBATCH --output=/home/zajac/vg.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global


k=$SLURM_ARRAY_TASK_ID
DIR=/home/zajac/pggb_allcommunities
BCFTOOLS_DIR=/data/biosoftware/bcftools/bcftools-1.21
ref=$(singularity run /home/zajac/scripts/Software/pggb/pggb_updated.sif odgi paths -L -i $DIR/community${k}.s5000/community${k}.sorted.og | grep 'rPodCre2.1')
REF=$DIR/../chromosome_sets/community${k}_unmasked/rPodCre2.1*.fasta
INX=$DIR/../chromosome_sets/community${k}_unmasked/rPodCre2.1*.fai

##Run vg deconstruct according to - https://github.com/pangenome/pggb-paper/blob/main/scripts/gfa2evaluation.sh

singularity run /home/zajac/scripts/Software/pggb/pggb_updated.sif vg deconstruct -P $ref -a -t 20 $DIR/community${k}.s5000/community${k}.sorted.gfa > $DIR/community${k}.s5000/community${k}.sorted.gfa.vcf

##Pop bubbles
singularity run /home/zajac/scripts/Software/pggb/pggb_updated.sif vcfbub -l 0 -a 100000 --input $DIR/community${k}.s5000/community${k}.sorted.gfa.vcf >  $DIR/community${k}.s5000/community${k}.sorted.gfa.bub.vcf #allele length no greater than 100k
bgzip $DIR/community${k}.s5000/community${k}.sorted.gfa.bub.vcf
tabix $DIR/community${k}.s5000/community${k}.sorted.gfa.bub.vcf.gz

##Split the vcf into chuncks
/data/biosoftware/bedtools2/bedtools2/bin/bedtools makewindows -b <(awk '{print $1"\t"1"\t"$2}' $INX) -w 1000000 | awk '{print $1":"$2"-"$3"\t""block_"$2}' > $DIR/community${k}.s5000/chuncks.txt
cat $DIR/community${k}.s5000/chuncks.txt | while read -r line;
do
  region=$(echo $line | cut -d" " -f1)
  prefix=$(echo $line | cut -d" " -f2)
  $BCFTOOLS_DIR/bcftools view -r $region --output-file /home/zajac/pggb_allcommunities/community${k}.s5000/${prefix}.vcf /home/zajac/pggb_allcommunities/community${k}.s5000/community${k}.sorted.gfa.bub.vcf.gz;
done

##Vcfwave each chnuck separately
for file in /home/zajac/pggb_allcommunities/community${k}.s5000/block_*.vcf;
do
  sbatch --job-name=vcfwave --partition=highmem2new --time=2-00:00:00 --mem=150G --ntasks=15 --error=/home/zajac/vcfwave.%J.err --output=/home/zajac/vcfwave.%J.out --wrap="echo $file; singularity run /home/zajac/scripts/Software/pggb/pggb_updated.sif vcfwave -I 1000 -t 15  $file  >  $file.waved.vcf";
done

##Check
for i in $(ls *waved.vcf | sed 's/.waved.vcf//g'); 
do 
  paste <(tail -n1 $i | awk '{print $1,$2}') <(tail -n1 $i.waved.vcf | awk '{print $1,$2}') | awk '{print $4-$2}' | awk '{$1 < 0}'; 
done

##Concatenate
VCFWAVE=$(awk '{print $2".vcf.waved.vcf"}' $DIR/community${k}.s5000/chuncks.txt)
for i in $VCFWAVE;
do
  $BCFTOOLS_DIR/bcftools sort --output  $DIR/community${k}.s5000/${i}.sorted.vcf.gz --write-index --output-type z $DIR/community${k}.s5000/${i};
done
VCFWAVE=$(awk -v dir=$DIR/community${k}.s5000/ '{print dir$2".vcf.waved.vcf.sorted.vcf.gz"}' $DIR/community${k}.s5000/chuncks.txt)
$BCFTOOLS_DIR/bcftools concat $VCFWAVE --allow-overlaps --output $DIR/community${k}.s5000/community${k}.pggb.waved.vcf




