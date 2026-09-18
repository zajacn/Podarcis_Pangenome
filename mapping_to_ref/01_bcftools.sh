#!/bin/bash
#SBATCH --job-name=snps
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-18
#SBATCH --time=24:00:00
#SBATCH --mem=45G
#SBATCH --error=/home/zajac/snps.%J.err
#SBATCH --output=/home/zajac/snps.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=standard

module load R/4.4.2
bcftools=/data/biosoftware/bcftools/bcftools-1.21/bcftools
k=$SLURM_ARRAY_TASK_ID 
reference=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/rPodCre2.1#1#*.fasta)
FILES=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*.sorted.bam | sed 's"/home/zajac/SYRI/community'${k}'/community'${k}'.rPodCre2.1_to_all/""g' | sed 's/.fasta.sorted.filtered.bam//g')
HAPS=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*.sorted.bam | sed 's"/home/zajac/SYRI/community'${k}'/community'${k}'.rPodCre2.1_to_all/""g' | sed 's/.fasta.sorted.filtered.bam//g' | grep '#2' | awk -F"#" '{print $1}')
OTHER=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*.sorted.bam | sed 's"/home/zajac/SYRI/community'${k}'/community'${k}'.rPodCre2.1_to_all/""g' | sed 's/.fasta.sorted.filtered.bam//g' | awk -F"#" '{print $1}' | sort | uniq | grep -vE "$(IFS='|'; echo "${HAPS[*]}")")

for i in $FILES;
do
  $bcftools mpileup -Ou -f $reference /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.fasta.sorted.filtered.bam > /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.pileup
  $bcftools call -mv -Oz -o /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.variants.bcftools.vcf.gz /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.pileup
  $bcftools index /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.variants.bcftools.vcf.gz
  $bcftools view -v snps -Oz -o /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${i}.snps_only.vcf.gz /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.variants.bcftools.vcf.gz
  $bcftools index /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${i}.snps_only.vcf.gz;
done

for x in $HAPS;
do
  HAP1=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${x}#1*snps_only.vcf.gz)
  HAP2=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${x}#2*snps_only.vcf.gz)
  Rscript /home/zajac/SYRI/01_bcftools.mergehaps.R $HAP1 $HAP2 /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${x}.snps_only.vcf;
done

for x in $OTHER;
do
  FILE=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${x}#1*snps_only.vcf.gz)
  mv $FILE /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${x}.snps_only.vcf.gz
  gunzip /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${x}.snps_only.vcf.gz;
done
