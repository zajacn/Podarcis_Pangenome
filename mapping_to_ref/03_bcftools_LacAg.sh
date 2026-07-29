#!/bin/bash
#SBATCH --job-name=snps
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-8,10-18
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --error=/home/zajac/snps.%J.err
#SBATCH --output=/home/zajac/snps.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

k=$SLURM_ARRAY_TASK_ID 

##  Call variants
bcftools=/data/biosoftware/bcftools/bcftools-1.21/bcftools
reference=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/rPodCre2.1#1#*.fasta)
FILES=$(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg*.sorted.bam | sed 's"/home/zajac/SYRI/community'${k}'/community'${k}'.rPodCre2.1_to_all/""g' | sed 's/fasta.sorted.bam//g')

for i in $FILES;
do
  $bcftools mpileup -Ou -f $reference /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}fasta.sorted.bam > /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.pileup
  $bcftools call -mv -Oz -o /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.variants.bcftools.vcf.gz /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.pileup
  $bcftools index /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.variants.bcftools.vcf.gz
  $bcftools view -v snps -Oz -o /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${i}.snps_only.vcf.gz /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${i}.variants.bcftools.vcf.gz
  $bcftools index /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${i}.snps_only.vcf.gz;
done

## Filter for regions where LacAg actually mapped
DIR=/home/zajac
samtools=/data/biosoftware/samtools/samtools-1.23/samtools
$samtools depth $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sorted.bam > $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sorted.depth

```R
for (i in seq(0,18,1)){
    df = read.delim(paste0("../SYRI/community", i, "/community", i, ".rPodCre2.1_to_all/LacAg#1fasta.sorted.depth"), header = F)
    df = df %>% arrange(V1, as.numeric(V2)) %>%
    group_by(V1) %>%
    mutate(block = cumsum(coalesce(V2 - lag(V2) != 1, TRUE))) %>%
    group_by(V1, block) %>%
    summarise(
    start = first(V2),
    end   = last(V2),
    V3    = first(V3),
    n     = n(),
    .groups = "drop")
  write_delim(df[,c(1,3,4)], paste0("../SYRI/community", i, "/community", i, ".rPodCre2.1_to_all/LacAg#1fasta.sorted.depth.bed"), delim = "\t")
}
```
sed -i '1d' */*/*depth.bed -> remove headers

## do bcftools merge all samples but only within regions where the outgroup actually mapped and filter for biallelic snps
bcftools=/data/biosoftware/bcftools/bcftools-1.21/bcftools
$bcftools merge /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*snps_only.vcf.gz --regions-file /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/LacAg#1fasta.sorted.depth.bed --missing-to-ref -Oz --output /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.forDsuite.vcf.gz
tabix /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.forDsuite.vcf.gz
$bcftools view --threads 10 -m2 -M2 --types snps /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.forDsuite.vcf.gz > /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.combined.forDsuite.vcf
