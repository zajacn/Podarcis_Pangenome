#!/bin/bash
 
#SBATCH --job-name=vcf_split
#SBATCH --ntasks=20
#SBATCH --array=0-18
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=5G
#SBATCH --error=/home/zajac/vcf_split.%J.err
#SBATCH --output=/home/zajac/vcf_split.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global

# Inspired by https://github.com/wwliao/pangenome-utils/blob/main/preprocess_vcf.sh
# Annotation added to be able to process files downstream with Jasmine - https://github.com/mkirsche/Jasmine/wiki/Jasmine-User-Manual

#Define params
k=$SLURM_ARRAY_TASK_ID
DIR=/home/zajac/pggb_allcommunities
SAMPLES=$(grep -v "##" $DIR/community${k}.s5000/community${k}.pggb.waved.vcf | head -n1 | sed 's/\t/\n/g' | tail -n14)
BCFTOOLS_DIR=/data/biosoftware/bcftools/bcftools-1.21
VCF=$DIR/community${k}.s5000/community${k}.pggb.waved.vcf
REF=$DIR/../chromosome_sets/community${k}_unmasked/rPodCre2.1*.fasta

##Split by sample and type and add all info
for SAMPLE in $SAMPLES;
do
  BASE=$DIR/community${k}.s5000/community${k}.${SAMPLE}.pggb.waved.norm
  $BCFTOOLS_DIR/bcftools view -s ${SAMPLE} $VCF \
  | $BCFTOOLS_DIR/bcftools norm -f ${REF} -c s -m -any --threads 12 \
  | $BCFTOOLS_DIR/bcftools view -f 'PASS,.' -e 'GT="." | GT="0" | GT="0|0" | GT=".|." | GT= "0|." | GT= ".|0"' --threads 12 \
  | $BCFTOOLS_DIR/bcftools sort -m 100G -T $DIR/community${k}.s5000/bcftools-sort.XXXXXX \
  | $BCFTOOLS_DIR/bcftools annotate -x ^INFO/LEN,^INFO/TYPE,^FORMAT/GT  \
  | $BCFTOOLS_DIR/bcftools norm -d exact -Oz --threads 12 -o ${BASE}.vcf.gz \
  && $BCFTOOLS_DIR/bcftools index -t ${BASE}.vcf.gz \
  && $BCFTOOLS_DIR/bcftools view -i "STRLEN(REF)>49 | STRLEN(ALT)>49" --threads 12 -Ov -o ${BASE}.SV.vcf ${BASE}.vcf.gz \
  && $BCFTOOLS_DIR/bcftools view -i "STRLEN(REF)==1 & STRLEN(ALT)==1" --threads 12 -Ov -o ${BASE}.SNP.vcf ${BASE}.vcf.gz \
  && $BCFTOOLS_DIR/bcftools view -e "STRLEN(REF)>49 | STRLEN(ALT)>49" --threads 12 -Oz -o ${BASE}.MNP1.vcf.gz ${BASE}.vcf.gz \
  && $BCFTOOLS_DIR/bcftools index -t ${BASE}.MNP1.vcf.gz \
  && $BCFTOOLS_DIR/bcftools view -e "STRLEN(REF)==1 & STRLEN(ALT)==1" --threads 12 -Ov -o ${BASE}.MNP.vcf ${BASE}.MNP1.vcf.gz;
  rm ${BASE}.MNP1.vcf.gz;
done

##Count unique SNPs and MNPs
###Per species per chr
for i in $(ls $DIR/community${k}.s5000/community${k}.*.MNP.vcf); 
do 
  echo $i >> $DIR/community${k}.s5000/samplesmnp;  
  grep -v "##" $i | awk '{print $1"\t"$2"\t"$8}' | grep -v "#CHROM" | sort | uniq | wc -l >> $DIR/community${k}.s5000/mnpcount; 
done

for i in $(ls $DIR/community${k}.s5000/community${k}.*.SNP.vcf); 
do 
  echo $i >> $DIR/community${k}.s5000/samplessnp;  
  grep -v "##" $i | awk '{print $1"\t"$2"\t"$8}' | grep -v "#CHROM" | sort | uniq | wc -l >> $DIR/community${k}.s5000/snpcount; 
done

###Per chr
cat $DIR/community${k}.s5000/community${k}.*.MNP.vcf | grep -v "##" | awk '{print $1"\t"$2"\t"$8}' | grep -v "#CHROM" | sort | uniq | wc -l >> $DIR/community${k}.s5000/mnpcount
cat $DIR/community${k}.s5000/community${k}.*.SNP.vcf | grep -v "##" | awk '{print $1"\t"$2}' | grep -v "#CHROM" | sort | uniq | wc -l >> $DIR/community${k}.s5000/snpcount

###Combine
for i in $(ls $DIR/community${k}.s5000/community${k}.*.MNP.vcf); 
do 
  grep -v "##" $i | awk '{print $1"\t"$2}' | grep -v "#CHROM" >> $DIR/community${k}.s5000/ALLMNPS.txt; done
cat $DIR/community${k}.s5000/ALLMNPS.txt | sort | uniq -c | awk '{count[$1]++} END {for (x in count) print x, count[x]}' > $DIR/community${k}.s5000/community${k}.mnp.summary
rm $DIR/community${k}.s5000/ALLMNPS.txt
