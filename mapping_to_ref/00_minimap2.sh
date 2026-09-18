#!/bin/bash
 
#SBATCH --job-name=minimap2
#SBATCH --ntasks=12
#SBATCH --array=0-18
#SBATCH --time=2-00:00:00
#SBATCH --mem=30G
#SBATCH --error=/home/zajac/minimap2.%J.err
#SBATCH --output=/home/zajac/minimap2.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=standard
 
DIR=/home/zajac
k=$SLURM_ARRAY_TASK_ID

##Map per chromosome
for i in $(ls $DIR/chromosome_sets/community${k}_unmasked/*.fasta | grep -v "combined" | grep -v "rPodCre2.1" | sed 's"/home/zajac/chromosome_sets/community'$k'_unmasked/""g');
do
  /data/biosoftware/minimap2/minimap2/minimap2 -ax asm10 -t 12 --eqx --cs -r2k $DIR/chromosome_sets/community${k}_unmasked/rPodCre2.1#1*.fasta $DIR/chromosome_sets/community${k}_unmasked/$i > $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/$i.sam
  /data/biosoftware/samtools/samtools-1.20/samtools sort -m4G -@12 -o $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/$i.sorted.bam $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/$i.sam
  /data/biosoftware/samtools/samtools-1.20/samtools view -q 30 --threads 4 -b $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/$i.sorted.bam > $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/$i.sorted.filtered.bam
  /data/biosoftware/samtools/samtools-1.20/samtools index $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/$i.sorted.filtered.bam;
done

##Coverage
for i in $(ls $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*sorted.filtered.bam);
do
  echo $i >> $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/files;
  /data/biosoftware/samtools/samtools-1.21/samtools coverage $i >> $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/coverage.temp.txt;
done
paste $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/files <(grep -v "#rname" $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/coverage.temp.txt)  > $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/coverage.summary.txt
rm $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/files $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/coverage.temp.txt

#Depth
for i in $(ls $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*sorted.filtered.bam); 
do 
  /data/biosoftware/samtools/samtools-1.21/samtools depth -aa $i > ${i}.depth; 
done

#Quality - just as a sanity check
for i in $(ls $DIR/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*sorted.filtered.bam); 
do
  /data/biosoftware/samtools/samtools-1.20/samtools mpileup $i --output-extra MAPQ --fasta-ref $DIR/chromosome_sets/community${k}_unmasked/rPodCre2.1#1*.fasta |  awk '{print $1"\t"$2"\t"$7}' > ${i}.mapq;
done
