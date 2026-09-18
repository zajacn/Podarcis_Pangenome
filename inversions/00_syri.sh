#!/bin/bash
#SBATCH --job-name=syri
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-18
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/syri.%J.err
#SBATCH --output=/home/zajac/syri.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate syri
k=$SLURM_ARRAY_TASK_ID 

#Run syri using bam files
dir=/home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all
fastadir=/home/zajac/chromosome_sets/community${k}_unmasked/
for i in $(ls $dir/*fasta.sorted.filtered.bam);
do
  A=$(ls $fastadir/rPodCre2.1#1*.fasta)
  sp1=$(echo $A | sed 's"'${fastadir}'/""g' | sed 's/.fasta//g')
  sp2=$(echo $i | sed 's"/home/zajac/SYRI/community'${k}'/community'${k}'.rPodCre2.1_to_all/""g' | sed 's/.fasta.sorted.filtered.bam//g')
  B=$(ls $fastadir/$sp2.fasta)
  syri -c $i -r $A -q $B -F B --nc 10 --dir $dir --prefix ${sp1}_${sp2} --samplename ${sp2}.;
done

#Select inversions
for i in $(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*syri.vcf);
do
  awk '$5 == "<INV>" {print}' $i | awk -F";" '{print $1,$2,$3,$4}' | sed 's/END=//g' | sed 's/ChrB=//g' | sed 's/StartB=//g' | sed 's/EndB=//g' > ${i}.temp.invs;
done

cat /home/zajac/SYRI/*/*temp.invs > /home/zajac/INVERSIONS/inversions.syri.txt
rm /home/zajac/SYRI/*/*temp.invs