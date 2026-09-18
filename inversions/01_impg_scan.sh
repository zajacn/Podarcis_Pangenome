#!/bin/bash
#SBATCH --job-name=inversions
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-18
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/invs.%J.err
#SBATCH --output=/home/zajac/invs.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global

k=$SLURM_ARRAY_TASK_ID 

for i in $(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*syri.vcf);
do
  awk '$5 == "<INV>" {print}' $i | awk -F";" '{print $1}' | sed 's/END=//g' | awk '{print $1":"$2"-"$8}' > ${i}.forimpg.txt;
done

for i in $(ls /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/*syri.vcf.forimpg.txt); 
do 
  name=$(basename $i)
  cat $i | while read -r line;
  do
    impg query -p /home/zajac/pggb_allcommunities/community${k}.s5000.paf -r $line --output-format bedpe | awk '$9 == "-" && $10 == "+" {print}' | awk -v name=$name '{print $0"\t"name}' > /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/${line}.${name}.temp.out;
  done
done

cat /home/zajac/SYRI/*/*/*.temp.out > /home/zajac/INVERSIONS/impg.scan.of.syri.inversions.txt
rm /home/zajac/SYRI/*/*/*.temp.out
