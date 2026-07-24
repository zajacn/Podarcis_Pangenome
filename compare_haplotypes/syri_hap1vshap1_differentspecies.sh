#!/bin/bash
#SBATCH --job-name=syri
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=199-336
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --error=/home/zajac/syri.%J.err
#SBATCH --output=/home/zajac/syri.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew


k=$SLURM_ARRAY_TASK_ID 

#Run minimap2
line=`sed -n ${k}p < /home/zajac/CompareHaplotypes/comparisons.txt`

genome2=$(echo $line | cut -d" " -f1)
genome2_name=$(echo $genome2 | cut -d"/" -f6 | sed 's/.fasta//g')
genome1=$(echo $line | cut -d" " -f2)
genome1_name=$(echo $genome1 | cut -d"/" -f6 | sed 's/.fasta//g')
community=$(echo $line | cut -d"/" -f5 | sed 's/_unmasked//g')

mkdir /home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/$community
/data/biosoftware/minimap2/minimap2/minimap2 -ax asm5 -t 12 --eqx --cs -r2k $genome1  $genome2 > /home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/$community/${genome1_name}_${genome2_name}.sam
/data/biosoftware/samtools/samtools-1.20/samtools sort -m4G -@12 -o /home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/$community/${genome1_name}_${genome2_name}.sorted.bam  /home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/$community/${genome1_name}_${genome2_name}.sam
/data/biosoftware/samtools/samtools-1.20/samtools index /home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/$community/${genome1_name}_${genome2_name}.sorted.bam

#Run syri
source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate syri
dir=/home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/community${k}/
fastadir=/home/zajac/chromosome_sets/community${k}_unmasked/
for i in $(ls $dir/*bam);
do
  sp1=$(echo $i | sed 's"/home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/community'${k}'//""g' | sed 's/.sorted.bam//g' | cut -d"_" -f1)
  A=$(ls $fastadir/$sp1*.fasta)
  sp2=$(echo $i | sed 's"/home/zajac/CompareHaplotypes/SYRI_Hap1vsHap1_DifferentSpecies/community'${k}'/""g' | sed 's/.sorted.bam//g' | cut -d"_" -f2)
  B=$(ls $fastadir/$sp2*.fasta)
  syri -c $i -r $A -q $B -F B --nc 10 --prefix ${sp1}_${sp2}. --samplename ${sp2};
done
