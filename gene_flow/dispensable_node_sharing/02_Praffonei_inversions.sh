#!/bin/bash
 
#SBATCH --job-name=impgtogfainv
#SBATCH --ntasks=25
#SBATCH --array=3,0,11,2,17,1,15,12,7,5,8,18,4
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --error=/home/zajac/impgtogfainv.%J.err
#SBATCH --output=/home/zajac/impgtogfainv.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

k=$SLURM_ARRAY_TASK_ID
DIR=/home/zajac/INVERSIONS/community${k}
community=community${k}
pggb=/home/zajac/scripts/Software/pggb/pggb_updated.sif
source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate wfmash
module load python/3.9.13

cat $DIR/inversions.raffonei.txt | while read -r line; 
do 
  samtools faidx $DIR/impg.${line}.fasta 
  wfmash -t 15 $DIR/impg.${line}.fasta -p 90 > $DIR/impg.${line}.paf
  f=$DIR/impg.${line}.fasta
  newname1=$(echo "$f" | tr ':' '_' | tr "#" "_" | tr "-" "_")
  cp "$f" "$newname1"
  f=$DIR/impg.${line}.paf
  newname2=$(echo "$f" | tr ':' '_' | tr "#" "_" | tr "-" "_")
  cp "$f" "$newname2"
  singularity run --bind $DIR $pggb seqwish -p $newname2 --seqs=$newname1 --gfa=${newname2}.gfa -k 0 # seqwish does not tolerate # in the name
  rm $newname2 $newname1
  l1=$(echo $line | awk -F ":" '{print $2}' | awk -F"-" '{print $2}')
  l2=$(echo $line | awk -F ":" '{print $2}' | awk -F"-" '{print $1}')
  diff=$((l1-l2))
  python 00_pipeline/pangenome_dispensable_windows.py -i ${newname2}.gfa -w $diff -o ${newname2}.dispnodes.colinear.csv -t 25 --colinear; # use the length of the inversion as a single window
done
