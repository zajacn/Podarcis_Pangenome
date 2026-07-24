#!/bin/bash

#SBATCH --job-name=minigraph
#SBATCH --ntasks=10
#SBATCH --array=0-18
#SBATCH --time=2-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/minigraph.%J.err
#SBATCH --output=/home/zajac/minigraph.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new

DIR=/home/zajac/minigraph
k=$SLURM_ARRAY_TASK_ID

name="rPodCre2.1#1"
i=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/*.fasta | grep ${name})
FILES=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/*.fasta | grep -v ${name} | grep -v "combined")
ALL_FILES=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/*.fasta | grep -v "combined" | tail -n8)
MINIGRAPH=/home/zajac/scripts/Software/minigraph-0.21

#Map sequences back to the graph
for sample in $ALL_FILES;
do
  sname=$(echo $sample | sed 's"/home/zajac/chromosome_sets/community'${k}'_unmasked/""g')
  echo $sname
  echo community${k}
  $MINIGRAPH/minigraph  -t 32 -cxasm --call -j 0.1 $DIR/community${k}/community${k}.graph.out  ${sample} > $DIR/community${k}/${sname}.bed 2> $DIR/community${k}/${sname}.bed.log
  $MINIGRAPH/minigraph  -t 10 -cxasm --cov  -j 0.1  $DIR/community${k}/community${k}.graph.out ${sample} > $DIR/community${k}/${sname}.cov.out 2> $DIR/community${k}/${sname}.cov.log
  cat $DIR/community${k}/${sname}.cov.out | awk '$1 == "S" {print}' | awk '{print $2"\t"$4"\t"$5"\t"$8}' > $DIR/community${k}/${sname}.cov.txt;
done

