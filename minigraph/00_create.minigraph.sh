#!/bin/bash

#SBATCH --job-name=minigraph
#SBATCH --ntasks=32
#SBATCH --array=0-18
#SBATCH --time=6-00:00:00
#SBATCH --mem=250G
#SBATCH --error=/home/zajac/minigraph.%J.err
#SBATCH --output=/home/zajac/minigraph.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new

DIR=/home/zajac/minigraph
k=$SLURM_ARRAY_TASK_ID

DIRECTORY=$DIR/community${k}
if [ -d "$DIRECTORY" ]; then
  mkdir $DIRECTORY 
fi

name="rPodCre2.1#1"
i=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/*.fasta | grep ${name})
FILES=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/*.fasta | grep -v ${name} | grep -v "combined")
ALL_FILES=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/*.fasta | grep -v "combined")

#Create a graph
echo $FILES
/home/zajac/scripts/Software/minigraph-0.21/minigraph -t 32 -j 0.1 -cxggs $i $FILES > $DIR/community${k}/community${k}.graph.out
grep 'SN:Z:'  $DIR/community${k}/community${k}.graph.out | awk '{print $4"\t"$5}' > $DIR/community${k}/contributing.genome



