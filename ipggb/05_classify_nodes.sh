#!/bin/bash
 
#SBATCH --job-name=classify_nodes
#SBATCH --ntasks=15
#SBATCH --array=3
#SBATCH --time=2-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/classify_nodes.%J.err
#SBATCH --output=/home/zajac/classify_nodes.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global
 
DIR=/home/zajac
k=$SLURM_ARRAY_TASK_ID

# Inspired by https://github.com/noecochetel/North_American_Vitis_Pangenome
# The same is done within 01_dispensable_node_sharing.py

##Extract nodes
awk '$1=="P"' $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.gfa | cut -f 2-3 > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.paths
cut -f 1 $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.paths > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.ids

##Per hap
tot_len=$(wc -l $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.ids | cut -f 1 -d ' ')
for line_nb in `seq 1 ${tot_len}`; 
do 
  name=$(sed -n "${line_nb}p" $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.ids); 
  echo "sed -n "${line_nb}p" $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.paths | cut -f 2 | tr ',' '\n' | sed 's:+:\\t+:g' | sed 's:-:\\t-:g' | sed 's:^:s:g'> $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.${name}.paths"; 
done > $DIR/pggb_allcommunities/community${k}.s5000/extract_path.sh
parallel -j 26 :::: $DIR/pggb_allcommunities/community${k}.s5000/extract_path.sh

##Per genome
genomes=$(cut -f 1 -d '#' $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.ids | sort | uniq)
for genome in $genomes; 
do 
  name=$(basename ${genome}); echo "cat $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.${name}*.paths | cut -f 1 | sort -u > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.${name}.paths"; 
done > $DIR/pggb_allcommunities/community${k}.s5000/path_concat_perGenome.sh
parallel -j 26 :::: $DIR/pggb_allcommunities/community${k}.s5000/path_concat_perGenome.sh

##Classify into core, dispensable and private
for genome in $genomes; 
do 
  name=$(basename ${genome}); cat $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.${name}.paths; 
done > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.paths
cat $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.paths | sort | uniq -c > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.freq

max=15
awk -v max=$max '$1 == max {print $2}' $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.freq > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.core
awk '$1 == 1 {print $2}' $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.freq > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.private
awk '$1 < 15 && $1 > 1 {print $2}' $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.freq > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.dispensable

##Check lengths with community$k.sorted.stats
grep -wFf $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.core $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.nodes.fasta.len  | awk '{sum += $2} END {print sum}' > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.categories.len
grep -wFf $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.private $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.nodes.fasta.len  | awk '{sum += $2} END {print sum}' >> $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.categories.len
grep -wFf $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.dispensable $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.nodes.fasta.len  | awk '{sum += $2} END {print sum}' >> $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.allSpecies.nodes.categories.len