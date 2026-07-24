#!/bin/bash
 
#SBATCH --job-name=ipggb_workflow_graph
#SBATCH --ntasks=15
#SBATCH --array=3
#SBATCH --time=2-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/ipggb_workflow.%J.err
#SBATCH --output=/home/zajac/ipggb_workflow.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=global
 
DIR=/home/zajac
k=$SLURM_ARRAY_TASK_ID

source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate wfmash
impg=~/.cargo/bin/impg
export PATH=$PATH:/home/zajac/scripts/Software/pggb/pggb/
pggb=/home/zajac/scripts/Software/pggb/pggb_updated.sif
gfalace=/home/zajac/scripts/Software/gfalace/target/release/gfalace

##All_to_All Mapping
wfmash -t 3 $DIR/chromosome_sets/community${k}_unmasked/community${k}.combined.fasta -p 90 -k 19 -J 0.001 -s 5000 > $DIR/pggb_allcommunities/community${k}.s5000.paf

##Partition
$impg partition -p $DIR/pggb_allcommunities/community${k}.s5000.paf -w 1000000 -o fasta --sequence-files $DIR/chromosome_sets/community${k}_unmasked/community${k}.combined.fasta --output-folder $DIR/pggb_allcommunities/community${k}.s5000 --threads 15

##Local map and seqwish
nb_partitions=$(ls $DIR/pggb_allcommunities/community${k}.s5000/partition*.fasta | wc -l)
for ((m=1; m<=nb_paritions; m++));
do
  /data/biosoftware/samtools/samtools-1.20/samtools faidx $DIR/pggb_allcommunities/community${k}.s5000/partition${m}.fasta
  wfmash -t 15  $DIR/pggb_allcommunities/community${k}.s5000/partition${m}.fasta -p 90 > $DIR/pggb_allcommunities/community${k}.s5000/partition${m}.fasta.paf
  singularity run $pggb seqwish -p $DIR/pggb_allcommunities/community${k}.s5000/partition${m}.fasta.paf --seqs=$DIR/pggb_allcommunities/community${k}.s5000/partition${m}.fasta --gfa=$DIR/pggb_allcommunities/community${k}.s5000/partition${m}.fasta.gfa -k 0;
done

##Lace
$gfalace --gfa-files $DIR/pggb_allcommunities/community${k}.s1000/*.gfa --output $DIR/pggb_allcommunities/community${k}.s1000/community${k}.gfa --num-threads 64

##Sort
singularity run $pggb odgi sort -i $DIR/pggb_allcommunities/community${k}.s5000/community${k}.gfa -P -Y --threads=64 -o $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.og

##Visualise
singularity run $pggb odgi layout -i $DIR/pggb_allcommunities/community${k}.s10000/community${k}.sorted.og -o $DIR/pggb_allcommunities/community${k}.s10000/community${k}.sorted.og.lay -T $DIR/pggb_allcommunities/community${k}.s10000/community${k}.sorted.og.lay.tsv -t 10 -P

singularity run $pggb odgi draw -i $DIR/pggb_allcommunities/community${k}.s10000/community${k}.sorted.og -c $DIR/pggb_allcommunities/community${k}.s10000/community${k}.sorted.og.lay --png=$DIR/pggb_allcommunities/community${k}.s10000/community${k}.sorted.2D.png
 
singularity run $pggb odgi viz -i $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.og -o $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.png -z --threads=64

##Stats
singularity run $pggb odgi stats -i $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.og -S -a '*',0 --path-statistics -s --threads=15 > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.stats 

singularity run $pggb odgi similarity -i $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.og --distances -t 16 -D '#' -p 2 -P > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.smooth.final.dist.tsv

##Convert to GFA
singularity run $pggb odgi view -i $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.og --to-gfa --threads=20 > $DIR/pggb_allcommunities/community${k}.s5000/community${k}.sorted.gfa

