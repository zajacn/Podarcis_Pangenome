#!/bin/bash
 
#SBATCH --job-name=protexclud
#SBATCH --ntasks=10
#SBATCH --array=1-15
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
#SBATCH --error=/home/zajac/protexclud.%J.err
#SBATCH --output=/home/zajac/protexclud.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=standard

module load perl/5.30.1

DIR=/home/zajac/Species_Specific_Masking
k=$SLURM_ARRAY_TASK_ID
name=`sed -n ${k}p < /home/zajac/busco_list.txt`

/data/biosoftware/ncbi-blast/ncbi-blast-2.16.0+/bin/blastx -num_threads 10 -query $DIR/${name}.db-families.fa -db /home/zajac/scripts/Software/uniprot_sprot.fasta  -max_target_seqs 5  -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen" -out $DIR/${name}.all_reapeats.blastx.out >  $DIR/${name}.all_reapeats.blastx.log 2> $DIR/${name}.all_repeats.blastx.err

# Eliminate all repeats that cover at least 50% of a protein from uniprot and have at least an 75% identity match
awk '$3 > 75 && $4 >= $13/2 {print $1}' $DIR/${name}.all_reapeats.blastx.out | sort | uniq > $DIR/${name}.temp
grep '>' $DIR/${name}.db-families.fa | awk '{print $1}' | sed 's/>//g' > $DIR/${name}.all
cat $DIR/${name}.all | grep -vf $DIR/${name}.temp > $DIR/${name}.temp2
perl /home/zajac/scripts/generally_useful/extractseqids.pl $DIR/${name}.temp2 $DIR/${name}.db-families.fa > $DIR/${name}.db-families.fanoprot
rm $DIR/${name}.temp $DIR/${name}.temp2 $DIR/${name}.all
