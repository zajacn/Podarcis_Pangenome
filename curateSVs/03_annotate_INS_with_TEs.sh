#!/bin/bash
 
#SBATCH --job-name=blast_INS
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --array=1-332%50
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
#SBATCH --error=/home/zajac/blast_INS.%J.err
#SBATCH --output=/home/zajac/blast_INS.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=standard

k=$SLURM_ARRAY_TASK_ID
i=`sed -n ${k}p < /home/zajac/Species_Specific_Masking/transposable_elements_fasta/species`
DIR=/home/zajac/Species_Specific_Masking/transposable_elements_fasta
conda activate seqkit

#Annotate INS with output from RepeatModeler/RepeatMasker
grep 'SVTYPE=INS' /home/zajac/SV_intersects/largeSVs/merged.DEFAULT.*.vcf | awk '{print $1"_"$2"\n"$5}' | sed 's"/home/zajac/SV_intersects/largeSVs/">"g' > $DIR/INS.fasta
awk '{print $1"#"$2"#"$3"\t"$4"\t"$5}' $DIR/../repeats.allhaplotypes.tsv > $DIR/../transposable_elements.bed
grep $i $DIR/../transposable_elements.bed > $DIR/${i}.bed
bedtools getfasta -fi /groups/mpistaff/Zajac/genomes/ncbi_dataset/data/ALL_corrected_masked/${i}.all.fasta -bed $DIR/${i}.bed > $DIR/${i}.fa
seqkit rmdup -s < $DIR/${i}.fa > $DIR/${i}.dedup.fa
makeblastdb -in $DIR/${i}.dedup.fa -parse_seqids -blastdb_version 5 -dbtype nucl
blastn -query $DIR/insertions/${k}.fasta -db $DIR/${i}.dedup.fa -max_target_seqs 1 -evalue 1e-10 -outfmt 6 -out $DIR/insertions/INS.${i}.${k}.dedup.fa
