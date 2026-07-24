#!/bin/bash
 
#SBATCH --job-name=blast
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --array=1-11
#SBATCH --error=/home/zajac/blast.err
#SBATCH --output=/home/zajac/blast.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new


k=$SLURM_ARRAY_TASK_ID 
ref=`sed -n ${k}p < /home/zajac/CompareHaplotypes/BLAST_GenicHemizygosity/species.txt`
DIR=/groups/mpistaff/Zajac/proteomes/data/ALL_Tiberius_annotated
OUTDIR=/home/zajac/CompareHaplotypes/BLAST_GenicHemizygosity
blast=/data/biosoftware/ncbi-blast/ncbi-blast-2.15.0+/bin

name1=$(echo ${ref}"#1#Z")
name2=$(echo ${ref}"#1#W")
awk -v name1=$name1 -v name2=$name2 '$1 != name1 && $1 != name2 {print}' $DIR/${ref}#1.gtf | awk '$3 == "gene" {print}' | awk '{print $9".t1"}' > $OUTDIR/${ref}.hap1_genes
perl /home/zajac/scripts/generally_useful/extractseqids.pl $OUTDIR/${ref}.hap1_genes $DIR/${ref}#1.all.aa > $OUTDIR/${ref}#1.all.subset.aa

$blast/makeblastdb -in $OUTDIR/${ref}#1.all.subset.aa -dbtype prot -title ${ref}#1
$blast/blastp -num_threads 8 -query $DIR/${ref}#2.all.aa  -db $OUTDIR/${ref}#1.all.subset.aa -max_target_seqs 5 -evalue 1e-10 -outfmt 6 -out $OUTDIR/${ref}.blastn.tsv

ln -s $DIR/${ref}#2.all.aa $OUTDIR/
$blast/makeblastdb -in $OUTDIR/${ref}#2.all.aa -dbtype prot -title ${ref}#2
$blast/blastp -num_threads 8 -query  $OUTDIR/${ref}#1.all.subset.aa  -db $OUTDIR/${ref}#2.all.aa -max_target_seqs 5 -evalue 1e-10 -outfmt 6 -out $OUTDIR/${ref}.blastn.RB.tsv
