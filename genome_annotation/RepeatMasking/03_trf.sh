#!/bin/bash
 
#SBATCH --job-name=trf
#SBATCH --ntasks=15
#SBATCH --array=1-15
#SBATCH --nodes=3
#SBATCH --time=6-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/trf.%J.err
#SBATCH --output=/home/zajac/trf.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

DIR=/home/zajac/Species_Specific_Masking
module load python/3.9.13
module load perl/5.30.1
export PATH=$PATH:/data/biosoftware/ncbi-blast/ncbi-blast-2.16.0+/bin/

k=$SLURM_ARRAY_TASK_ID
name=`sed -n ${k}p < /home/zajac/busco_list.2.txt`


#Run TRF
name_new=$(basename $name .fasta)
ln -s $DIR/${name}.masked $DIR/TRF/${name_new}
perl /home/zajac/scripts/Software/Augustus_Gaia_Sup_Scripts/splitMfasta.pl --minsize=25000000 $DIR/TRF/${name_new}
ls $DIR/TRF/${name_new}.split.*.fa | parallel '/data/biosoftware/trf/bin/trf {} 2 7 7 80 10 50 500 -d -m -h &> {}.log'
ls $DIR/TRF/${name_new}.split.*.fa.2.7.7.80.10.50.500.dat | parallel 'python /home/zajac/scripts/Software/Augustus_Gaia_Sup_Scripts/parseTrfOutput.py {} --minCopies 1 \
    --statistics {}.STATS > {}.raw.gff 2> {}.parsedLog'
ls $DIR/TRF/${name_new}.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff | parallel 'sort -k1,1 -k4,4n -k5,5n {} \
   > {}.sorted 2> {}.sortLog'
export FILES=$DIR/TRF/${name_new}.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff.sorted
for f in $FILES;
do
  /data/biosoftware/bedtools2/bedtools-2.30.0/bin/bedtools merge -i $f | awk 'BEGIN{OFS="\t"} {print $1,"trf","repeat",$2+1,$3,".",".",".","."}' \
  > $f.merged.gff 2> $f.bedtools_merge.log;
done
ls $DIR/TRF/${name_new}.split.*.fa | parallel '/data/biosoftware/bedtools2/bedtools-2.30.0/bin/bedtools maskfasta -fi {} -bed \
  {}.2.7.7.80.10.50.500.dat.raw.gff.sorted.merged.gff -fo {}.combined.masked -soft &> {}.bedools_mask.log'
cat $DIR/TRF/${name_new}.split.*.fa.combined.masked > $DIR/TRF/${name_new}.fa.combined.masked
