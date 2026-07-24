#!/bin/bash
 
#SBATCH --job-name=jasmine
#SBATCH --ntasks=20
#SBATCH --array=0-18
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=80G
#SBATCH --error=/home/zajac/jasmine.%J.err
#SBATCH --output=/home/zajac/jasmine.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmem2new


source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate jasmine
DIR=/home/zajac
k=$SLURM_ARRAY_TASK_ID
s="s5000"
outdir=/home/zajac/SV_intersects/community${k}.${s}
SAMPLES=$(grep -v "##" $DIR/pggb_allcommunities/community${k}.${s}/community${k}.pggb.waved.vcf | head -n1 | sed 's/\t/\n/g')

if [ ! -d "$outdir" ]; then
  mkdir ${outdir}
fi


for SAMPLE in $SAMPLES;
do
  ###Merge any SVs within each sample if differ by 10bp
  echo /home/zajac/minigraph/community${k}/community${k}.${SAMPLE}.graph.waved.norm.SV.vcf > ${outdir}/list_of_files
  jasmine file_list=${outdir}/list_of_files --allow_intrasample max_dist=10 --output_genotypes out_file=${outdir}/temp.mng.vcf
  awk 'BEGIN{OFS="\t"} /VARCALLS=2/ {gsub(/\t0\|1:/,"\t1|1:") gsub(/\t1\|0:/,"\t1|1:") } {print}' ${outdir}/temp.mng.vcf > ${outdir}/temp2.mng.vcf
  echo /home/zajac/pggb_allcommunities/community${k}.${s}/community${k}.${SAMPLE}.pggb.waved.norm.SV.vcf > ${outdir}/list_of_files
  jasmine file_list=${outdir}/list_of_files --allow_intrasample max_dist=10 --output_genotypes out_file=${outdir}/temp.pggb.vcf
  awk 'BEGIN{OFS="\t"} /VARCALLS=2/ {gsub(/\t0\|1:/,"\t1|1:") gsub(/\t1\|0:/,"\t1|1:") } {print}' ${outdir}/temp.pggb.vcf > ${outdir}/temp2.pggb.vcf
  echo /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.${SAMPLE}.svimasm.norm.SV.vcf  > ${outdir}/list_of_files
  jasmine file_list=${outdir}/list_of_files --allow_intrasample max_dist=10 --output_genotypes out_file=${outdir}/temp.svim.vcf
  awk 'BEGIN{OFS="\t"} /VARCALLS=2/ {gsub(/\t0\|1:/,"\t1|1:") gsub(/\t1\|0:/,"\t1|1:") } {print}' ${outdir}/temp.svim.vcf > ${outdir}/temp2.svim.vcf

  echo ${outdir}/temp2.mng.vcf  > ${outdir}/list_of_files
  echo ${outdir}/temp2.pggb.vcf >> ${outdir}/list_of_files
  echo ${outdir}/temp2.svim.vcf >> ${outdir}/list_of_files

  ###Merge between methods - also tested max dist 100 and 10k 
  jasmine file_list=${outdir}/list_of_files  --normalize_type --allow_intrasample max_dist=1000 max_dist_linear=0.5 min_support=2 out_file=/home/zajac/SV_intersects/community${k}.${s}.minigraph.${SAMPLE}.DEFAULT.comp.vcf
  
  rm -R ${outdir}/*;
done

##Merge across all samples
genome_file=$(ls /home/zajac/chromosome_sets/community${k}_unmasked/rPodCre2.1#1#*.fasta)

ls /home/zajac/SV_intersects/community${k}.${s}.minigraph.*.DEFAULT.comp.vcf > /home/zajac/SV_intersects/list_of_files_DEFAULT_${k}
jasmine file_list=/home/zajac/SV_intersects/list_of_files_DEFAULT_${k} --ignore_merged_inputs min_seq_id=0.5 --output_genotypes --default_zero_genotype  --normalize_type out_file=/home/zajac/SV_intersects/largeSVs/merged.DEFAULT.${k}.vcf
