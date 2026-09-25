#!/bin/bash
#SBATCH --job-name=Dsuite
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --array=0-8,10-18
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --error=/home/zajac/Dsuite.%J.err
#SBATCH --output=/home/zajac/Dsuite.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zajac@evolbio.mpg.de
#SBATCH --partition=highmemnew

k=$SLURM_ARRAY_TASK_ID 

## run Dsuite Dtrios and Combine 
echo $k
/home/zajac/scripts/Software/Dsuite/Build/Dsuite Dtrios /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.combined.forDsuite.vcf  /home/zajac/SYRI/Dsuite/SETS.Dsuite.txt -t /home/zajac/SYRI/Dsuite/tree_for_Dsuite.nwk  --ABBAclustering --out-prefix community${k}

/home/zajac/scripts/Software/Dsuite/Build/Dsuite DtriosCombine community0 community1 community2 community3 community4 community5 community6 community7 community8 community10 community11 community12 community13 community14 community15 community16 community17 community18 --tree=/home/zajac/SYRI/Dsuite/tree_for_Dsuite.nwk --out-prefix=communities.combined

## run Fbrach based on Patterson's D and visualise
/home/zajac/scripts/Software/Dsuite/Build/Dsuite Fbranch /home/zajac/SYRI/Dsuite/tree_for_Dsuite.nwk /home/zajac/SYRI/Dsuite/communities.combined_combined_tree.txt > /home/zajac/SYRI/Dsuite/Fbranch_output.txt

module load python/3.9.13
python /home/zajac/scripts/Software/Dsuite/utils/dtools.py /home/zajac/SYRI/Dsuite/Fbranch_output.txt /home/zajac/SYRI/Dsuite/tree_for_Dsuite.nwk --outgroup LacAg

## run Dinvestigate 
/home/zajac/scripts/Software/Dsuite/Build/Dsuite Dinvestigate --run-name community${k}.Dinvestigate --window=100,100 /home/zajac/SYRI/community${k}/community${k}.rPodCre2.1_to_all/community${k}.combined.forDsuite.vcf /home/zajac/SYRI/Dsuite/SETS.Dsuite.txt /home/zajac/SYRI/Dsuite/test_trios.txt
