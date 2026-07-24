##Community detection
source /data/modules/python/python-miniconda3-2025-02/etc/profile.d/conda.sh
conda activate wfmash

###Mapping
wfmash all_Podarcis_genomes.fasta.gz -p 90 -n 15 -t 4 -m > all_Podarcis_genomes.mapping.paf

###Mapping to network
python3 /home/zajac/scripts/Software/pggb/pggb/scripts/paf2net.py -p all_Podarcis_genomes.mapping.paf

###Network to communities
python3 /home/zajac/scripts/Software/pggb/pggb/scripts/net2communities.py -e all_Podarcis_genomes.mapping.paf.edges.list.txt -w all_Podarcis_genomes.mapping.paf.edges.weights.txt -n all_Podarcis_genomes.mapping.paf.vertices.id2name.txt

###Summarise
seq 0 19 | while read i; do
    chromosomes=$(cat all_Podarcis_genomes.mapping.paf.edges.weights.txt.community.$i.txt | cut -f 3 -d '#' | sort | uniq | tr '\n' ' ');
    echo "community $i --> $chromosomes";
done

###Plot
python3 /home/zajac/scripts/Software/pggb/pggb/scripts/net2communities.py -e scerevisiae7.mapping.paf.edges.list.txt -w scerevisiae7.mapping.paf.edges.weights.txt -n scerevisiae7.mapping.paf.vertices.id2name.txt --plot