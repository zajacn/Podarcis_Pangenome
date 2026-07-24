for name in rPodSic rPodBoc rPodVau rPodLio rPodMel rPodErh rPodGai rPodPit rPodTil rPodFil rPodMur119; 
do 
  for x in {0..18}; 
  do 
    singularity exec --cleanenv /data/biosoftware/cactus/cactus.sif halSynteny --queryGenome "${name}1_2" --targetGenome "${name}1_1" /home/zajac/cactus/community${x}/community${x}.hal community${x}.${name}.psl; 
  done;
done

for name in rPodSic rPodBoc rPodVau rPodLio rPodMel rPodErh rPodGai rPodPit rPodTil rPodFil rPodMur119; 
do  
  cat *.${name}.psl > ${name}_all.psl;
done

