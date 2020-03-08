### compute matrix
file='balanced_peaklist_for_compare.stable.highsig.uniq'
tail -n+2 all_pk/LMqPCRnorm_decay_folder/$file'.txt' | cut -f1,2,3 > all_pk/LMqPCRnorm_decay_folder/$file'.bed'
time computeMatrix reference-point -S logIScorrespond240min_240min_merge_10000_qnorm_genomewide_12_1LDv3.bw -R all_pk/LMqPCRnorm_decay_folder/$file'.bed' -a 100000 -b 100000 -o all_pk/LMqPCRnorm_decay_folder/$file'.bed.txt.gz' --binSize 4000
gunzip all_pk/LMqPCRnorm_decay_folder/$file'.bed.txt.gz'

file='balanced_peaklist_for_compare.unstable.highsig.subs'
tail -n+2 all_pk/LMqPCRnorm_decay_folder/$file'.txt' | cut -f1,2,3 > all_pk/LMqPCRnorm_decay_folder/$file'.bed'
time computeMatrix reference-point -S logIScorrespond240min_240min_merge_10000_qnorm_genomewide_12_1LDv3.bw -R all_pk/LMqPCRnorm_decay_folder/$file'.bed' -a 100000 -b 100000 -o all_pk/LMqPCRnorm_decay_folder/$file'.bed.txt.gz' --binSize 4000
gunzip all_pk/LMqPCRnorm_decay_folder/$file'.bed.txt.gz'

file='balanced_peaklist_for_compare.unstable.lowsig.subs'
tail -n+2 all_pk/LMqPCRnorm_decay_folder/$file'.txt' | cut -f1,2,3 > all_pk/LMqPCRnorm_decay_folder/$file'.bed'
time computeMatrix reference-point -S logIScorrespond240min_240min_merge_10000_qnorm_genomewide_12_1LDv3.bw -R all_pk/LMqPCRnorm_decay_folder/$file'.bed' -a 100000 -b 100000 -o all_pk/LMqPCRnorm_decay_folder/$file'.bed.txt.gz' --binSize 4000
gunzip all_pk/LMqPCRnorm_decay_folder/$file'.bed.txt.gz'

time Rscript plotIS.R

#for i in {1..9}
#do
#echo $i
### get bed
#tail -n+2 all_pk/LMqPCRnorm__folder/'raw_kmeans.'$i'.pk.clean.txt' | cut -f1,2,3 > all_pk/LMqPCRnorm_decay_folder/'raw_kmeans.'$i'.pk.clean.bed'
### compute matrix
#time computeMatrix reference-point -S logIScorrespond240min_240min_merge_10000_qnorm_genomewide_12_1LDv3.bw -R all_pk/LMqPCRnorm_decay_folder/'raw_kmeans.'$i'.pk.clean.bed' -a 100000 -b 100000 -o all_pk/LMqPCRnorm_decay_folder/'raw_kmeans.'$i'.pk.clean.bed.IS.txt.gz' --binSize 4000
#gunzip all_pk/LMqPCRnorm_decay_folder/'raw_kmeans.'$i'.pk.clean.bed.IS.txt.gz'
#done



