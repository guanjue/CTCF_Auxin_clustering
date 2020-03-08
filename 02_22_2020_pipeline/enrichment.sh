### 
cd ~/group/projects/vision/figure_7_state_transition_ct_0807
~/group/software/ucsc/liftOver lsk13_eryad13.bed ~/group/genome/mm10/mm10ToMm9.over.chain.gz lsk13_eryad13.mm9.bed lsk13_eryad13.notmm9.bed
~/group/software/ucsc/liftOver lsk13_eryad7.bed ~/group/genome/mm10/mm10ToMm9.over.chain.gz lsk13_eryad7.mm9.bed lsk13_eryad7.notmm9.bed
cd ~/group/projects/vision/state13_to_7_revision_HPC7_vs_ER4
~/group/software/ucsc/liftOver ER4.OnTad.joint.boundary.bed ~/group/genome/mm10/mm10ToMm9.over.chain.gz ER4.OnTad.joint.boundary.mm9.bed ER4.OnTad.joint.boundary.notmm9.bed

### sort independent peak
sort -k1,1 -k2,2n lsk13_eryad13.mm9.bed > lsk13_eryad13.mm9.sort.bed
sort -k1,1 -k2,2n lsk13_eryad7.mm9.bed > lsk13_eryad7.mm9.sort.bed
sort -k1,1 -k2,2n ER4.OnTad.joint.boundary.mm9.bed > ER4.OnTad.joint.boundary.mm9.sort.bed


cd /storage/home/gzx103/scratch/ctcf_auxin
### sort peak
for i in {1..3}
do
echo $i
sort -k1,1 -k2,2n ~/scratch/ctcf_auxin/all_pk/LMqPCRnorm_decay_folder/'balanced_peaklist_for_compare.'$i'.pk.clean.bed' > ~/scratch/ctcf_auxin/all_pk/LMqPCRnorm_decay_folder/'balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed'
done

for i in {1..3}
do
echo $i
bedtools intersect -a independent_pk/lsk13_eryad13.mm9.sort.bed -b all_pk/LMqPCRnorm_decay_folder/'balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed' -wa -u > all_pk/LMqPCRnorm_decay_folder/'lsk13_eryad13.balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed'
bedtools intersect -a independent_pk/lsk13_eryad7.mm9.sort.bed -b all_pk/LMqPCRnorm_decay_folder/'balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed' -wa -u > all_pk/LMqPCRnorm_decay_folder/'lsk13_eryad7.balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed'
done

for i in {1..5}
do
	echo $i
	cat independent_pk/ER4.OnTad.joint.boundary.mm9.sort.bed | awk -F '\t' -v OFS='\t' -v tl=$i '{if ($4==tl) print $0}' > independent_pk/'ER4.OnTad.joint.boundary.mm9.sort.l'$i'.bed'
done


for i in {1..5}
do
echo $i
for j in {1..3}
do
bedtools intersect -a independent_pk/'ER4.OnTad.joint.boundary.mm9.sort.l'$i'.bed' -b all_pk/LMqPCRnorm_decay_folder/'balanced_peaklist_for_compare.'$j'.pk.clean.sort.bed' -wa -u > independent_pk/'ER4.OnTad.joint.boundary.mm9.sort.l'$i'.c'$j'.bed'
done
done

for i in {1..5}
do
echo $i
wc -l independent_pk/'ER4.OnTad.joint.boundary.mm9.sort.l'$i'.c'*'.bed'
done

for i in {1..3}
do
sort -u all_pk/LMqPCRnorm_decay_folder/'balanced_peaklist_for_compare.'$i'.pk.clean.bed' | wc -l
done


for i in {1..3}
do
wc -l all_pk/LMqPCRnorm_decay_folder/'lsk13_eryad13.balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed'
done

for i in {1..3}
do
wc -l all_pk/LMqPCRnorm_decay_folder/'lsk13_eryad7.balanced_peaklist_for_compare.'$i'.pk.clean.sort.bed'
done

