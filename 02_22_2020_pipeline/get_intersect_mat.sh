cd /Users/universe/Downloads/compare_two_methods

for i in {1..6}
do
echo $i
tail -n+2 'raw_mclust.'$i'.pk.clean.txt' | cut -f1,2,3 > 'raw_mclust.'$i'.pk.clean.bed'
done

tail -n+2 stable_peaks.txt.unstable.txt | cut -f1,2,3,4 > unstable.bed
tail -n+2 stable_peaks.txt | cut -f1,2,3,4 > stable_peaks.bed

cp unstable.bed unstable.mat.txt
cp stable_peaks.bed stable.mat.txt

for i in {1..6}
do
echo $i
### unstable.mat.txt
bedtools intersect -a unstable.bed -b 'raw_mclust.'$i'.pk.clean.bed' -c > tmp1.txt
cut -f5 tmp1.txt > tmp1a.txt
paste unstable.mat.txt tmp1a.txt > unstable.mat.txt.tmp && mv unstable.mat.txt.tmp unstable.mat.txt
### stable.mat.txt
bedtools intersect -a stable_peaks.bed -b 'raw_mclust.'$i'.pk.clean.bed' -c > tmp1.txt
cut -f5 tmp1.txt > tmp1a.txt
paste stable.mat.txt tmp1a.txt > stable.mat.txt.tmp && mv stable.mat.txt.tmp stable.mat.txt
done

time Rscript get_cluster_Decay_curve_intersect.R
