#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

source ~/.bash_profile

cd /storage/home/gzx103/group/projects/ctcf_auxin/
cd /storage/home/gzx103/scratch/ctcf_auxin

############################################################
#######(((1)))###### get bed file of peak lists
### (1) get merged_peak list
cat peaks/*_CTCF*.bed | sort -u | sort -k1,1 -k2,2n > all.bed
### get mid point
cat all.bed | awk -F '\t' -v OFS='\t' '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1}' | sort -k1,1 -k2,2n > all.mid.bed
### merge mid point if distance <= 250-bp
bedtools merge -i all.mid.bed -d 250 > all.merged.mid.bed
### expand mid point to peak (+- 250-bp)
cat all.merged.mid.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-250, int(($2+$3)/2)+250, "ctcf_"$1"_"$2"_"$3}' > allmerged.bed
sort -k1,1 -k2,2n allmerged.bed > allmerged.sorted.bed
rm allmerged.bed
rm all.merged.mid.bed
rm all.mid.bed
rm all.bed

### (2) get qPCR peak
cat 20180910_CTCF_ChIPqPCR_7pts.txt | awk -F '\t' -v OFS='\t' '{if($6!="NA") print $0}' > 20180910_CTCF_ChIPqPCR_7pts.rmNA.txt
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.rmNA.txt | awk -F '\t' -v OFS='\t' '{if($6!="NA" && $5!="Weakadasd") print $2, int(($3+$4)/2)-250, int(($3+$4)/2)+250,"qPCR_"$2"_"$3"_"$4}' | sort -k1,1 -k2,2n > qtPCR_regions.sorted.6TP.0729.withID.bed
### get qPCR signal mat
tail -n+2 20180910_CTCF_ChIPqPCR_7pts.rmNA.txt | awk -F '\t' -v OFS='\t' '{if($6!="NA" && $5!="Weakadasd") print $2, int(($3+$4)/2)-250, int(($3+$4)/2)+250,"qPCR_"$2"_"$3"_"$4, $8+$9, $10+$11, $12+$13, $14+$15, $16+$17, $18+$19}' | sort -k1,1 -k2,2n > 20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.sig.txt


### (3) get non_peak
bedtools makewindows -g ~/group/genome/mm9/mm9.chrom.sizes -w 500 > mm9.500bp.randomwin.bed
bedtools intersect -a mm9.500bp.randomwin.bed -b allmerged.sorted.bed -v > mm9.500bp.randomwin.bg.bed
shuf -n 120000 mm9.500bp.randomwin.bg.bed > mm9.500bp.randomwin.120000.bg.bed
cat mm9.500bp.randomwin.120000.bg.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,"randbg_"$1"_"$2"_"$3}' | head -100000 | sort -k1,1 -k2,2n > mm9.500bp.randomwin.100000.bg.sort.bed
rm mm9.500bp.randomwin.bed
rm mm9.500bp.randomwin.bg.bed
rm mm9.500bp.randomwin.120000.bg.bed

### get all peak list pooled
cat allmerged.sorted.bed qtPCR_regions.sorted.6TP.0729.withID.bed mm9.500bp.randomwin.100000.bg.sort.bed > ctcf.qPCR.randbg.bed
bedtools intersect -a ctcf.qPCR.randbg.bed -b blacklist/mm9.Input.blacklist.sort.merged.bed -v > ctcf.qPCR.randbg.blackrm.bed
rm ctcf.qPCR.randbg.bed

sort -k1,1 ctcf.qPCR.randbg.blackrm.bed > ctcf.qPCR.randbg.blackrm.idsort.bed
cat ctcf.qPCR.randbg.blackrm.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500, $1"_"$2"_"$3}' | sort -k1,1 > ctcf.qPCR.randbg.blackrm.idsort.1kb.bed
cat ctcf.qPCR.randbg.blackrm.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' | sort -k1,1 > ctcf.qPCR.randbg.blackrm.idsort.5kb.bed
cat ctcf.qPCR.randbg.blackrm.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' | sort -k1,1 > ctcf.qPCR.randbg.blackrm.idsort.10kb.bed

rm ctcf.qPCR.randbg.blackrm.bed

mkdir ctcf_master_pk
mv ctcf.qPCR.randbg.blackrm.idsort.bed ctcf_master_pk/
mv ctcf.qPCR.randbg.blackrm.idsort.1kb.bed ctcf_master_pk/
mv ctcf.qPCR.randbg.blackrm.idsort.5kb.bed ctcf_master_pk/
mv ctcf.qPCR.randbg.blackrm.idsort.10kb.bed ctcf_master_pk/

mv mm9.500bp.randomwin.100000.bg.sort.bed ctcf_master_pk/
mv qtPCR_regions.sorted.6TP.0729.withID.bed ctcf_master_pk/
mv allmerged.sorted.bed ctcf_master_pk/


