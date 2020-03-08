bedtools makewindows -g ~/group/genome/mm9/mm9.chrom.sizes -w 500 > mm9.500bp.randomwin.bed
cat mm9.500bp.randomwin.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' > mm9.500bp.randomwin.withid.bed

### get IP bedgraph
while read LINE
do
	set -- $LINE
	echo $1
	~/group/software/ucsc/bigWigAverageOverBed bam_bw_bedgraph/$1'IP.rawRC.1bp.bw' mm9.500bp.randomwin.withid.bed tmp.tab
	cut -f1,5 tmp.tab | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' > bam_bw_bedgraph/$1'IP.rawRC.1bp.bedgraph'
done < bam_bw_bedgraph/file_list.all.txt

### get control bedgraph 
while read LINE
do
	set -- $LINE
	echo $2
	~/group/software/ucsc/bigWigAverageOverBed bam_bw_bedgraph/$2'CTRL.rawRC.1bp.bw' mm9.500bp.randomwin.withid.bed tmp.tab
	cut -f1,5 tmp.tab | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' > bam_bw_bedgraph/$2'CTRL.rawRC.1bp.bedgraph'
done < bam_bw_bedgraph/file_list.all.txt

### get bgsub bedgraph
cut -f1,2,3 bam_bw_bedgraph/1579_0A_CTCF_mm9_duprm.bamIP.rawRC.1bp.bedgraph > bam_bw_bedgraph/IP.rawRC.bgsub.mat.txt
cut -f1,2,3 bam_bw_bedgraph/1578_0A_Input_mm9_duprm.bamCTRL.rawRC.1bp.bedgraph > bam_bw_bedgraph/CTRL.rawRC.bgsub.mat.txt
while read LINE
do
	set -- $LINE
	echo $1
	cut -f4 bam_bw_bedgraph/$1'IP.rawRC.1bp.bedgraph' > tmp1.txt
	paste bam_bw_bedgraph/IP.rawRC.bgsub.mat.txt tmp1.txt > bam_bw_bedgraph/IP.rawRC.bgsub.mat.txt.tmp && mv bam_bw_bedgraph/IP.rawRC.bgsub.mat.txt.tmp bam_bw_bedgraph/IP.rawRC.bgsub.mat.txt
	cut -f4 bam_bw_bedgraph/$2'CTRL.rawRC.1bp.bedgraph' > tmp1.txt
	paste bam_bw_bedgraph/CTRL.rawRC.bgsub.mat.txt tmp1.txt > bam_bw_bedgraph/CTRL.rawRC.bgsub.mat.txt.tmp && mv bam_bw_bedgraph/CTRL.rawRC.bgsub.mat.txt.tmp bam_bw_bedgraph/CTRL.rawRC.bgsub.mat.txt
done < bam_bw_bedgraph/file_list.all.txt

### get signal bgsubtracted
time Rscript get_bg_subtracted_sigmat_bedgraph.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.0kb.sigmat.txt bam_bw_bedgraph/IP.rawRC.bgsub.mat.txt bam_bw_bedgraph/CTRL.rawRC.bgsub.mat.txt bam_bw_bedgraph/IP.rawRC.bgsub0.mat.txt

### get qPCR_LMnorm bedgraph
time Rscript qPCR_LMnorm_bw.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt bam_bw_bedgraph/file_list.all.timepoints.txt 20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.sig.txt 20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.timepoints.txt bam_bw_bedgraph/IP.rawRC.bgsub0.mat.txt bam_bw_bedgraph/file_list.all.txt

### get each bw file
while read LINE
do
	set -- $LINE
	echo $1
	time sort -k1,1 -k2,2n $1'.bgsub.qpcrnorm.bedgraph' > $1'.bgsub.qpcrnorm.sort.bedgraph'
	time ~/group/software/ucsc/bedGraphToBigWig $1'.bgsub.qpcrnorm.sort.bedgraph' ~/group/genome/mm9/mm9.chrom.sizes $1'.bgsub.qpcrnorm.sort.bw'
	#rm $1'.bgsub.qpcrnorm.bedgraph' 
	#rm $1'.bgsub.qpcrnorm.sort.bedgraph' 
	#mv $1'.bgsub.qpcrnorm.sort.bw' bam_bw_bedgraph/
done < bam_bw_bedgraph/file_list.all.txt








