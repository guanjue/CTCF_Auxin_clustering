### pks signal mat
cp ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.bed ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt

while read LINE
do
	set -- $LINE
	echo $1
	~/group/software/ucsc/bigWigAverageOverBed bam_bw_bedgraph/$1'IP.rawRC.1bp.bw' ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.bed tmp.tab
	sort -k1,1 tmp.tab | cut -f5 > tmp.tab.txt
	paste ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt tmp.tab.txt > ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt.tmp
	mv ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt.tmp ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt
done < bam_bw_bedgraph/file_list.all.txt

### 1kb signal mat
cp ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.1kb.bed ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt

while read LINE
do
	set -- $LINE
	echo $2
	~/group/software/ucsc/bigWigAverageOverBed bam_bw_bedgraph/$2'CTRL.rawRC.1bp.bw' ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.1kb.bed tmp.tab
	sort -k1,1 tmp.tab | cut -f5 > tmp.tab.txt
	paste ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt tmp.tab.txt > ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt.tmp
	mv ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt.tmp ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt
done < bam_bw_bedgraph/file_list.all.txt

### 5kb signal mat
cp ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.5kb.bed ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt

while read LINE
do
	set -- $LINE
	echo $2
	~/group/software/ucsc/bigWigAverageOverBed bam_bw_bedgraph/$2'CTRL.rawRC.1bp.bw' ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.5kb.bed tmp.tab
	sort -k1,1 tmp.tab | cut -f5 > tmp.tab.txt
	paste ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt tmp.tab.txt > ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt.tmp
	mv ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt.tmp ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt
done < bam_bw_bedgraph/file_list.all.txt

### 10kb signal mat
cp ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.10kb.bed ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt

while read LINE
do
	set -- $LINE
	echo $2
	~/group/software/ucsc/bigWigAverageOverBed bam_bw_bedgraph/$2'CTRL.rawRC.1bp.bw' ctcf_master_pk/ctcf.qPCR.randbg.blackrm.idsort.10kb.bed tmp.tab
	sort -k1,1 tmp.tab | cut -f5 > tmp.tab.txt
	paste ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt tmp.tab.txt > ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt.tmp
	mv ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt.tmp ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt
done < bam_bw_bedgraph/file_list.all.txt



