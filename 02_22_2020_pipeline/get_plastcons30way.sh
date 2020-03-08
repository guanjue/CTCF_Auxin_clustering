#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/gzx103/group/genome/mm9/plastcons

### download files from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=715959001_rhLWLGxtdFx7ANtoa4a4r8RrqCbn&c=chr12&g=cons30way
#rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/mm9/phastCons30way/euarchontoglires/ ./

### unzip files
#gunzip *.gz

### get file list
#ls *.pp.data > plastcons30way_list.txt

### get bigwig files
for file in $(cat plastcons30way_list.txt)
do
	echo $file
	#time ~/group/software/ucsc/wigToBigWig $file ../mm9.1to19_X.genome $file'.bw'
done

### merge bigwig
time ~/group/software/ucsc/bigWigMerge *bw plastcons30way.euarchontoglires.bedgraph

time sort -k1,1 -k2,2n plastcons30way.euarchontoglires.bedgraph > plastcons30way.euarchontoglires.sort.bedgraph
#time ~/group/software/ucsc/bedGraphToBigWig plastcons30way.euarchontoglires.sort.bedgraph ../mm9.1to19_X.genome phyloP.bw
time ~/group/software/ucsc/bedGraphToBigWig plastcons30way.euarchontoglires.sort.bedgraph ../mm9.1to19_X.genome plastcons30way.euarchontoglires.bw



