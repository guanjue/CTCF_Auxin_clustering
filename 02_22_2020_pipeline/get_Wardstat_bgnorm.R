library(pheatmap)
library(mclust)
library(DESeq2)
library(MASS)
library(data.table)

### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_mat_file = args[1]
et_mat_output_file = args[2]

#signal_mat_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.txt'
#et_mat_output_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.Ward.txt'


### get input mat
chipseq_mat = as.data.frame(fread(signal_mat_file))


### get R2
get_R2 = function(x1, x2){
	R2 = 1 - mean((x1-x2)^2)/mean((x1-mean(x1))^2)
	return(R2)
}
###### Use R2 to remove replciates
### R2 0A
get_R2(chipseq_mat[,5],chipseq_mat[,6])
get_R2(chipseq_mat[,5],chipseq_mat[,7])
get_R2(chipseq_mat[,6],chipseq_mat[,7])
### R2 4A
get_R2(chipseq_mat[,8],chipseq_mat[,9])
### R2 6A
get_R2(chipseq_mat[,10],chipseq_mat[,11])
### R2 12A
get_R2(chipseq_mat[,12],chipseq_mat[,13])
get_R2(chipseq_mat[,12],chipseq_mat[,14])
get_R2(chipseq_mat[,13],chipseq_mat[,14])
### R2 18A
get_R2(chipseq_mat[,15],chipseq_mat[,16])
get_R2(chipseq_mat[,15],chipseq_mat[,17])
get_R2(chipseq_mat[,16],chipseq_mat[,17])
### R2 24A
get_R2(chipseq_mat[,18],chipseq_mat[,19])
get_R2(chipseq_mat[,18],chipseq_mat[,20])
get_R2(chipseq_mat[,19],chipseq_mat[,20])

### select good replicates
used_replicates = c(5,6, 8,9, 10,11, 12,13, 15,17, 18,20)
chipseq_mat_colused = cbind(chipseq_mat[,1:4], chipseq_mat[,used_replicates])
chipseq_mat_used = chipseq_mat_colused[(grepl("ctcf", chipseq_mat_colused[,4])),]

### DESeq2
d1_0 = chipseq_mat_used[,-c(1:4)]
d1_TP_condition = cbind(c('0A','0A','4A','4A','6A','6A','12A','12A','18A','18A','24A','24A'))
rownames(d1_TP_condition) = colnames(d1_0)
colnames(d1_TP_condition) = c('condition')

### get DEseq object
dds <- DESeqDataSetFromMatrix(countData = round(d1_0,0),
                              colData = d1_TP_condition,
                              design= ~ condition)

### set reference
dds$condition <- relevel(dds$condition, ref = "0A")

### run DEseq2
dds <- DESeq(dds)

### set sizeFactors to 1
sizeFactors(dds) = rep(1.0, 12)

### run DEseq2 without normalization
dds <- DESeq(dds, test='Wald')

### get Ward stat
res_0_4 <- results(dds, contrast=c("condition","4A","0A"))
res_0_6 <- results(dds, contrast=c("condition","6A","0A"))
res_0_12 <- results(dds, contrast=c("condition","12A","0A"))
res_0_18 <- results(dds, contrast=c("condition","18A","0A"))
res_0_24 <- results(dds, contrast=c("condition","24A","0A"))

### matrix for clustering
et_mat = cbind(
        res_0_4[,4],
        res_0_6[,4],
        res_0_12[,4],
        res_0_18[,4],
        res_0_24[,4]
)

### output Ward statistics
et_mat_output = cbind(chipseq_mat_used[,c(1:4)], et_mat)
colnames(et_mat_output) = c('chr', 'start', 'end', 'id', '4vs0', '6vs0', '12vs0', '18vs0', '24vs0')
write.table(et_mat_output, et_mat_output_file, quote=F, col.names=T, row.names=F, sep='\t')


