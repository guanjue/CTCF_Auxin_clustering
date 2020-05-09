library(pheatmap)
library(mclust)
library(data.table)

### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_mat_file = args[1]
signal_mat_tp_file = args[2]
qPCR_mat_file = args[3]
qPCR_mat_tp_file = args[4]
output_file_name = args[5]

signal_mat_file_raw = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt'
signal_mat_file_norm = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.txt'

used_replicates = c(5,6, 8,9, 10,11, 12,13, 15,17, 18,20)


### function
getR2 = function(x1,x2){
	r2 = 1-mean((x1-x2)^2)/mean((x1-mean(x1))^2)
	return(r2)
}

### get input mat
chipseq_sig_mat_raw = as.data.frame(fread(signal_mat_file_raw))
chipseq_sig_mat_norm = as.data.frame(fread(signal_mat_file_norm))
used_row = (grepl("ctcf", chipseq_sig_mat_raw[,4]))

chipseq_sig_mat_raw = chipseq_sig_mat_raw[,used_replicates]
chipseq_sig_mat_norm = chipseq_sig_mat_norm[,used_replicates]

r2_mat = c()
### raw R2
r2_vec = c()
for (i in seq(1,11,2)){
	r2_tmp = getR2(chipseq_sig_mat_raw[used_row,i],chipseq_sig_mat_raw[used_row,i+1])
	print(r2_tmp)
	r2_vec = c(r2_vec, r2_tmp)
}
r2_mat = cbind(r2_mat, r2_vec)

### qPCR R2
r2_vec = c()
for (i in seq(1,11,2)){
	r2_tmp = getR2(chipseq_sig_mat_norm[used_row,i],chipseq_sig_mat_norm[used_row,i+1])
	print(r2_tmp)
	r2_vec = c(r2_vec, r2_tmp)
}
r2_mat = cbind(r2_mat, r2_vec)

### BGnorm
r2_vec = c()
for (i in seq(1,11,2)){
	r2_tmp = getR2(chipseq_sig_mat_raw[used_row,i],chipseq_sig_mat_raw[used_row,i+1]/mean(chipseq_sig_mat_raw[,i+1])*mean(chipseq_sig_mat_raw[,11]/2+chipseq_sig_mat_raw[,12]/2))
	print(r2_tmp)
	r2_vec = c(r2_vec, r2_tmp)
}
r2_mat = cbind(r2_mat, r2_vec)

colnames(r2_mat) = c('RAW', 'qPCR_LMnorm', 'BGnorm')
rownames(r2_mat) = c('0A', '4A', '6A', '12A', '18A', '24A')


write.table(r2_mat, 'between_replicates_R2.txt', quote=F, sep='\t', col.names=T, row.names=T)

library(ggplot2)
library(RColorBrewer)


x_2c = c()
for (i in c(1,3,2)){
x1 = as.data.frame(rep(colnames(r2_mat)[i], dim(r2_mat)[1]))
x2 = r2_mat[,i]
x12 = cbind(x1,x2)
colnames(x12) = c('Method','R2')
x_2c = rbind(x_2c, x12)
}

### all R2
pdf('R2_boxplot.pdf', width=5, height=4)
x_2c$R2 = as.numeric(x_2c$R2)
x_2c$Method = factor(x_2c$Method, levels = unique(x_2c$Method),ordered = TRUE)
p = ggplot(data = x_2c, aes(x=Method, y=(R2)))
p = p + geom_boxplot(aes(fill = Method))
p = p + geom_point(aes(y=(R2), group=Method), position = position_dodge(width=0.1))
p = p + scale_fill_manual(values=rep(c('gray', 'orange1', 'dodgerblue1', 'dodgerblue'),each = 1))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
#p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
#p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10))
p = p + ylim(0, 1)
plot(p)
dev.off()

