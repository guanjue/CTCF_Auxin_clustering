library(pheatmap)
library(mclust)
library(data.table)
library(mixtools)
library(LSD)
library(changepoint)

### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_mat_file = args[1]
et_mat_file = args[2]
new_folder = args[3]
plot_lim = as.numeric(args[4])

#signal_mat_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt'
#et_mat_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.Ward.txt'
#new_folder='raw_folder'
#plot_lim = 80

signal_mat_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.txt'
new_folder='/Users/universe/Downloads/'
output_name = 'stable_peaks.txt'
pnew = 1e-4
used_replicates = c(1:4, 5,6, 8,9, 10,11, 12,13, 15,17, 18,20)
tp = c(0,0,4,4,6,6,12,12,18,18,24,24)

### function
colMedian = function(x){
	return(apply(x, 2, median))
}

### get input mat
chipseq_sig_mat = as.data.frame(fread(signal_mat_file))
chipseq_sig_mat = chipseq_sig_mat[,used_replicates]
### na is all 0s

### get signal mat
used_row = (grepl("ctcf", chipseq_sig_mat[,4]))
chipseq_sig_mat_ctcf = chipseq_sig_mat[used_row,]
d1_0_OD0 = chipseq_sig_mat[order(chipseq_sig_mat_ctcf[,4]),]
tp_plot = c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5)

### get decay curve
get_decay_curve_pre = function(y,x){
a = lm(y~x)
print(summary(a))
print(summary(a$coefficients))
return(c(a$coefficients[1], a$coefficients[2], mean(a$residuals)))
}

get_decay_curve = function(y,x){
A = (y[1]+y[2])/2
y_sub = y - A
a = lm(y_sub~x-1)
#print(summary(a))
#print(summary(a$coefficients))
return(c(A, a$coefficients[1], summary(a)$r.squared))
}

get_linear_curve = function(y,x){
a = lm(y~x)
#print(summary(a))
#print(summary(a$coefficients))
return(summary(a)$r.squared)
}



### get decay curve coefficients
d1_0_OD_sig = d1_0_OD0[,-c(1:4)]
#d1_0_OD_sig = t(apply(d1_0_OD_sig, 1, function(x) (x+1)/(x[1]+x[2]+1)*50*2))
d1_0_OD_sig_decay0 = t(apply(as.matrix((d1_0_OD_sig)), 1, function(x) get_decay_curve(((log(x+1))), tp)))
d1_0_OD_sig_decay = d1_0_OD_sig_decay0[,1:2]
d1_0_OD_sig_decay_residual = d1_0_OD_sig_decay0[,3]

set.seed(2019)
mod1 <- Mclust(d1_0_OD_sig_decay, G = 3)
png('test.png', width=1000, height=1000)
plot(mod1, what = "classification")
dev.off()

mod1_mean = c()
for (i in 1:3){
	print(i)
	mean_tmp = mean(d1_0_OD_sig_decay0[mod1$classification==i,1])
	print(mean_tmp)
	mod1_mean = c(mod1_mean, mean_tmp)
}


class_num = order(mod1_mean)
kk = 0
for (i in class_num){
	kk = kk+1
	data_tmp = d1_0_OD0[mod1$classification==i,]
	write.table(data_tmp, paste('decay_curve.', kk,'.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
}



