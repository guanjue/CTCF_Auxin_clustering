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
new_folder='/storage/home/gzx103/scratch/ctcf_auxin/all_pk/LMqPCRnorm_decay_folder/'
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

d1_0_OD_sig_linear_residual = t(apply(as.matrix((d1_0_OD_sig)), 1, function(x) get_linear_curve(x, tp)))

### get mean vs var
seq_mean_var = seq(0,max(d1_0_OD_sig_decay[,1]), length.out=1000)
mat_mean_var = c()
for (i in 1:(length(seq_mean_var)-1)){
used_id = (d1_0_OD_sig_decay[,1]>seq_mean_var[i]) & (d1_0_OD_sig_decay[,1]<=seq_mean_var[i+1])
d1_0_OD_sig_decay_B_tmp = d1_0_OD_sig_decay[used_id,2]
d1_0_OD_sig_decay_A_tmp = d1_0_OD_sig_decay[used_id,1]
#print(mean(d1_0_OD_sig_decay_B_tmp))
mat_mean_var_tmp = c(mean(d1_0_OD_sig_decay_A_tmp), sd(d1_0_OD_sig_decay_B_tmp))
mat_mean_var = rbind(mat_mean_var, mat_mean_var_tmp)
}

### get changing point
# rm na
mat_mean_var_narm = mat_mean_var[(!is.na(rowSums(mat_mean_var))),]
# sort by A mean (time series)
mat_mean_var_narm_sort = mat_mean_var_narm[order(mat_mean_var_narm[,1]),]
# get change point thresh
used_id1 = (mat_mean_var_narm_sort[,1]<2) #& (mat_mean_var_narm_sort[,1]>1)
lmmodel1 = lm(mat_mean_var_narm_sort[used_id1,2]~mat_mean_var_narm_sort[used_id1,1])
residuals_all = mat_mean_var_narm_sort[,2] - mat_mean_var_narm_sort[,1]*lmmodel1$coefficients[2] - lmmodel1$coefficients[1]

### split residuals
set.seed(2019)
mixmdl = normalmixEM(residuals_all, k=2)
pdf(paste(new_folder, 'Bsd_vs_Amean.residuals.pdf', sep=''))
plot(mixmdl,which=2)
lines(density(residuals_all), lty=2, lwd=2)
dev.off()

### get Amean lim
gmm_posterior = mixmdl$posterior
print('mixmdl$mu')
print(mixmdl$mu)
used_id = (gmm_posterior[,1]>0.5) & (mat_mean_var_narm_sort[,1]<3)
Amean_lim = max(mat_mean_var_narm_sort[used_id,1])

# get change point thresh
used_id2 = (mat_mean_var_narm_sort[,1]<Amean_lim) & (mat_mean_var_narm_sort[,1]>1)
lmmodel2 = lm(mat_mean_var_narm_sort[used_id1,2]~mat_mean_var_narm_sort[used_id1,1])

### get mean vs sd linear relationships
pdf(paste(new_folder, 'B.mean_vs_sd.pdf', sep=''), width=5, height=5)
heatscatter(mat_mean_var_narm_sort[,1],mat_mean_var_narm_sort[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation')
abline(lmmodel2)
abline(v=Amean_lim)
dev.off()

### fit loess to decay parameter curve
d1_0_OD_sig_decay_df = as.data.frame(d1_0_OD_sig_decay)
p = 1e-3
z = qnorm(1-p)
span_win = 1
colnames(d1_0_OD_sig_decay_df) = c('A', 'B')
loessMod10 <- loess(B ~ A, data=d1_0_OD_sig_decay_df, span=span_win)
smoothed10 <- predict(loessMod10) 

pdf(paste(new_folder, 'A_vs_B.ci.pdf', sep=''), width=5, height=5)
heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation')
d1_0_OD_sig_decay_df_A = d1_0_OD_sig_decay_df$A
bg_sd = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A
d1_0_OD_sig_decay_df_A_order = order(d1_0_OD_sig_decay_df_A)
lines(smoothed10[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="black")
lines(smoothed10[d1_0_OD_sig_decay_df_A_order]+z*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
lines(smoothed10[d1_0_OD_sig_decay_df_A_order]-z*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
dev.off()


pdf(paste(new_folder, 'A_vs_B.ci.residual.hist.pdf', sep=''), width=5, height=5)
hist(d1_0_OD_sig_decay_residual, breaks=50)
dev.off()
pdf(paste(new_folder, 'A_vs_B.linear.residual.hist.pdf', sep=''), width=5, height=5)
hist(d1_0_OD_sig_linear_residual, breaks=50)
dev.off()





### iterative change loess
d1_0_OD_sig_decay_df_iter = d1_0_OD_sig_decay_df
smoothed10_iter = smoothed10
loessMod10_iter = loessMod10
# initialize pk number
pk_pre = dim(d1_0_OD_sig_decay_df_iter)[1]
### remove peaks with different dynamics
for (i in 1:100){
	print(i)
	smoothed10_exp_up = smoothed10_iter+z*bg_sd
	smoothed10_exp_down = smoothed10_iter-z*bg_sd
	pk_pre = length(smoothed10_exp_up)
	d1_0_OD_sig_decay_df_A = d1_0_OD_sig_decay_df_iter$A
	d1_0_OD_sig_decay_df_B = d1_0_OD_sig_decay_df_iter$B
	used_id = (d1_0_OD_sig_decay_df_B<=smoothed10_exp_up) | (d1_0_OD_sig_decay_df_A==max(d1_0_OD_sig_decay_df_A)) | (d1_0_OD_sig_decay_df_A==min(d1_0_OD_sig_decay_df_A)) #& (d1_0_OD_sig_decay_df_B>=smoothed10_exp_down)
	d1_0_OD_sig_decay_df_iter = d1_0_OD_sig_decay_df_iter[used_id,]
	bg_sd = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_iter$A
	print(dim(d1_0_OD_sig_decay_df_iter))
	loessMod10_iter <- loess(B ~ A, data=d1_0_OD_sig_decay_df_iter, span=span_win)
	smoothed10_iter <- predict(loessMod10_iter)
	pk_new = dim(d1_0_OD_sig_decay_df_iter)[1]
	if ((pk_pre-pk_new)<1){
		break
	} 
}

### binarize CTCF by iterative loess curve
z_new = qnorm(1-pnew)
d1_0_OD_sig_decay_df_A_all = d1_0_OD_sig_decay_df$A
bg_sd_all = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A_all
### binarize it 
smoothed10_alldata = predict(loessMod10_iter, newdata = d1_0_OD_sig_decay_df)
upperlims = smoothed10_alldata+z_new*bg_sd_all
stable_peak_binary = d1_0_OD_sig_decay_df$B > upperlims

### write stable peaks
d1_0_OD0_stable_peak = d1_0_OD0[stable_peak_binary,]
write.table(d1_0_OD0_stable_peak, paste(new_folder, output_name, sep=''), quote=F, col.names=T, row.names=F, sep='\t')
d1_0_OD0_unstable_peak = d1_0_OD0[!stable_peak_binary,]
write.table(d1_0_OD0_unstable_peak, paste(new_folder, output_name, '.unstable.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')

pdf(paste(new_folder, 'A_vs_B.ci.iter.pdf', sep=''), width=5, height=5)
heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation')
d1_0_OD_sig_decay_df_A = d1_0_OD_sig_decay_df_iter$A
bg_sd = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A
d1_0_OD_sig_decay_df_A_order = order(d1_0_OD_sig_decay_df_A)
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="black")
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order]+z_new*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order]-z_new*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
dev.off()

rbPal <- colorRampPalette(c('black','red'))
Col_all <- rbPal(101)[as.numeric(cut(d1_0_OD_sig_decay_residual,breaks = 100))]
pdf(paste(new_folder, 'A_vs_B.ci.iter.r2.pdf', sep=''), width=5, height=5)
plot(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation', col=Col_all, pch=20)
d1_0_OD_sig_decay_df_A = d1_0_OD_sig_decay_df_iter$A
bg_sd = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A
d1_0_OD_sig_decay_df_A_order = order(d1_0_OD_sig_decay_df_A)
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="black")
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order]+z_new*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order]-z_new*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
dev.off()



