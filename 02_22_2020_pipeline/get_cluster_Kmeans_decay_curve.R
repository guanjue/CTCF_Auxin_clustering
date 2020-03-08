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
new_folder='LMqPCRnorm_decay_folder'
output_name = 'stable_peaks.txt'
plot_lim = 80
output_folder = ''

used_replicates = c(1:4, 5,6, 8,9, 10,11, 12,13, 15,17, 18,20)


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
tp = c(0,0,4,4,6,6,12,12,18,18,24,24)
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
return(c(A, a$coefficients[1], mean(a$residuals)))
}


### get decay curve coefficients
d1_0_OD_sig = d1_0_OD0[,-c(1:4)]
#d1_0_OD_sig = t(apply(d1_0_OD_sig, 1, function(x) (x+1)/(x[1]+x[2]+1)*50*2))
d1_0_OD_sig_decay0 = t(apply(as.matrix((d1_0_OD_sig)), 1, function(x) get_decay_curve(((log(x+1))), tp)))
d1_0_OD_sig_decay = d1_0_OD_sig_decay0[,1:2]

### get raw A(x)_vs_B(y)
png('A_vs_B.png')
heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2])
dev.off()

### get mean vs var
seq_mean_var = seq(0,max(d1_0_OD_sig_decay[,1]), length.out=1000)
mat_mean_var = c()
for (i in 1:(length(seq_mean_var)-1)){
used_id = (d1_0_OD_sig_decay[,1]>seq_mean_var[i]) & (d1_0_OD_sig_decay[,1]<=seq_mean_var[i+1])
d1_0_OD_sig_decay_B_tmp = d1_0_OD_sig_decay[used_id,2]
d1_0_OD_sig_decay_A_tmp = d1_0_OD_sig_decay[used_id,1]
print(mean(d1_0_OD_sig_decay_B_tmp))
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
pdf(paste('Bsd_vs_Amean.residuals.pdf', sep=''))
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
pdf('B.mean_vs_sd.pdf', width=5, height=5)
heatscatter(mat_mean_var_narm_sort[,1],mat_mean_var_narm_sort[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation')
abline(lmmodel2)
abline(v=Amean_lim)
dev.off()

### fit loess to decay parameter curve
d1_0_OD_sig_decay_df = as.data.frame(d1_0_OD_sig_decay)
span_win = 1
p = 1e-3
z = qnorm(1-p)
colnames(d1_0_OD_sig_decay_df) = c('A', 'B')
loessMod10 <- loess(B ~ A, data=d1_0_OD_sig_decay_df, span=span_win)
smoothed10 <- predict(loessMod10) 

pdf('A_vs_B.ci.pdf', width=5, height=5)
heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation')
d1_0_OD_sig_decay_df_A = d1_0_OD_sig_decay_df$A
bg_sd = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A
d1_0_OD_sig_decay_df_A_order = order(d1_0_OD_sig_decay_df_A)
lines(smoothed10[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="black")
lines(smoothed10[d1_0_OD_sig_decay_df_A_order]+z*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
lines(smoothed10[d1_0_OD_sig_decay_df_A_order]-z*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
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
pnew = 1e-3
z_new = qnorm(1-pnew)
d1_0_OD_sig_decay_df_A_all = d1_0_OD_sig_decay_df$A
bg_sd_all = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A_all
### binarize it 
smoothed10_alldata = predict(loessMod10_iter, newdata = d1_0_OD_sig_decay_df)
upperlims = smoothed10_alldata+z_new*bg_sd_all
stable_peak_binary = d1_0_OD_sig_decay_df$B > upperlims

### write stable peaks
d1_0_OD0_stable_peak = d1_0_OD0[stable_peak_binary,]
write.table(d1_0_OD0_stable_peak, output_name, quote=F, col.names=T, row.names=F, sep='\t')

pdf('A_vs_B.ci.iter.pdf', width=5, height=5)
heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2], xlab='0A mean signal (log)', ylab='-DR standard deviation')
d1_0_OD_sig_decay_df_A = d1_0_OD_sig_decay_df_iter$A
bg_sd = lmmodel2$coefficients[1] + lmmodel2$coefficients[2] * d1_0_OD_sig_decay_df_A
d1_0_OD_sig_decay_df_A_order = order(d1_0_OD_sig_decay_df_A)
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="black")
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order]+z_new*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
lines(smoothed10_iter[d1_0_OD_sig_decay_df_A_order]-z_new*bg_sd[d1_0_OD_sig_decay_df_A_order], x=d1_0_OD_sig_decay_df_A[d1_0_OD_sig_decay_df_A_order], col="cyan2")
dev.off()

















chr14:69,620,667-70,484,256






### Kmeans
nr=9
iter_num = 10
repeat_num = iter_num/2
set.seed(seed)

fit_cluster_reorder_mat0 = c()
for (iter in 1:iter_num){
	print(iter)
	### cluster
	#fit = Mclust(log(d1_0_OD_sig+1), nr)
	fit = Mclust(d1_0_OD_sig_decay, nr)
	#fit = kmeans((d1_0_OD_sig_decay), center = nr)
	fit_label = fit$classification
	### sort cluster
	nr1= dim(table(fit$classification))
	fit_mean = c()
	for (i in c(1:nr1)){
		sig_tmp = log(d1_0_OD_sig+1)[fit$classification==i,2]
		sig_tmp_mean = mean(as.matrix(sig_tmp))
		fit_mean[i] = sig_tmp_mean
	}
	plot_cluster_rank = order(fit_mean)
	plot_cluster_rank
	### relabel cluster label
	fit_cluster_reorder = fit$classification
	for (i in 1:nr){
		label = plot_cluster_rank[i]
		fit_cluster_reorder[fit$classification==label] = i
	}
	fit_cluster_reorder_mat0 = cbind(fit_cluster_reorder_mat0, fit_cluster_reorder)
}
colnames(fit_cluster_reorder_mat0) = 1:iter_num


##################
### align center
fit_cluster_forcenter = fit_cluster_reorder_mat0[,1]
Mclust_center = c()
k=0
for (i in 1:nr){
	k = k+1
	print(sum(fit_cluster_forcenter==i))
	pk_k = log(d1_0_OD_sig+1)[fit_cluster_forcenter==i,]
	Mclust_center = rbind(Mclust_center, colMeans(pk_k[!is.na(pk_k[,1]),]))
}

fit_cluster_reorder_mat_centeralign = c()
for (j in 1:iter_num){
fit_cluster_realign = fit_cluster_reorder_mat[,1]
fit_cluster_realign_new_label = fit_cluster_realign
k=0
for (i in 1:nr){
	k = k+1
	pk_k = log(d1_0_OD_sig+1)[fit_cluster_realign==i,]
	new_c_center = colMeans(pk_k[!is.na(pk_k[,1]),])
	dist_to_center = t(apply(Mclust_center, 1, function(x) sum((x-new_c_center)^2)))
	new_label = which.min(dist_to_center)
	fit_cluster_realign_new_label[fit_cluster_realign==i] = new_label
}
fit_cluster_reorder_mat_centeralign = cbind(fit_cluster_reorder_mat_centeralign, fit_cluster_realign_new_label)
}

fit_cluster_reorder_mat = fit_cluster_reorder_mat_centeralign


### shared cluster member
fit_cluster_reorder_vec = rep(0, dim(fit_cluster_reorder_mat)[1])
for (i in 1:nr){
	used_label = i
	fit_cluster_reorder_mat_binary = apply(fit_cluster_reorder_mat, 1, function(x) sum(x==used_label)>repeat_num)
	#fit_cluster_reorder_mat_num = apply(fit_cluster_reorder_mat, 1, function(x) sum(x==used_label))
	#print(summary(fit_cluster_reorder_mat_num))
	fit_cluster_reorder_vec[fit_cluster_reorder_mat_binary] = used_label
}

### get cluster peaks
dr_mclust_plot = c()
dr_mclust_plot_label = c()
dr_mclust_plot_cluster_num = c()
pk_k_all = c()
k=0

for (i in 1:nr){
	k = k+1
	if (sum(fit_cluster_reorder_vec==i)!=0){
	dr_mclust_plot_cluster_num = c(dr_mclust_plot_cluster_num, rep(k, sum(fit_cluster_reorder_vec==i)))
	print(sum(fit_cluster_reorder_vec==i))
	dr_mclust_plot = rbind(dr_mclust_plot, d1_0_OD[fit_cluster_reorder_vec==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	pk_k = d1_0_OD[fit_cluster_reorder_vec==i,]
	#pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit_cluster_reorder_vec==i))))
	write.table(pk_k, paste('all_pk/', new_folder, '/raw_kmeans.', toString(k), '.pk.clean.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
	### plot sig boxplot
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0_OD[fit_cluster_reorder_vec==i,-c(1:4)]))
	sigmat = sigmat0#apply(sigmat0, 2, function(x) x-sigmat0[,1])
	boxplot((sigmat), ylim=(c(0.1, plot_lim)) ,outline=FALSE, at =tp_plot, log='')
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(tp_plot, (sigmat[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	#
	lines(tp_plot, colMedian((sigmat)), col='black')
	dev.off()
}
}

i = 8
j = 9
sigmat0 = as.matrix((d1_0_OD[fit_cluster_reorder_vec==i,-c(1:4)]))
sigmat0_0A_8 = sigmat0[,1]/2+sigmat0[,2]/2

sigmat0_9 = as.matrix((d1_0_OD[fit_cluster_reorder_vec==j,-c(1:4)]))
sigmat0_0A_9 = sigmat0[,1]/2+sigmat0[,2]/2

sigmat0_9_modified = c()
for (i in 1:dim(sigmat0)[1]){
	tmp_sig = round(sigmat0[i,1]/2+sigmat0[i,2]/2)
	tmp = sigmat0_9[round(sigmat0_9[,1]/2+sigmat0_9[,2]/2)==tmp_sig,]
	sigmat0_9_modified = rbind(sigmat0_9_modified, tmp)
}

i = 9
pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', toString(i), '.modified.clean.pdf', sep=''))
sigmat = sigmat0_9_modified
boxplot((sigmat), ylim=(c(0.1, plot_lim)) ,outline=FALSE, at =tp_plot, log='y')
if (dim(sigmat)[1]>1000){line_num = sample(dim(sigmat)[1], 1000)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
	lines(tp_plot, (sigmat[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(tp_plot, colMedian((sigmat)), col='black')
dev.off()












d1_0_OD_sig_decay_df = as.data.frame(d1_0_OD_sig_decay)
colnames(d1_0_OD_sig_decay_df) = c('A', 'B')
loessMod10 <- loess(B ~ A, data=d1_0_OD_sig_decay_df, span=0.10)
smoothed10 <- predict(loessMod10) 
reg1 <- lm(B ~ A,data=d1_0_OD_sig_decay_df[d1_0_OD_sig_decay_df[,1]<1.5,]) 


png('hist.sig.A.png')
hist(d1_0_OD_sig_decay[d1_0_OD_sig_decay[,1]>2,1], breaks=50, xlim=c(0,4.5))
dev.off()

png('hist.sig.B.png')
hist(d1_0_OD_sig_decay[d1_0_OD_sig_decay[,1]>2,2], breaks=50, xlim=c(-0.2,0.05))
dev.off()

png('hist.sig.A.all.png')
hist(d1_0_OD_sig_decay[d1_0_OD_sig_decay[,1]>-200,1], breaks=50, xlim=c(0,4.5))
dev.off()

png('hist.sig.B.all.png')
hist(d1_0_OD_sig_decay[d1_0_OD_sig_decay[,1]>-200,2], breaks=50, xlim=c(-0.2,0.05))
dev.off()




png('A_vs_B.ci.1.png')
library(LSD)
#heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2]/smoothed10, ylim=c(-5,5))
heatscatter(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2])
lines(smoothed10[order(d1_0_OD_sig_decay_df$A)], x=d1_0_OD_sig_decay_df$A[order(d1_0_OD_sig_decay_df$A)], col="red")

#lines(smoothed10[order(d1_0_OD_sig_decay_df$A)], x=d1_0_OD_sig_decay_df$A[order(d1_0_OD_sig_decay_df$A)], col="red")
#abline(reg1)
dev.off()


set.seed(2019)
AB_mclust = Mclust(d1_0_OD_sig_decay[,], 2)
png('A_vs_B.gmm.png')
plot(AB_mclust, what = "classification", main = "Mclust Classification")
dev.off()

mclust_1r = apply(AB_mclust$z, 1, function(x) which.max(x))
set.seed(2019)
AB_mclust = Mclust(d1_0_OD_sig_decay[d1_0_OD_sig_decay[,1]>2,], 3)
png('A_vs_B.gmm.2.png')
plot(AB_mclust, what = "classification", main = "Mclust Classification", xlim=c(0,4.5),ylim=c(-0.2,0.15))
dev.off()




set.seed(2019)
AB_mclust = Mclust(d1_0_OD_sig_decay[,], 4)
png('A_vs_B.gmm.png')
plot(AB_mclust, what = "classification", main = "Mclust Classification")
dev.off()


colors = c('black','gray', 'purple','blue', 'cyan', 'green', 'orange', 'red', 'brown')
png('A_vs_B.Mclust.png')
library(LSD)
plot(d1_0_OD_sig_decay[,1],d1_0_OD_sig_decay[,2])
for (i in 1:9){
points(d1_0_OD_sig_decay[fit_cluster_reorder_vec==i,1],d1_0_OD_sig_decay[fit_cluster_reorder_vec==i,2], col=colors[i])
}
dev.off()






### get cluster peaks
dr_mclust_plot = c()
dr_mclust_plot_label = c()
dr_mclust_plot_cluster_num = c()
pk_k_all = c()
k=0

for (i in 1:nr){
	k = k+1
	dr_mclust_plot_cluster_num = c(dr_mclust_plot_cluster_num, rep(k, sum(fit_cluster_reorder_vec==i)))
	print(sum(fit_cluster_reorder_vec==i))
	dr_mclust_plot = rbind(dr_mclust_plot, d1_0_OD[fit_cluster_reorder_vec==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	pk_k = d1_0_OD[fit_cluster_reorder_vec==i,]
	#pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit_cluster_reorder_vec==i))))
	write.table(pk_k, paste('all_pk/', new_folder, '/raw_kmeans.', toString(k), '.pk.clean.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
	### plot sig boxplot
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0_OD[fit_cluster_reorder_vec==i,-c(1:4)]))
	sigmat = sigmat0#apply(sigmat0, 2, function(x) x-sigmat0[,1])
	boxplot(log(sigmat), ylim=log(c(1.1, plot_lim)) ,outline=FALSE, at =tp_plot, log='y')
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(tp_plot, log(sigmat[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	#
	lines(tp_plot, colMedian(log(sigmat)), col='black')
	dev.off()
}










##################
### peak rescue
fit_cluster_reorder_vec0 = fit_cluster_reorder_vec
#fit_cluster_reorder_vec = fit_cluster_reorder_vec0
dr_mclust_mat = c()
dr_mclust_label = c()
kmeans_center = c()
k=0
for (i in 1:nr){
	k = k+1
	print(sum(fit_cluster_reorder_vec==i))
	dr_mclust_label = rbind(dr_mclust_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	#pk_k = chipseq_sig_mat_ctcf[fit_cluster_reorder_vec==i,]
	pk_k = d1_0[fit_cluster_reorder_vec==i,-c(1:4)]
	dr_mclust_mat = rbind(dr_mclust_mat, pk_k)
	kmeans_center = rbind(kmeans_center, colMeans(pk_k[!is.na(pk_k[,1]),]))
}

### rescuing
kmeans_rescue = function(xnew){
	if(!is.na(xnew[1])){
		ed = apply(kmeans_center, 1, function(x) sqrt(sum((x-xnew)^2)))
		cn = which.min(ed)
	} else{
		cn = 1
	}
	return(cn)
}
rescue_id = fit_cluster_reorder_vec==0
dr_mclust_label_rescued = unlist(apply(d1_0[rescue_id,-c(1:4)], 1, function(x) kmeans_rescue(x)))
fit_cluster_reorder_vec[rescue_id] = dr_mclust_label_rescued
table(fit_cluster_reorder_vec)
table(fit_cluster_reorder_vec0)

### get cluster peaks
dr_mclust_plot = c()
dr_mclust_plot_label = c()
dr_mclust_plot_cluster_num = c()
pk_k_all = c()
k=0

for (i in 1:nr){
	k = k+1
	dr_mclust_plot_cluster_num = c(dr_mclust_plot_cluster_num, rep(k, sum(fit_cluster_reorder_vec==i)))
	print(sum(fit_cluster_reorder_vec==i))
	dr_mclust_plot = rbind(dr_mclust_plot, d1_0_OD[fit_cluster_reorder_vec==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	pk_k = d1_0_OD[fit_cluster_reorder_vec==i,]
	#pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit_cluster_reorder_vec==i))))
	write.table(pk_k, paste('all_pk/', new_folder, '/raw_kmeans.', toString(k), '.pk.clean.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
	### plot sig boxplot
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0_OD[fit_cluster_reorder_vec==i,-c(1:4)]))
	sigmat = sigmat0#apply(sigmat0, 2, function(x) x-sigmat0[,1])
	boxplot(sigmat, ylim=c(0.1, plot_lim) ,outline=FALSE, at =c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5), log='')
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	#
	lines(c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5), colMedian(sigmat), col='black')
	dev.off()
}





i = 8
new_folder = 'LMqPCRnorm_folder'
d1 = read.table(paste('all_pk/LMqPCRnorm_folder/raw_kmeans.', i, '.pk.clean.txt', sep=''), header=T)
d2 = read.table(paste('all_pk/LMqPCRnorm_folder/raw_kmeans.', i+1, '.pk.clean.txt', sep=''), header=T)
sigmat0 = rbind(d1[,-c(1:4)], d2[,-c(1:4)])
sigmat = sigmat0

plot_lim = 80
pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', i, '_', i+1, '.clean.pdf', sep=''))
boxplot(sigmat, ylim=c(0, plot_lim) ,outline=FALSE)
if (dim(sigmat)[1]>1000){used_line = sample(dim(sigmat)[1], 1000)} else {used_line=1:dim(sigmat)[1]}
for (j in used_line){
	lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
dev.off()


i = 6
new_folder = 'LMqPCRnorm_folder'
d1 = read.table(paste('all_pk/LMqPCRnorm_folder/raw_kmeans.', i, '.pk.clean.txt', sep=''), header=T)
d2 = read.table(paste('all_pk/LMqPCRnorm_folder/raw_kmeans.', i+1, '.pk.clean.txt', sep=''), header=T)
sigmat0 = rbind(d1[,-c(1:4)], d2[,-c(1:4)])
sigmat = sigmat0

plot_lim = 80
pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', i, '_', i+1, '.clean.pdf', sep=''))
boxplot(sigmat, ylim=c(0, plot_lim) ,outline=FALSE)
if (dim(sigmat)[1]>1000){used_line = sample(dim(sigmat)[1], 1000)} else {used_line=1:dim(sigmat)[1]}
for (j in used_line){
	lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
dev.off()


