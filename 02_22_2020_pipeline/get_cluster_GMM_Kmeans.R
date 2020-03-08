library(pheatmap)
library(mclust)
library(data.table)
library(mixtools)

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
et_mat_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.Ward.txt'
new_folder='LMqPCRnorm_GMM_Kmeans_folder'
plot_lim = 50

used_replicates = c(1:4, 5,6, 8,9, 10,11, 12,13, 15,17, 18,20)


### function
colMedian = function(x){
	return(apply(x, 2, median))
}

### get input mat
chipseq_sig_mat = as.data.frame(fread(signal_mat_file))
chipseq_sig_mat = chipseq_sig_mat[,used_replicates]
chipseq_Ward_mat = as.data.frame(fread(et_mat_file))
### na is all 0s

### get signal mat
used_row = (grepl("ctcf", chipseq_sig_mat[,4]))
chipseq_sig_mat_ctcf = chipseq_sig_mat[used_row,]
d1_0 = chipseq_sig_mat[used_row,-c(1:4)]
et_mat = chipseq_sig_mat[used_row,-c(1:4)]

### 
sig_0A = chipseq_sig_mat_ctcf[,5]/2+chipseq_sig_mat_ctcf[,6]/2
mixmdl = normalmixEM(sig_0A, k=3)
pdf(paste('all_pk/', new_folder, '/gmm_hist.pdf', sep=''))
plot(mixmdl,which=2)
lines(density(sig_0A), lty=2, lwd=2)
dev.off()

### with/across cluster distance
nr = 20
seed = 19
set.seed(seed)
var_within = c()
withinss_all = c()
for (nr_tmp in 2:nr){
	fit = kmeans(et_mat,centers=nr_tmp)
	withinss_all = c(withinss_all, sum(fit$withinss))
}
center0 = colMeans(et_mat)
dist0 = sum(apply(et_mat, 1, function(x) ((x-center0)^2)))
### plot with/across cluster distance
pdf(paste('all_pk/', new_folder, '/raw_kmeans.dist.clean.pdf', sep=''))
plot(1:nr, c(1, withinss_all/dist0))
lines(1:nr, c(1, withinss_all/dist0))
dev.off()
### get GMM predict
gmm_cluster = apply(mixmdl$posterior, 1, which.max)


### Kmeans
nr=9
nr_gi = nr/3
iter_num = 100
repeat_num = 50
set.seed(seed)

fit_cluster_reorder_mat = c()
d1_0_reorder = c()
chipseq_sig_mat_ctcf_reorder = c()
for (iter in 1:iter_num){
print(iter)
fit_cluster_reorder_vec_gi = c()
for (gi in 1:3){
gi_c_id = c(1:3)+3*(gi-1)
### get GMM i
et_mat_gi = et_mat[gmm_cluster==gi,]
d1_0_gi = d1_0[gmm_cluster==gi,]
if (iter==1){
d1_0_reorder = rbind(d1_0_reorder, d1_0_gi)
chipseq_sig_mat_ctcf_reorder_gi = chipseq_sig_mat_ctcf[gmm_cluster==gi,]
chipseq_sig_mat_ctcf_reorder = rbind(chipseq_sig_mat_ctcf_reorder, chipseq_sig_mat_ctcf_reorder_gi)
}
### cluster
fit = kmeans(et_mat_gi,centers=nr_gi)
fit_label = fit$cluster
### sort cluster
nr1= dim(table(fit$cluster))
fit_mean = c()
for (i in c(1:nr_gi)){
	sig_tmp = d1_0_gi[fit$cluster==i,]
	sig_tmp_mean = mean(as.matrix(sig_tmp))
	fit_mean[i] = sig_tmp_mean
}
plot_cluster_rank = order(fit_mean)
plot_cluster_rank
### relabel cluster label
fit_cluster_reorder = fit$cluster
for (i in 1:nr_gi){
	label = plot_cluster_rank[i]
	fit_cluster_reorder[fit$cluster==label] = gi_c_id[i]
}
fit_cluster_reorder_vec_gi = c(fit_cluster_reorder_vec_gi, fit_cluster_reorder)
}
fit_cluster_reorder_mat = cbind(fit_cluster_reorder_mat, fit_cluster_reorder_vec_gi)
}
colnames(fit_cluster_reorder_mat) = 1:iter_num


### shared cluster member
fit_cluster_reorder_vec = rep(0, dim(fit_cluster_reorder_mat)[1])
for (i in 1:nr){
used_label = i
fit_cluster_reorder_mat_binary = apply(fit_cluster_reorder_mat, 1, function(x) sum(x==used_label)>repeat_num)
fit_cluster_reorder_vec[fit_cluster_reorder_mat_binary] = used_label
}

### get new order
d1_0 = d1_0_reorder
chipseq_sig_mat_ctcf = chipseq_sig_mat_ctcf_reorder

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
	pk_k = chipseq_Ward_mat[fit_cluster_reorder_vec==i,-c(1:4)]
	dr_mclust_mat = rbind(dr_mclust_mat, pk_k)
	kmeans_center = rbind(kmeans_center, colMeans(pk_k[!is.na(pk_k[,1]),]))
}
ward_stat = dr_mclust_mat[,-c(1:4)]
ward_stat_c = cbind(dr_mclust_label[,1], ward_stat)
colnames(ward_stat_c)[1] = 'label'

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
dr_mclust_label_rescued = unlist(apply(chipseq_Ward_mat[rescue_id,-c(1:4)], 1, function(x) kmeans_rescue(x)))
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
	dr_mclust_plot = rbind(dr_mclust_plot, d1_0[fit_cluster_reorder_vec==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	pk_k = chipseq_sig_mat_ctcf[fit_cluster_reorder_vec==i,]
	#pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit_cluster_reorder_vec==i))))
	write.table(pk_k, paste('all_pk/', new_folder, '/raw_kmeans.', toString(k), '.pk.clean.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
	### plot sig boxplot
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0[fit_cluster_reorder_vec==i,]))
	sigmat = sigmat0#apply(sigmat0, 2, function(x) x-sigmat0[,1])
	boxplot(sigmat, ylim=c(0, plot_lim) ,outline=FALSE)
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	#
	lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
	dev.off()
	### plot FC
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.fc.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0[fit_cluster_reorder_vec==i,]))
	sigmat = apply(sigmat0, 2, function(x) (x*2/(sigmat0[,1]+sigmat0[,2])) )
	boxplot(sigmat, ylim=c(0, 2) ,outline=FALSE)
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
	dev.off()
}






i = 8
new_folder = 'LMqPCRnorm_folder'
d1 = read.table(paste('all_pk/LMqPCRnorm_folder/raw_kmeans.', i, '.pk.clean.txt', sep=''), header=T)
d2 = read.table(paste('all_pk/LMqPCRnorm_folder/raw_kmeans.', i+1, '.pk.clean.txt', sep=''), header=T)
sigmat0 = rbind(d1[,-c(1:4)], d2[,-c(1:4)])
sigmat = sigmat0

plot_lim = 50
pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', i, '_', i+1, '.clean.pdf', sep=''))
boxplot(sigmat, ylim=c(0, plot_lim) ,outline=FALSE)
if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
for (j in c(1:line_num)){
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

plot_lim = 50
pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', i, '_', i+1, '.clean.pdf', sep=''))
boxplot(sigmat, ylim=c(0, plot_lim) ,outline=FALSE)
if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
for (j in c(1:line_num)){
	lines(c(1:(dim(sigmat)[2])), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(c(1:(dim(sigmat)[2])), colMedian(sigmat), col='black')
dev.off()



