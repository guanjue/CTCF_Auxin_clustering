library(pheatmap)
library(mclust)
library(data.table)

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
new_folder='LMqPCRnorm_folder'
plot_lim = 80

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
d1_0_OD = chipseq_sig_mat[order(chipseq_sig_mat_ctcf[,4]),]
et_mat_OD = chipseq_Ward_mat[order(chipseq_Ward_mat[,4]),-c(1:4)]

### rm NA
rm_na = !is.na(et_mat_OD[,1])
d1_0 = d1_0_OD[rm_na,]
et_mat = et_mat_OD[rm_na,]

### add missing time by mean
#et_mat2 = d1_0_OD[,1]/2
#et_mat4 = d1_0_OD[,1]
#et_mat6 = d1_0_OD[,2]
#et_mat8 = (d1_0_OD[,3]-d1_0_OD[,2])/3*1 + d1_0_OD[,2]
#et_mat10 = (d1_0_OD[,3]-d1_0_OD[,2])/3*2 + d1_0_OD[,2]
#et_mat12 = d1_0_OD[,3]
#et_mat14 = (d1_0_OD[,4]-d1_0_OD[,3])/3*2 + d1_0_OD[,2]
#et_mat16 = (d1_0_OD[,4]-d1_0_OD[,3])/3*2 + d1_0_OD[,2]
#et_mat18 = d1_0_OD[,4]
#et_mat20 = (d1_0_OD[,5]-d1_0_OD[,4])/3*2 + d1_0_OD[,2]
#et_mat22 = (d1_0_OD[,5]-d1_0_OD[,4])/3*2 + d1_0_OD[,2]
#et_mat24 = d1_0_OD[,5]
#et_mat = cbind(et_mat2, et_mat4, et_mat6, et_mat8, et_mat10, et_mat12, et_mat14, et_mat16, et_mat18, et_mat20, et_mat22, et_mat24)

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

### Kmeans
nr=9
iter_num = 50
repeat_num = iter_num/2
set.seed(seed)

fit_cluster_reorder_mat = c()
for (iter in 1:iter_num){
print(iter)
### cluster
fit = kmeans((et_mat),centers=nr)
fit_label = fit$cluster
### sort cluster
nr1= dim(table(fit$cluster))
fit_mean = c()
for (i in c(1:nr1)){
	sig_tmp = d1_0[fit$cluster==i,5:6]
	sig_tmp_mean = mean(as.matrix(sig_tmp))
	fit_mean[i] = sig_tmp_mean
}
plot_cluster_rank = order(fit_mean)
plot_cluster_rank
### relabel cluster label
fit_cluster_reorder = fit$cluster
for (i in 1:nr){
	label = plot_cluster_rank[i]
	fit_cluster_reorder[fit$cluster==label] = i
}
fit_cluster_reorder_mat = cbind(fit_cluster_reorder_mat, fit_cluster_reorder)
}
colnames(fit_cluster_reorder_mat) = 1:iter_num

### shared cluster member
fit_cluster_reorder_vec = rep(0, dim(fit_cluster_reorder_mat)[1])
for (i in 1:nr){
used_label = i
fit_cluster_reorder_mat_binary = apply(fit_cluster_reorder_mat, 1, function(x) sum(x==used_label)>repeat_num)
fit_cluster_reorder_vec[fit_cluster_reorder_mat_binary] = used_label
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
	dr_mclust_plot = rbind(dr_mclust_plot, d1_0[fit_cluster_reorder_vec==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	pk_k = d1_0[fit_cluster_reorder_vec==i,]
	#pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit_cluster_reorder_vec==i))))
	write.table(pk_k, paste('all_pk/', new_folder, '/raw_kmeans.', toString(k), '.pk.clean.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
	### plot sig boxplot
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0[fit_cluster_reorder_vec==i,-c(1:4)]))
	sigmat = sigmat0#apply(sigmat0, 2, function(x) x-sigmat0[,1])
	boxplot(sigmat, ylim=c(0.1, plot_lim) ,outline=FALSE, at =c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5), log='')
	if (dim(sigmat)[1]>1000){line_num = 1000} else {line_num=dim(sigmat)[1]}
	for (j in c(1:line_num)){
		lines(c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5), sigmat[j,], col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	#
	lines(c(0,0.5, 4,4.5, 6,6.5, 12,12.5, 18,18.5, 24,24.5), colMedian(sigmat), col='black')
	dev.off()
	### plot FC
	pdf(paste('all_pk/', new_folder, '/raw_kmeans.box.fc.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0[fit_cluster_reorder_vec==i,-c(1:4)]))
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


