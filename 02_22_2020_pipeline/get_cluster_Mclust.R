library(pheatmap)
library(mclust)
library(data.table)
library(MASS)

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
new_folder='LMqPCRnorm_folder_Mclust'
plot_lim = 50

used_replicates = c(1:4, 5,6, 8,9, 10,11, 12,13, 15,17, 18,20)

tp = c(0,0,4,4,6,6,12,12,18,18,24,24)

### function
colMedian = function(x){
	return(apply(x, 2, median))
}


### get input mat
chipseq_sig_mat = as.data.frame(fread(signal_mat_file))
chipseq_sig_mat = chipseq_sig_mat[,used_replicates]
chipseq_Ward_mat = as.data.frame(fread(et_mat_file))

### get signal mat
used_row = (grepl("ctcf", chipseq_sig_mat[,4]))
chipseq_sig_mat_ctcf0 = chipseq_sig_mat[used_row,]
d1_00 = chipseq_sig_mat[used_row,-c(1:4)]
et_mat0 = as.matrix(chipseq_Ward_mat[used_row,-c(1:4)])


used_id_narm = !is.na(et_mat0[,1])
chipseq_sig_mat_ctcf = chipseq_sig_mat_ctcf0[used_id_narm,]
d1_0 = d1_00[used_id_narm,]
et_mat = et_mat0[used_id_narm,]

chipseq_sig_mat_ctcf_FPpk = chipseq_sig_mat_ctcf0[!used_id_narm,]
write.table(chipseq_sig_mat_ctcf_FPpk, 'FPpk.txt', quote=F, col.names=F, row.names=F, sep='\t')

### with/across cluster distance
nr = 20
seed = 1

### Kmeans
iter_num = 30
repeat_num = 15
set.seed(seed)

model_list = c("EII", "VII", "EEI", "EVI", "VEI", "VVI", 'EEE', 'EVE', 'VEE', 'VVE', 'EEV', 'VEV', 'EVV', 'VVV')
set.seed(seed)

BIC = mclustBIC(et_mat, modelNames=model_list, prior = priorControl())
pdf(paste('all_pk/', new_folder, '/all_sig_hist.pdf', sep=''))
plot(BIC)
dev.off()

fit0 = Mclust(et_mat, x = BIC)
fit0_label = fit0$classification
nr0 = length(unique(fit0_label))


### Mclust
fit_cluster_reorder_mat = c()
fit_label = fit0$classification
### sort cluster
nr1= dim(table(fit0$classification))
fit_mean = c()
for (i in c(1:nr1)){
	sig_tmp = d1_0[fit0$classification==i,]
	#fit_mean = rbind(fit_mean, c(mean(sig_tmp[,1]/2+sig_tmp[,2]/2)))
	fit_mean = rbind(fit_mean, mean(sig_tmp[!is.na(sig_tmp)]))
	#fit_mean = c(fit_mean, mean(decay_tmp[,1]))
}
plot_cluster_rank = order(fit_mean)
plot_cluster_rank
### relabel cluster label
fit_cluster_reorder = fit0$classification
for (i in 1:nr){
	label = plot_cluster_rank[i]
	fit_cluster_reorder[fit0$classification==label] = i
}
fit_cluster_reorder_mat = cbind(fit_cluster_reorder_mat, fit_cluster_reorder)


for (iter in 1:iter_num){
print(iter)
### cluster
fit = Mclust(et_mat, G = nr1)
fit_label = fit$classification
### sort cluster
nr1= dim(table(fit$classification))
fit_mean = c()
for (i in c(1:nr1)){
	sig_tmp = d1_0[fit$classification==i,]
	fit_mean = rbind(fit_mean, c(mean(sig_tmp[,1]/2+sig_tmp[,2]/2)))
	#fit_mean = rbind(fit_mean, mean(sig_tmp[!is.na(sig_tmp)]))
	#fit_mean = c(fit_mean, mean(decay_tmp[,1]))
}
plot_cluster_rank = order(fit_mean)
plot_cluster_rank
### relabel cluster label
fit_cluster_reorder = fit$classification
for (i in 1:nr1){
	label = plot_cluster_rank[i]
	fit_cluster_reorder[fit$classification==label] = i
}
fit_cluster_reorder_mat = cbind(fit_cluster_reorder_mat, fit_cluster_reorder)
}
colnames(fit_cluster_reorder_mat) = 1:(iter_num+1)

repeat_num = as.integer(dim(fit_cluster_reorder_mat)[2]/2)

### shared cluster member
fit_cluster_reorder_vec = rep(0, dim(fit_cluster_reorder_mat)[1])
for (i in 1:nr1){
used_label = i
fit_cluster_reorder_mat_binary = apply(fit_cluster_reorder_mat, 1, function(x) sum(x==used_label)>repeat_num)
fit_cluster_reorder_vec[fit_cluster_reorder_mat_binary] = used_label
}
#fit_cluster_reorder_vec = fit_cluster_reorder_mat[,1]
### rescue
#qda_model = qda(et_mat[fit_cluster_reorder_vec!=0,], fit_cluster_reorder_vec[fit_cluster_reorder_vec!=0])
#fit_cluster_reorder_vec_after_rescue = predict(qda_model,et_mat[fit_cluster_reorder_vec==0,])$class
#fit_cluster_reorder_vec[fit_cluster_reorder_vec==0] = fit_cluster_reorder_vec_after_rescue
##################
### peak rescue
fit_cluster_reorder_vec0 = fit_cluster_reorder_vec
#fit_cluster_reorder_vec = fit_cluster_reorder_vec0
kmeans_center = c()
k=0
for (i in 1:nr0){
	k = k+1
	print(sum(fit_cluster_reorder_vec==i))
	#pk_k = chipseq_sig_mat_ctcf[fit_cluster_reorder_vec==i,]
	pk_k = d1_0[fit_cluster_reorder_vec==i,]
	kmeans_center = rbind(kmeans_center, colMeans(log(pk_k[!is.na(pk_k[,1]),]+1)))
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
dr_mclust_label_rescued = unlist(apply(d1_0[rescue_id,], 1, function(x) kmeans_rescue(log(x+1))))
fit_cluster_reorder_vec[rescue_id] = dr_mclust_label_rescued
table(fit_cluster_reorder_vec)
table(fit_cluster_reorder_vec0)


fit_mean = c()
for (i in c(1:nr0)){
	sig_tmp = d1_0[fit_cluster_reorder_vec==i,]
	#fit_mean = rbind(fit_mean, c(mean(sig_tmp[,1]/2+sig_tmp[,2]/2)))
	fit_mean = rbind(fit_mean, mean(sig_tmp[!is.na(sig_tmp)]))
	#fit_mean = c(fit_mean, mean(decay_tmp[,1]))
}

plot_cluster_rank = order(fit_mean)
plot_cluster_rank




### get cluster peaks
dr_mclust_plot = c()
dr_mclust_plot_label = c()
dr_mclust_plot_cluster_num = c()
pk_k_all = c()
k=0

tp_plot = c(-0.25,0.25, 3.75,4.25, 5.75,6.25, 11.75,12.25, 17.25,18.25, 23.75,24.25)
tp_plot_forcurve = c(0,4,6,12,18,24)

table(fit_cluster_reorder_vec)


for (i in plot_cluster_rank){
	if (sum(fit_cluster_reorder_vec==i)>1){
	k = k+1
	dr_mclust_plot_cluster_num = c(dr_mclust_plot_cluster_num, rep(k, sum(fit_cluster_reorder_vec==i)))
	print(sum(fit_cluster_reorder_vec==i))
	dr_mclust_plot = rbind(dr_mclust_plot, d1_0[fit_cluster_reorder_vec==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit_cluster_reorder_vec==i)), rep(i, sum(fit_cluster_reorder_vec==i))) )
	#
	pk_k = chipseq_sig_mat_ctcf[fit_cluster_reorder_vec==i,]
	#pk_k_all = rbind(pk_k_all, cbind(pk_k, rep(k, sum(fit_cluster_reorder_vec==i))))
	write.table(pk_k, paste('all_pk/', new_folder, '/raw_mclust.', toString(k), '.pk.clean.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
	### plot sig boxplot
	pdf(paste('all_pk/', new_folder, '/raw_mclust.box.', toString(k), '.clean.pdf', sep=''))
	sigmat0 = as.matrix((d1_0[fit_cluster_reorder_vec==i,]))
	sigmat = sigmat0#apply(sigmat0, 2, function(x) x-sigmat0[,1])
	sigmat = sigmat0+1
	sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))
	boxplot(sigmat, ylim=c(1, 60) ,outline=FALSE, at =tp_plot, log='y')
	if (dim(sigmat)[1]>2000){line_num = sample(dim(sigmat)[1], 2000)} else {line_num=1:dim(sigmat)[1]}
	for (j in (line_num)){
		lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
	}
	dev.off()
	}
}












