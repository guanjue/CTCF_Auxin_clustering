library(mixtools)

### function
colMedian = function(x){
	return(apply(x, 2, median))
}

d1 = read.table('all_pk/LMqPCRnorm_decay_folder/stable_peaks.txt.unstable.txt', header=T)
d2 = read.table('all_pk/LMqPCRnorm_decay_folder/stable_peaks.txt', header=T)
new_folder='/storage/home/gzx103/scratch/ctcf_auxin/all_pk/LMqPCRnorm_decay_folder/'
output_name = 'balanced_peaklist_for_compare'

d10A = d1[,5]/2+d1[,6]/2
d20A = d2[,5]/2+d2[,6]/2


###### remove relatively low signals
### unstable 0A signal
all_0A_signal = c(d10A)
### split residuals
set.seed(2019)
mixmdl = normalmixEM(log(all_0A_signal+1), k=2)
pdf(paste(new_folder, 'extract_peaks_for_compare.signal.hist.0A.b.pdf', sep=''))
plot(mixmdl,which=2)
lines(density(log(all_0A_signal+1)), lty=2, lwd=2)
dev.off()

### get Amean lim
gmm_posterior = mixmdl$posterior
print('mixmdl$mu')
print(mixmdl$mu)
used_id = (gmm_posterior[,2]>0.5)
Amean_lim = min((all_0A_signal)[used_id])  #min(log(all_0A_signal+1)[used_id])
print(Amean_lim)

### get peak list
d10A_matched_stable_pk_id = rep(NA, 1000000)
d20A_matched_stable_pk_id = rep(NA, 1000000)
d10A_matched_stable_pk_allmatched_id = rep(NA, 1000000)
d20A_nomatch_pk_id = rep(NA, 1000000)

a1 = 0 ### d20A unmatched
a2 = 0 ### d10A only uniq
a3 = 0 ### d10A all 
set.seed(2019)
for (i in 1:length(d20A)){
if (i%%1000==0){
	print(i)
}
d20A_tmp = d20A[i]
dif_value = 1

used_id_binary = (d10A>(d20A_tmp-dif_value)) & (d10A<(d20A_tmp+dif_value))
if (sum(used_id_binary)==0){
a1 = a1+1
d20A_nomatch_pk_id[a1] = i
} else if (sum(used_id_binary)==1){
a2 = a2+1
a3 = a3+1
d10A_matched_stable_pk_id[a2] = which(used_id_binary)
d20A_matched_stable_pk_id[a2] = i
d10A_matched_stable_pk_allmatched_id[a3] = which(used_id_binary)
} else {
used_id_pos_od = which(used_id_binary)
used_id_pos <- used_id_pos_od[!used_id_pos_od %in% d10A_matched_stable_pk_id]
if (length(used_id_pos)==0){
used_id_pos=used_id_pos_od
}
used_id_pos_s1 = sample(length(used_id_pos), 1)
used_id_pos1 = used_id_pos[used_id_pos_s1]
### add id
a2 = a2+1
d10A_matched_stable_pk_id[a2] = used_id_pos1
d20A_matched_stable_pk_id[a2] = i
a3a = a3+1
a3b = a3+length(used_id_pos_od)
d10A_matched_stable_pk_allmatched_id[a3a:a3b] = used_id_pos_od
a3 = a3b
}
}
### clean ids
d10A_matched_stable_pk_allmatched_id_uniq = unique(d10A_matched_stable_pk_allmatched_id)
d10A_matched_stable_pk_allmatched_id_uniq_narm = d10A_matched_stable_pk_allmatched_id_uniq[!is.na(d10A_matched_stable_pk_allmatched_id_uniq)]
d10A_matched_stable_pk_id_narm = d10A_matched_stable_pk_id[!is.na(d10A_matched_stable_pk_id)]
d20A_matched_stable_pk_id_narm = d20A_matched_stable_pk_id[!is.na(d20A_matched_stable_pk_id)]
d20A_nomatch_pk_id_narm = d20A_nomatch_pk_id[!is.na(d20A_nomatch_pk_id)]
### get id to mat
d10A_matched_stable_pk_allmatched = d1[d10A_matched_stable_pk_allmatched_id_uniq_narm,]
d10A_matched_stable_pk = d1[d10A_matched_stable_pk_id_narm,]
d20A_matched_stable_pk = d2[d20A_matched_stable_pk_id_narm,]
d20A_nomatch_pk = d2[d20A_nomatch_pk_id_narm,]
###### write all split results
write.table(d10A_matched_stable_pk, paste(new_folder, output_name, '.unstable.submatch.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
write.table(d10A_matched_stable_pk_allmatched, paste(new_folder, output_name, '.unstable.allmatch.uniq.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
write.table(d20A_matched_stable_pk, paste(new_folder, output_name, '.stable.submatch.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
write.table(d20A_nomatch_pk, paste(new_folder, output_name, '.stable.unmatch.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')

d10A_matched_stable_pk_0Ameansig = d10A_matched_stable_pk[,5]/2+d10A_matched_stable_pk[,6]/2
d20A_matched_stable_pk_0Ameansig = d20A_matched_stable_pk[,5]/2+d20A_matched_stable_pk[,6]/2
write.table(d10A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig>Amean_lim,], paste(new_folder, output_name, '.unstable.submatch.highsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
write.table(d20A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig>Amean_lim,], paste(new_folder, output_name, '.stable.submatch.highsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
#write.table(d10A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig<=Amean_lim,], paste(new_folder, output_name, '.unstable.submatch.lowsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
write.table(d20A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig<=Amean_lim,], paste(new_folder, output_name, '.stable.submatch.lowsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')

d1_signal = d1[d10A<=Amean_lim,]
set.seed(2019)
used_id = sample(dim(d1_signal)[1], sum(d20A_matched_stable_pk_0Ameansig>Amean_lim))
write.table(d1_signal[used_id,], paste(new_folder, output_name, '.unstable.submatch.lowsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')



d1_signal = d1[d10A<=Amean_lim,]
write.table(d1_signal, paste(new_folder, output_name, '.unstable.all.lowsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
d1_signal = d1[d10A>Amean_lim,]
write.table(d1_signal, paste(new_folder, output_name, '.unstable.all.highsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')

d2_signal = d2[d20A<=Amean_lim,]
write.table(d2_signal, paste(new_folder, output_name, '.stable.all.lowsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')
d2_signal = d2[d20A>Amean_lim,]
write.table(d2_signal, paste(new_folder, output_name, '.stable.all.highsig.txt', sep=''), quote=F, col.names=T, row.names=F, sep='\t')




#balanced_peaklist_for_compare.stable.submatch.txt
#balanced_peaklist_for_compare.unstable.submatch.txt
#balanced_peaklist_for_compare.stable.allmatch.uniq.txt
#balanced_peaklist_for_compare.unstable.unmatch.txt
#balanced_peaklist_for_compare.stable.highsig.txt
#balanced_peaklist_for_compare.unstable.highsig.txt


### for plotting
tp_plot = c(-0.25,0.25, 3.75,4.25, 5.75,6.25, 11.75,12.25, 17.75,18.25, 23.75,24.25)
tp_plot_forcurve = c(0,4,6,12,18,24)
line_num_plot = 2000
#
pdf(paste(new_folder, 'stable.pk.pdf', sep=''))
sigmat = as.matrix(d20A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig>Amean_lim,][,-c(1:4)])+1
print(dim(sigmat))
#sigmat = sigmat[!duplicated(sigmat),]
sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))
boxplot(sigmat, ylim=c(1, 60) ,outline=FALSE, at =tp_plot, log='y')
if (dim(sigmat)[1]>line_num_plot){line_num = sample(dim(sigmat)[1], line_num_plot)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
	lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(tp_plot_forcurve, colMedian((sigmat_mean_for_curve)), col='black')
dev.off()
#
pdf(paste(new_folder, 'unstable.pk.pdf', sep=''))
sigmat = as.matrix(d10A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig>Amean_lim,][,-c(1:4)])+1
print(dim(sigmat))
#sigmat = sigmat[!duplicated(sigmat),]
sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))
boxplot(sigmat, ylim=c(1, 60) ,outline=FALSE, at =tp_plot, log='y')
if (dim(sigmat)[1]>line_num_plot){line_num = sample(dim(sigmat)[1], line_num_plot)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
	lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(tp_plot_forcurve, colMedian((sigmat_mean_for_curve)), col='black')
dev.off()
#
pdf(paste(new_folder, 'stable.lowsig.pk.pdf', sep=''))
sigmat = as.matrix(d20A_matched_stable_pk[d20A_matched_stable_pk_0Ameansig<=Amean_lim,][,-c(1:4)])+1
print(dim(sigmat))
#sigmat = sigmat[!duplicated(sigmat),]
sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))
boxplot(sigmat, ylim=c(1, 60) ,outline=FALSE, at =tp_plot, log='y')
if (dim(sigmat)[1]>line_num_plot){line_num = sample(dim(sigmat)[1], line_num_plot)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
        lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(tp_plot_forcurve, colMedian((sigmat_mean_for_curve)), col='black')
dev.off()
#
pdf(paste(new_folder, 'unstable.lowsig.pk.pdf', sep=''))
d1_signal = d1[d10A<=Amean_lim,]
set.seed(2019)
used_id = sample(dim(d1_signal)[1], sum(d20A_matched_stable_pk_0Ameansig>Amean_lim))
sigmat = as.matrix(d1_signal[used_id,][,-c(1:4)])+1
print(dim(sigmat))
#sigmat = sigmat[!duplicated(sigmat),]
sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))
boxplot(sigmat, ylim=c(1, 60) ,outline=FALSE, at =tp_plot, log='y')
if (dim(sigmat)[1]>line_num_plot){line_num = sample(dim(sigmat)[1], line_num_plot)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
        lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(tp_plot_forcurve, colMedian((sigmat_mean_for_curve)), col='black')
dev.off()




