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

small_num = 1e-5
###### remove relatively low signals
### unstable 0A signal
all_0A_signal = c(d10A)
### split residuals
set.seed(2019)
mixmdl = normalmixEM(log(all_0A_signal+small_num), k=2)
pdf(paste(new_folder, 'extract_peaks_for_compare.signal.hist.0A.b.pdf', sep=''))
plot(mixmdl,which=2)
lines(density(log(all_0A_signal+small_num)), lty=2, lwd=2)
dev.off()

### get Amean lim
gmm_posterior = mixmdl$posterior
print('mixmdl$mu')
print(mixmdl$mu)
used_id = (gmm_posterior[,2]>0.5)
Amean_lim = min((all_0A_signal)[used_id])  #min(log(all_0A_signal+1)[used_id])
print(Amean_lim)


### high signal
d2_high = d2[d20A>Amean_lim,]


### for plotting
tp_plot = c(-0.25,0.25, 3.75,4.25, 5.75,6.25, 11.75,12.25, 17.75,18.25, 23.75,24.25)
tp_plot_forcurve = c(0,4,6,12,18,24)
tp_plot_forcurve_all = c(0,0,4,4,6,6,12,12,18,18,24,24)
line_num_plot = 3000

### get sig
sigmat = as.matrix(d2_high[,-c(1:4)])
print(dim(sigmat))
sigmat_colMedian = colMedian(sigmat)
#sigmat = sigmat[!duplicated(sigmat),]
sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))


### get curve
y_x = c()
for (i in 1:dim(sigmat)[2]){
y_tmp = as.data.frame(sigmat[,i])
x_tmp = as.data.frame(rep(tp_plot_forcurve_all[i], dim(sigmat)[1]))
y_x_tmp = cbind(y_tmp,x_tmp)
colnames(y_x_tmp) = c('y','x')
y_x = rbind(y_x, y_x_tmp)
}

### get curve
get_decay_curve = function(y,x){
A = mean(y[x==0])
y_sub = y - A
a = lm(y_sub~x-1)
#print(summary(a))
#print(summary(a$coefficients))
return(c(A, a$coefficients[1], summary(a)$r.squared))
}

### fit model
expdecay_curve_average = get_decay_curve(log(y_x[,1]+small_num), y_x[,2])
expdecay_curve_median = get_decay_curve(log(sigmat_colMedian+1), tp_plot_forcurve_all)

expdecay_curve_mat = t(apply(sigmat, 1, function(x) get_decay_curve(log(x+1), tp_plot_forcurve_all)))
expdecay_curve_mat_median = expdecay_curve_mat[expdecay_curve_mat[,1]==median(expdecay_curve_mat[,1]),]

curve_x = seq(0,24)
curve_y = exp(expdecay_curve_average[1] + curve_x * expdecay_curve_average[2])-1

# plotting
pdf(paste(new_folder, 'allhighsig.stable.pk.expdecay.pdf', sep=''), width=5, height=5)
boxplot(sigmat, ylim=c(1, 60) ,outline=FALSE, at =tp_plot, log='')
if (dim(sigmat)[1]>line_num_plot){line_num = sample(dim(sigmat)[1], line_num_plot)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
	lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
#lines(tp_plot_forcurve, colMedian((sigmat_mean_for_curve[!is.na(sigmat_mean_for_curve[,1]),])), col='black', lwd=2)
lines(curve_x, curve_y, col='dodgerblue1', lwd=2)
dev.off()

A = 3.3784207
B = -0.0909538

