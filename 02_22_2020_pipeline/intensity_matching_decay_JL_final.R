library(mixtools)
library(dplyr)
library(ggplot2)
library(data.table)


### function
colMedian = function(x){
  return(apply(x, 2, median))
}

d1 = read.table('all.unstable.txt', header=T) #rapidly.degraded; 34747
d2 = read.table('all.stable.txt', header=T) #persistent; 4097

output_name = 'matched'

d10A = d1[,5]/2+d1[,6]/2
d20A = d2[,5]/2+d2[,6]/2



###### remove relatively low signals
### rapidly.degraded 0A signal
all_0A_signal = c(d10A)

### split residuals
set.seed(2019)
mixmdl = normalmixEM(log(all_0A_signal+1), k=2)
pdf(paste('extract_peaks_for_compare.signal.hist.0A.b.pdf', sep=''))
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


## Separte rapidly degraded peaks into lowsig (cluster I) and highsig (cluster II) subsets
d1_all_lowsig = d1[d10A<=Amean_lim,]
write.table(d1_all_lowsig, 'Unstable.all.lowsig.txt', quote=F, col.names=T, row.names=F, sep='\t')

d1_all_highsig = d1[d10A>Amean_lim,]
write.table(d1_all_highsig, 'Unstable.all.highsig.txt', quote=F, col.names=T, row.names=F, sep='\t') #14694

pdf('Intensity_clusterI.pdf')
sigmat = as.matrix(d1_all_lowsig[,-c(1:4)])+1
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
pdf('Intensity_clusterII.pdf')
sigmat = as.matrix(d1_all_highsig[,-c(1:4)])+1
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
pdf('Intensity_ClusterIII.pdf')
sigmat = as.matrix(d2[,-c(1:4)])+1
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



### Start intensity matching ###
## Read in files
d1_all_highsig<-read.table('Unstable.all.highsig.txt', header=T) #14694; cluster II
d2<-read.table("all.stable.txt",header=T) #4097; cluster III

# Intensity matching betweeen Cluster II ("unstable") and Cluster III ("persistent")

var=0.01

for (n in 1:3){
  print(n)
  persistent.sampled=NULL
  unstable.sampled=NULL
  
  persistent=d2
  unstable=d1_all_highsig[sample(nrow(d1_all_highsig)),]

  for (i in 1:nrow(persistent)){
    a=unstable%>%filter(((unstable[,5]+unstable[,6])/2>=(1-var)*(persistent[i,5]+persistent[i,6])/2) & ((unstable[,5]+unstable[,6])/2<=(1+var)*(persistent[i,5]+persistent[i,6])/2))
    if (nrow(a)>0){
      persistent.sampled <- rbindlist(list(persistent.sampled,persistent[i,]))
      unstable.sampled <- rbindlist(list(unstable.sampled,a[1,]))
      unstable=anti_join(unstable,a[1,],by=c("id")) %>% as.data.table()
    }
  }
  write.table(persistent.sampled,paste0("persistent.sampled.run",n,".txt",sep=""),sep="\t",eol="\n",quote=F,row.names=F,col.names=T)
  write.table(unstable.sampled,paste0("unstable_highsig.sampled.run",n,".txt",sep=""),sep="\t",eol="\n",quote=F,row.names=F,col.names=T)
  print (nrow(persistent.sampled))
  print (nrow(unstable.sampled))
}



## plot binding density for sampled lists
currentdir="/Volumes/Jing_Luan/Lab_work/TAD_CTCF/Computation/CTCF_Dynamics_PSU/2020_Guanjue/Decay/Sub-stratification/intensity_matching"
files<-list.files(path=currentdir,pattern="run[[:digit:]].txt", full.names=TRUE)
files

# a few plot settings
tp_plot = c(-0.25,0.25, 3.75,4.25, 5.75,6.25, 11.75,12.25, 17.75,18.25, 23.75,24.25)
tp_plot_forcurve = c(0,4,6,12,18,24)
line_num_plot = 2000

for (i in 1:length(files)) {
  # 
  file<-read.table(files[i],header=T,sep="\t")
  #
  pdf(paste0(gsub("txt","",files[i]), 'sig.pdf', sep=''))
  sigmat = as.matrix(file[,-c(1:4)])+1
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
}




