library(mixtools)
library(dplyr)
library(data.table)


colMedian = function(x){
  return(apply(x, 2, median))
}

mclust6 = read.table('raw_mclust.6.pk.clean.txt', header=T) #unstable; 6804
mclust5 = read.table('raw_mclust.5.pk.clean.txt', header=T) #stable; 5125


### Start intensity matching ###

var=0.01

for (n in 1:3){
  print(n)
  mclust6.sampled=NULL
  mclust5.sampled=NULL
  
  persistent=mclust6
  unstable=mclust5[sample(nrow(mclust5)),]
  
  for (i in 1:nrow(persistent)){
    a=unstable%>%filter(((unstable[,5]+unstable[,6])/2>=(1-var)*(persistent[i,5]+persistent[i,6])/2) & ((unstable[,5]+unstable[,6])/2<=(1+var)*(persistent[i,5]+persistent[i,6])/2))
    if (nrow(a)>0){
      mclust6.sampled <- rbindlist(list(mclust6.sampled,persistent[i,]))
      mclust5.sampled <- rbindlist(list(mclust5.sampled,a[1,]))
      unstable=anti_join(unstable,a[1,],by=c("id")) %>% as.data.table()
    }
  }
  write.table(mclust6.sampled,paste0("raw_mclust.6.sampled.run",n,".txt",sep=""),sep="\t",eol="\n",quote=F,row.names=F,col.names=T)
  write.table(mclust5.sampled,paste0("raw_mclust.5.sampled.run",n,".txt",sep=""),sep="\t",eol="\n",quote=F,row.names=F,col.names=T)
  print (nrow(mclust6.sampled))
  print (nrow(mclust5.sampled))
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




