colMedians = function(x){
	return(apply(x, 2, function(x) median(x, na.rm=T)))
}

robust_mean = function(x){
	x = x[!is.na(x)] 
	robust_mean_tmp = mean(x[(x<quantile(x,0.95))&(x>quantile(x,0.05))])
	return(robust_mean_tmp)
}

colMeans_robust = function(x){
	return(apply(x, 2, function(x) robust_mean(x)) )
}

color_list = c('black', 'gray','purple','blue','green', 'orange', 'white', 'red')
pdf(paste('all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.pdf',sep=''))
file = 'all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.unstable.lowsig.subs.bed.txt'
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean[1:length(dIS_mean)/2])
plot(1:length(dIS_mean), dIS_mean, col='black', ylim=c(-0.1,0.3), type='l')
### l2
file = 'all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.unstable.highsig.subs.bed.txt'
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean[1:length(dIS_mean)/2])
lines(1:length(dIS_mean), dIS_mean, col='blue')
### l3
file = 'all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.stable.highsig.uniq.bed.txt'
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean[1:length(dIS_mean)/2])
lines(1:length(dIS_mean), dIS_mean, col='red')
dev.off()


color_list = c('black', 'gray','purple','blue','green', 'orange', 'white', 'red')
pdf(paste('all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.relative.pdf',sep=''))
file = 'all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.unstable.lowsig.subs.bed.txt'
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean[1:length(dIS_mean)/2])
plot(1:length(dIS_mean), dIS_mean, col='black', ylim=c(-0.1,0.3), type='l')
### l2
file = 'all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.unstable.highsig.subs.bed.txt'
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean[1:length(dIS_mean)/2])
lines(1:length(dIS_mean), dIS_mean, col='blue')
### l3
file = 'all_pk/LMqPCRnorm_decay_folder/balanced_peaklist_for_compare.stable.highsig.uniq.bed.txt'
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean[1:length(dIS_mean)/2])
lines(1:length(dIS_mean), dIS_mean, col='red')
dev.off()































for (i in 1:9){
	print(i)
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
print(i)
print(dim(dIS))
}

for (i in 1:9){
	print(i)
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
png(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt.png',sep=''))
boxplot(dIS, ylim=c(-1,1))
dev.off()
}

for (i in c(1:9)){
	print(i)
if ((i!=60)&(i!=80)){
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
png(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt.png',sep=''))
boxplot(dIS, ylim=c(-1,1))
dev.off()
} else{
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i+1,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = rbind(dIS, d[,-c(1:6)])
png(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i, '_', i+1,'.pk.clean.bed.IS.txt.png',sep=''))
boxplot(dIS, ylim=c(-1,1))
dev.off()
}
}


colMedians = function(x){
	return(apply(x, 2, function(x) median(x, na.rm=T)))
}

robust_mean = function(x){
	x = x[!is.na(x)] 
	robust_mean_tmp = mean(x[(x<quantile(x,0.95))&(x>quantile(x,0.05))])
	return(robust_mean_tmp)
}

colMeans_robust = function(x){
	return(apply(x, 2, function(x) robust_mean(x)) )
}

color_list = c('black', 'gray','purple','blue','green','yellow', 'orange', 'pink', 'red')
pdf(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.all.pk.clean.bed.IS.txt.pdf',sep=''))
i=1
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean)
plot(1:length(dIS_mean), dIS_mean, col=paste('gray', i*9-9, sep=''), ylim=c(-0.1,0.5), type='l')
for (i in 1:9){
print(i)
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean)
print(head(dIS_mean))
lines(1:length(dIS_mean), dIS_mean, col=color_list[i])
}
dev.off()

color_list = c('black', 'gray','purple','blue','green','yellow', 'orange', 'pink', 'red')
pdf(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.all.pk.clean.bed.IS.txt.dif.pdf',sep=''))
i=1
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean)#[1:length(dIS_mean)/2])
plot(1:length(dIS_mean), dIS_mean, col=paste('gray', i*9-9, sep=''), ylim=c(-0.5,0.1), type='l')
for (i in 1:9){
print(i)
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean)#[1:length(dIS_mean)/2])
print(head(dIS_mean))
lines(1:length(dIS_mean), dIS_mean, col=color_list[i])
}
dev.off()




color_list = c('black', 'gray','purple','blue','green', 'orange', 'white', 'red')
pdf(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.all.pk.clean.bed.IS.txt.67_89merged.pdf',sep=''))
i=1
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean[1:length(dIS_mean)/2])
plot(1:length(dIS_mean), dIS_mean, col=paste('gray', i*9-9, sep=''), ylim=c(-0.1,0.3), type='l')
for (i in c(1:5,6,8)){
print(i)
if ((i!=6)&(i!=8)){
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
} else{
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i+1,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = rbind(dIS, d[,-c(1:6)])	
}
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean#-max(dIS_mean)#[1:length(dIS_mean)/2])
print(head(dIS_mean))
lines(1:length(dIS_mean), dIS_mean, col=color_list[i])
}
dev.off()



color_list = c('black', 'gray','purple','blue','green', 'orange', 'white', 'red')
pdf(paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.all.pk.clean.bed.IS.txt.dif.67_89merged.pdf',sep=''))
i=1
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean[1:length(dIS_mean)/2])
plot(1:length(dIS_mean), dIS_mean, col=paste('gray', i*9-9, sep=''), ylim=c(-0.3,0.1), type='l')
for (i in c(1:5,6,8)){
print(i)
if ((i!=6)&(i!=8)){
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
} else{
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = d[,-c(1:6)]
file = paste('all_pk/LMqPCRnorm_decay_folder/raw_kmeans.',i+1,'.pk.clean.bed.IS.txt',sep='')
d = read.table(file, header=F,  skip = 1)
dIS = rbind(dIS, d[,-c(1:6)])	
}
dIS_mean = colMeans_robust(dIS)#,na.rm=T)
dIS_mean = dIS_mean-max(dIS_mean)#[1:length(dIS_mean)/2])
print(head(dIS_mean))
lines(1:length(dIS_mean), dIS_mean, col=color_list[i])
}
dev.off()




