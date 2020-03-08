setwd('/storage/home/gzx103/scratch/ctcf_auxin/all_pk/LMqPCRnorm_decay_folder')

library(ggplot2)
library(RColorBrewer)
library(ggpubr)

col_sd = function(x){
	return(sd(x))
}

d_mean = c()
d_sd = c()
d_mean_ci_up = c()
d_mean_ci_down = c()
z = 1.96 # 95% CI
d_df = as.data.frame(c())
d_mean_sd_df = as.data.frame(c())
colors = c("#999999", "#FF8C00", "#FF4500")
cluster_name = c('UN-stable-Low', 'UN-stable-High', 'Stable-High')

for (i in c(1:3)){
print(i)
d1 = read.table(paste('ctcf2_c', i, '_fimo.bestmatch.5exp.bed.matrix.plastcons.txt', sep=''), header=F)[,-c(1:6)]
d_mean_i = as.data.frame(colMeans(d1))
d_sd_i = apply(d1, 2, col_sd)
d_ci_i = z * d_sd_i / sqrt(dim(d1)[1])
d_mean_sd_df = rbind(d_mean_sd_df, cbind(d_mean_i, d_ci_i, c(1:length(d_sd_i))-length(d_sd_i)/2, rep(cluster_name[i], length(d_sd_i)), rep(colors[4-i], length(d_sd_i)) ))
for (j in c(1:dim(d1)[2])){
info_tmp = cbind(d1[,i], rep(j,dim(d1)[1]), rep(i,dim(d1)[1]))
d_df = rbind(d_df, info_tmp)
}
}

colnames(d_df) = c('PC', 'P', 'C')
colnames(d_mean_sd_df) = c('PC', 'ci', 'P', 'C', 'colors')

pdf('phastCon.pdf', width=7, height=4)#, width=500, height=500)
p = ggplot(d_mean_sd_df, aes(x=P, y=PC, group=C, color=C)) 
p = p + geom_line()
p = p + geom_point()
p = p + geom_errorbar(aes(ymin=PC-ci, ymax=PC+ci), width=2,
                 position=position_dodge(0.05))
p = p + scale_color_manual(values=c("#999999", "#FF8C00", "#FF4500"))
plot(p)
dev.off()


