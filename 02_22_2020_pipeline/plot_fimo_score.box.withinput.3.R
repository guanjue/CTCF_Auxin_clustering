library(ggplot2)
library(RColorBrewer)
library(ggpubr)

### get parameters
args = commandArgs(trailingOnly=TRUE)
header_name = args[1]
output = args[2]

time Rscript ~/scratch/ctcf_auxin/Mclust_pk/plot_fimo_score.box.R 'Decay_a/Matched/match_c' ED_match.fimo.pdf
time Rscript ~/scratch/ctcf_auxin/Mclust_pk/plot_fimo_score.box.R 'Decay_a/All/all_c' ED_all.fimo.pdf

d_all = c()
for (i in c(1:3)){
d1 = read.table(paste(header_name, i, '_fimo.bestmatch.txt', sep=''), header=F)
d_all = rbind(d_all, cbind(as.data.frame(-log10(d1[,8])), as.factor(rep(i, dim(d1)[1]))))
}


d_all = as.data.frame(d_all)
colnames(d_all) = c('FIMO_log10_pval', 'cluster')

pdf(output, width=4, height=4)#, width=500, height=500)
p = ggplot(data = d_all, aes(x=cluster, y=FIMO_log10_pval)) 
p = p + geom_boxplot(aes(fill = cluster))
#p = p + geom_point(aes(y=FIMO_log10_pval, group=cluster), position = position_dodge(width=0.75))
p = p + scale_fill_manual(values=c("#999999", "#FF8C00", "#FF4500"))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size = 12))
my_comparisons = list( c("2", "3"), c("1", "3") )
p = p + stat_compare_means(aes(group = cluster), comparisons = my_comparisons, label = "p.format", paired = F, method = "t.test", method.args=list(alternative = "less")) 
p = p + coord_cartesian(ylim = c(0, 13))
plot(p)
dev.off()


