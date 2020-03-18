args = commandArgs(trailingOnly=TRUE)

input_matrix = args[1]
output_name = args[2]
ydw = as.numeric(args[3])
yup = as.numeric(args[4])

#input_matrix = 'input_matrix.txt'
#output_name = 'input_matrix.pdf'
#yup = 60
#ydw = 1

sigmat0 = read.table(input_matrix, header=T)
sigmat = as.matrix(sigmat0[,-c(1:4)])+1

### for plotting
set.seed(2019)
tp_plot = c(-0.25,0.25, 3.75,4.25, 5.75,6.25, 11.75,12.25, 17.75,18.25, 23.75,24.25)
tp_plot_forcurve = c(0,4,6,12,18,24)
line_num_plot = 2000
#
pdf(output_name)
print(dim(sigmat))
sigmat_mean_for_curve = as.matrix(cbind(sigmat[,1]/2+sigmat[,2]/2, sigmat[,3]/2+sigmat[,4]/2, sigmat[,5]/2+sigmat[,6]/2, sigmat[,7]/2+sigmat[,8]/2, sigmat[,9]/2+sigmat[,10]/2, sigmat[,11]/2+sigmat[,12]/2))
boxplot(sigmat, ylim=c(ydw, yup) ,outline=FALSE, at =tp_plot, log='y')
if (dim(sigmat)[1]>line_num_plot){line_num = sample(dim(sigmat)[1], line_num_plot)} else {line_num=1:dim(sigmat)[1]}
for (j in (line_num)){
	lines(tp_plot_forcurve, (sigmat_mean_for_curve[j,]), col=rgb(255/255,0/255,0/255,alpha=0.05) )
}
#
lines(tp_plot_forcurve, colMedian((sigmat_mean_for_curve)), col='black')
dev.off()
