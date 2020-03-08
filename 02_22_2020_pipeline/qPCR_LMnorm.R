### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_mat_file = args[1]
signal_mat_tp_file = args[2]
qPCR_mat_file = args[3]
qPCR_mat_tp_file = args[4]
output_file_name = args[5]


#signal_mat_file = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt'
#signal_mat_tp_file = 'bam_bw_bedgraph/file_list.all.timepoints.txt'
#qPCR_mat_file = '20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.sig.txt'
#qPCR_mat_tp_file = '20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.timepoints.txt'
#output_file_name = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.txt'
library(data.table)

### get input mat
chipseq_mat = as.data.frame(fread(signal_mat_file))
qPCR_mat = as.data.frame(fread(qPCR_mat_file))
rm_rows = !((qPCR_mat[,4]=="qPCR_chr8_92095917_92095999") | (qPCR_mat[,4]=="qPCR_chrX_98417947_98418017"))
qPCR_mat_modifed = qPCR_mat[rm_rows,]
qPCR_mat_modifed = qPCR_mat_modifed[order(qPCR_mat_modifed[,4]),]
qPCR_mat_onlysig = as.matrix(qPCR_mat_modifed[,-c(1:4)])

### add tp
chipseq_sig_tp = read.table(signal_mat_tp_file, header=F)
chipseq_sig_tp = chipseq_sig_tp[,3]
qPCR_sig_tp = read.table(qPCR_mat_tp_file, header=F)

### get chip sig at qPCR regions
qPCR_region_id = (grepl("qPCR", chipseq_mat[,4]))
chipseq_sig_qPCR_regions = chipseq_mat[qPCR_region_id,]
rm_rows = !((chipseq_sig_qPCR_regions[,4]=="qPCR_chr8_92095917_92095999") | (chipseq_sig_qPCR_regions[,4]=="qPCR_chrX_98417947_98418017"))
chipseq_sig_qPCR_regions_modified = chipseq_sig_qPCR_regions[rm_rows,]
chipseq_sig_qPCR_regions_modified = chipseq_sig_qPCR_regions_modified[order(chipseq_sig_qPCR_regions_modified[,4]),]
chipseq_sig_qPCR_regions_onlysig = as.matrix(chipseq_sig_qPCR_regions_modified[,-c(1:4)])

### change qPCR signal to ChIP-seq level
qPCR_mat_onlysig_log2 = log2(qPCR_mat_onlysig[qPCR_mat_onlysig!=0])
chipseq_sig_qPCR_regions_onlysig_log2 = log2(chipseq_sig_qPCR_regions_onlysig[chipseq_sig_qPCR_regions_onlysig!=0])
B = sd(chipseq_sig_qPCR_regions_onlysig_log2) / sd(qPCR_mat_onlysig_log2)
A = 2^(mean(chipseq_sig_qPCR_regions_onlysig_log2) - mean(qPCR_mat_onlysig_log2) * B)
qPCR_mat_onlysig_rescale = qPCR_mat_onlysig ^ B * A

### get LMnorm
chipseq_mat_only_sig_rescale = chipseq_mat[,-c(1:4)]
for (i in 1:length(chipseq_sig_tp)){
	print(i)
	chip_tp_tmp = chipseq_sig_tp[i]
	chip_tp_sig_tmp = chipseq_sig_qPCR_regions_onlysig[,i]
	qPCR_used_col = qPCR_sig_tp==toString(chip_tp_tmp)
	qPCR_tp_sig_tmp = qPCR_mat_onlysig_rescale[,qPCR_used_col]
	chip_tp_sig_tmp_log2 = log2(chip_tp_sig_tmp+0.1)
	qPCR_tp_sig_tmp_log2 = log2(qPCR_tp_sig_tmp+0.1)
	### get A & B
	LM_model_tmp = lm(qPCR_tp_sig_tmp~chip_tp_sig_tmp-1)
	chip_tp_sig_tmp_A = LM_model_tmp$coefficients[1]
	### rescale chip
	print(LM_model_tmp$coefficients)
	#print(chip_tp_sig_tmp_A)
	chip_tp_sig_all_tmp = chipseq_mat_only_sig_rescale[,i]
	chip_tp_sig_tmp_rescale = chip_tp_sig_all_tmp * chip_tp_sig_tmp_A
	chipseq_mat_only_sig_rescale[,i] = chip_tp_sig_tmp_rescale
}

chipseq_mat_only_sig_rescale = as.matrix(chipseq_mat_only_sig_rescale)
chipseq_mat_only_sig_rescale[!is.finite(chipseq_mat_only_sig_rescale)] = 0
colnames(chipseq_mat_only_sig_rescale) = chipseq_sig_tp


####################################################
### use pheatmap to plot heatmaps
color_heatmap = function(color_matrix, high_color, low_color, format, outputname){
	library(pheatmap)
	format(outputname, width = 500, height = 1000) ### output name
	par() ### set heatmap margins
	### plot pheatmap
	my_colorbar=colorRampPalette(c(low_color, high_color))(n = 128)
	col_breaks = c(seq(0, 2000,length=33))
	pheatmap(color_matrix, color=my_colorbar, cluster_cols = FALSE,cluster_rows=TRUE, clustering_method = 'average',annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE)
	dev.off()
}


chipseq_mat_only_sig_rescale_plot = chipseq_mat_only_sig_rescale[c(1:(38844+2000)),]

set.seed(2019)
used_id = sample(38844, 5000)
chipseq_mat_only_sig_rescale_plot = chipseq_mat_only_sig_rescale[used_id,]
chipseq_mat_only_sig_rescale_plot = rbind(chipseq_mat_only_sig_rescale_plot, chipseq_mat_only_sig_rescale[38844:40000,])

format = png
color_heatmap(chipseq_mat_only_sig_rescale_plot, 'red', 'white', format, 'rescaled_mat.png')

output_matrix = cbind(chipseq_mat[,1:4], chipseq_mat_only_sig_rescale)
colnames(output_matrix)[1:4] = c('chr','start','end','id')
write.table(output_matrix, output_file_name, quote=F, col.names=T, row.names=F, sep='\t')

