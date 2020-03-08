### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_mat = args[1]
signal_mat_1kb = args[2]
signal_mat_5kb = args[3]
signal_mat_10kb = args[4]
output_mat = args[5]

library(data.table)

signal_mat = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt'
signal_mat_1kb = 'ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt'
signal_mat_5kb = 'ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt'
signal_mat_10kb = 'ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt'
output_mat = 'ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt'

### get input sig
sig = as.data.frame(fread(signal_mat))
sig_1kb_bg = as.data.frame(fread(signal_mat_1kb))
sig_5kb_bg = as.data.frame(fread(signal_mat_5kb))
sig_10kb_bg = as.data.frame(fread(signal_mat_10kb))

### get local bg sig
sig_mat = sig[,-c(1:4)]
sig_1kb = sig_1kb_bg[,-c(1:4)]
sig_5kb = sig_5kb_bg[,-c(1:4)]
sig_10kb = sig_10kb_bg[,-c(1:4)]

### normalize to bg
bg_used_id = (grepl("randbg", sig[,4]))
# get sf
sf_1kb = colMeans(sig_mat[bg_used_id,]) / colMeans(sig_1kb[bg_used_id,])
sf_5kb = colMeans(sig_mat[bg_used_id,]) / colMeans(sig_5kb[bg_used_id,])
sf_10kb = colMeans(sig_mat[bg_used_id,]) / colMeans(sig_10kb[bg_used_id,])
# bgnorm
sf_1kb_norm = t(apply(sig_1kb, 1, function(x) x*sf_1kb))
sf_5kb_norm = t(apply(sig_5kb, 1, function(x) x*sf_5kb))
sf_10kb_norm = t(apply(sig_10kb, 1, function(x) x*sf_10kb))

### get wg bg sig
sig_mat_wg_bgsig = colMeans( (sf_1kb_norm[bg_used_id,]+sf_5kb_norm[bg_used_id,]+sf_10kb_norm[bg_used_id,])/3 )

### get MACS bg sig
sig_MACS_bgsig = sf_1kb_norm
for (i in 1:dim(sf_1kb_norm)[2]){
	wh_mean_tmp = sig_mat_wg_bgsig[i]
	sig_MACS_bgsig[sig_MACS_bgsig[,i]<wh_mean_tmp,i] = wh_mean_tmp
	bgsig_5kb = sf_5kb_norm[,i]
	sig_MACS_bgsig[sig_MACS_bgsig[,i]<bgsig_5kb,i] = bgsig_5kb[sig_MACS_bgsig[,i]<bgsig_5kb]
	bgsig_10kb = sf_10kb_norm[,i]
	sig_MACS_bgsig[sig_MACS_bgsig[,i]<bgsig_10kb,i] = bgsig_10kb[sig_MACS_bgsig[,i]<bgsig_10kb]
}
sig_MACS_bgsig[sig_MACS_bgsig<0] = 0

### get sig - MACS_bgsig
sig_mat_bg_sub = sig_mat - sig_MACS_bgsig
sig_mat_bg_sub[sig_mat_bg_sub<0] = 0
sig_mat_bg_sub_mat = cbind(sig[,c(1:4)], sig_mat_bg_sub)
write.table(sig_mat_bg_sub_mat, output_mat, quote=F, sep='\t', col.names=F, row.names=F)


