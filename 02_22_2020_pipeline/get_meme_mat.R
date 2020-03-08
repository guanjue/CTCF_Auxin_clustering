args = commandArgs(trailingOnly=TRUE)

input = args[1]
output_file = paste(input, '.meme.mat.txt', sep='')

d1 = read.table(input, header=F)
d1_c = d1[,-1]

toMEME_col = function(x){
	meme_col = c(0,0,0,0)
	meme_col[1] = sum(x=='A')
	meme_col[2] = sum(x=='C')
	meme_col[3] = sum(x=='G')
	meme_col[4] = sum(x=='T')
	return(meme_col)
}

###
meme_mat = apply(d1_c, 2, toMEME_col)
write.table(meme_mat, output_file, quote=F, sep='\t', col.names=F, row.names=F)

for i in {3..1}
do
time Rscript ~/scratch/ctcf_auxin/get_meme_mat.R 'ctcf2_c'$i'_fimo.bestmatch.5exp.fa.categoricalvec.txt'
done
