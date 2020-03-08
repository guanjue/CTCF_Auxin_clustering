rm ctcf2_c*_fimo.bestmatch.5exp.bed.matrix.plastcons.txt
for i in {3..1}
do
        echo $i
        time computeMatrix reference-point --referencePoint center --sortRegions keep -S ~/group/genome/mm9/plastCons/plastcons30way.euarchontoglires.bw -R 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed' -a 30 -b 30 --binSize 1 --missingDataAsZero --numberOfProcessors 4 --outFileName 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed.matrix.plastcons.txt.gz'
        time gunzip 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed.matrix.plastcons.txt.gz'
        tail -n+2 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed.matrix.plastcons.txt' > 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed.matrix.plastcons.txt.tmp' && mv 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed.matrix.plastcons.txt.tmp' 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed.matrix.plastcons.txt'
done

time Rscript plot_plastcon.R
