cd /storage/home/gzx103/scratch/ctcf_auxin/Mclust_pk/Mclust_a/All

### bed file list of each cluster
ls -ltrh bed_list.txt

### get CTCF motif centered peak
for i in {6..1}
do
        echo $i
        time bedtools getfasta -fi ~/group/genome/mm9/mm9.fasta -bed 'all.'$i'.pk.bed' > 'all.'$i'.pk.fa'
        time ~/group/software/meme/bin/fimo --o 'all_c'$i'_fimo' --text --thresh 1 ~/scratch/ctcf_auxin/target_ctcf_motif_full.meme 'all.'$i'.pk.fa' > 'all_c'$i'_fimo.txt'
        time python ~/scratch/ctcf_auxin/find_best_match_fimo.py -i 'all_c'$i'_fimo.txt' -o 'all_c'$i'_fimo.bestmatch.txt'
        time python ~/scratch/ctcf_auxin/fimo2bed.py -i 'all_c'$i'_fimo.bestmatch.txt' -o 'all_c'$i'_fimo.bestmatch.bed' -e 23
        time bedtools getfasta -s -name -fi ~/group/genome/mm9/mm9.fasta -bed 'all_c'$i'_fimo.bestmatch.bed' > 'all_c'$i'_fimo.bestmatch.fa'
done < bed_list.txt

Rscript ~/scratch/ctcf_auxin/Mclust_pk/plot_fimo_score.box.6.R 'Mclust_a/All/all_c' Mclust_all.fimo.pdf
