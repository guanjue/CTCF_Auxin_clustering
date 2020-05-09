cd /storage/home/gzx103/scratch/ctcf_auxin/Mclust_pk/Decay_a/All

cp Unstable.all.lowsig.bed all.1.pk.bed
cp Unstable.all.highsig.bed all.2.pk.bed
cp all.stable.bed all.3.pk.bed

### get CTCF motif centered peak
for i in {3..1}
do
        echo $i
        time bedtools getfasta -fi ~/group/genome/mm9/mm9.fasta -bed 'all.'$i'.pk.bed' > 'all.'$i'.pk.fa'
        time ~/group/software/meme/bin/fimo --o 'all_c'$i'_fimo' --text --thresh 1 ~/scratch/ctcf_auxin/target_ctcf_motif_full.meme 'all.'$i'.pk.fa' > 'all_c'$i'_fimo.txt'
        time python ~/scratch/ctcf_auxin/find_best_match_fimo.py -i 'all_c'$i'_fimo.txt' -o 'all_c'$i'_fimo.bestmatch.txt'
        time python ~/scratch/ctcf_auxin/fimo2bed.py -i 'all_c'$i'_fimo.bestmatch.txt' -o 'all_c'$i'_fimo.bestmatch.bed' -e 23
        time bedtools getfasta -s -name -fi ~/group/genome/mm9/mm9.fasta -bed 'all_c'$i'_fimo.bestmatch.bed' > 'all_c'$i'_fimo.bestmatch.fa'
done


cd /storage/home/gzx103/scratch/ctcf_auxin/Mclust_pk/Decay_a/Matched

cp unstable_highsig.sampled.run1.bed match.1.pk.bed
cp persistent.sampled.run1.bed match.2.pk.bed
cp all.stable.bed match.3.pk.bed

### get CTCF motif centered peak
for i in {3..1}
do
        echo $i
        time bedtools getfasta -fi ~/group/genome/mm9/mm9.fasta -bed 'match.'$i'.pk.bed' > 'match.'$i'.pk.fa'
        time ~/group/software/meme/bin/fimo --o 'match_c'$i'_fimo' --text --thresh 1 ~/scratch/ctcf_auxin/target_ctcf_motif_full.meme 'match.'$i'.pk.fa' > 'match_c'$i'_fimo.txt'
        time python ~/scratch/ctcf_auxin/find_best_match_fimo.py -i 'match_c'$i'_fimo.txt' -o 'match_c'$i'_fimo.bestmatch.txt'
        time python ~/scratch/ctcf_auxin/fimo2bed.py -i 'match_c'$i'_fimo.bestmatch.txt' -o 'match_c'$i'_fimo.bestmatch.bed' -e 23
        time bedtools getfasta -s -name -fi ~/group/genome/mm9/mm9.fasta -bed 'match_c'$i'_fimo.bestmatch.bed' > 'match_c'$i'_fimo.bestmatch.fa'
done



cd /storage/home/gzx103/scratch/ctcf_auxin/Mclust_pk/Mclust_a/Matched

cp unstable_highsig.sampled.run1.bed match.1.pk.bed
cp persistent.sampled.run1.bed match.2.pk.bed
cp all.stable.bed match.3.pk.bed

### get CTCF motif centered peak
for i in {3..1}
do
        echo $i
        time bedtools getfasta -fi ~/group/genome/mm9/mm9.fasta -bed 'match.'$i'.pk.bed' > 'match.'$i'.pk.fa'
        time ~/group/software/meme/bin/fimo --o 'match_c'$i'_fimo' --text --thresh 1 ~/scratch/ctcf_auxin/target_ctcf_motif_full.meme 'match.'$i'.pk.fa' > 'match_c'$i'_fimo.txt'
        time python ~/scratch/ctcf_auxin/find_best_match_fimo.py -i 'match_c'$i'_fimo.txt' -o 'match_c'$i'_fimo.bestmatch.txt'
        time python ~/scratch/ctcf_auxin/fimo2bed.py -i 'match_c'$i'_fimo.bestmatch.txt' -o 'match_c'$i'_fimo.bestmatch.bed' -e 23
        time bedtools getfasta -s -name -fi ~/group/genome/mm9/mm9.fasta -bed 'match_c'$i'_fimo.bestmatch.bed' > 'match_c'$i'_fimo.bestmatch.fa'
done



