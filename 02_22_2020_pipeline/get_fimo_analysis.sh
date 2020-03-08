cd ~/scratch/ctcf_auxin/all_pk/LMqPCRnorm_decay_folder

cp balanced_peaklist_for_compare.unstable.submatch.lowsig.txt balanced_peaklist_for_compare.1.pk.clean.txt
cp balanced_peaklist_for_compare.unstable.submatch.highsig.txt balanced_peaklist_for_compare.2.pk.clean.txt
cp balanced_peaklist_for_compare.stable.submatch.highsig.txt balanced_peaklist_for_compare.3.pk.clean.txt

### get CTCF motif centered peak
for i in {3..1}
do
        echo $i
        tail -n+2 'balanced_peaklist_for_compare.'$i'.pk.clean.txt' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3}' > 'balanced_peaklist_for_compare.'$i'.pk.clean.bed'
        time bedtools getfasta -fi ~/group/genome/mm9/mm9.fasta -bed 'balanced_peaklist_for_compare.'$i'.pk.clean.bed' > 'balanced_peaklist_for_compare.'$i'.pk.fa'
        time ~/group/software/meme/bin/fimo --o 'ctcf2_c'$i'_fimo' --text --thresh 1 ~/scratch/ctcf_auxin/target_ctcf_motif_full.meme 'balanced_peaklist_for_compare.'$i'.pk.fa' > 'ctcf2_c'$i'_fimo.txt'
        time python ~/scratch/ctcf_auxin/find_best_match_fimo.py -i 'ctcf2_c'$i'_fimo.txt' -o 'ctcf2_c'$i'_fimo.bestmatch.txt'
        time python ~/scratch/ctcf_auxin/fimo2bed.py -i 'ctcf2_c'$i'_fimo.bestmatch.txt' -o 'ctcf2_c'$i'_fimo.bestmatch.bed' -e 23
        time bedtools getfasta -s -name -fi ~/group/genome/mm9/mm9.fasta -bed 'ctcf2_c'$i'_fimo.bestmatch.bed' > 'ctcf2_c'$i'_fimo.bestmatch.fa'
done


### get CTCF motif catevec
for i in {3..1}
do
        echo $i
        time cat 'ctcf2_c'$i'_fimo.bestmatch.bed' | awk -F '\t' -v OFS='\t' -v group=$i '{print $1":"int($2/2+$3/2), group}' > 'ctcf2_c'$i'_fimo.mid.txt'
        time python ~/scratch/ctcf_auxin/fimo2bed.py -i 'ctcf2_c'$i'_fimo.bestmatch.txt' -o 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed' -e 23
        time bedtools getfasta -s -name -fi ~/group/genome/mm9/mm9.fasta -bed 'ctcf2_c'$i'_fimo.bestmatch.5exp.bed' > 'ctcf2_c'$i'_fimo.bestmatch.5exp.fa'
        time python ~/scratch/ctcf_auxin/fasta2categoricalvec.py -i 'ctcf2_c'$i'_fimo.bestmatch.5exp.fa'
done


### get random peak result
rand_bed='/storage/home/gzx103/scratch/ctcf_auxin/ctcf_master_pk/mm9.500bp.randomwin.100000.bg.sort.bed'
shuf -n 100000 $rand_bed | head -10000 > $rand_bed'.rand10k.bed'
time bedtools getfasta -fi ~/group/genome/mm9/mm9.fasta -bed $rand_bed'.rand10k.bed' > $rand_bed'.rand10k.fa'
time ~/group/software/meme/bin/fimo --o ctcf2_fimo_randbed --text --thresh 1 ~/scratch/ctcf_auxin/target_ctcf_motif_full.meme $rand_bed'.rand10k.fa' > ctcf2_fimo_randbed.txt


### count the random pk results
for i in {1..3}
do
cat 'ctcf2_c'$i'_fimo.txt' | awk -F '\t' -v OFS='\t' '{if ($8<0.00001) print $0}' | wc -l
done
cat ctcf2_fimo_randbed.txt | awk -F '\t' -v OFS='\t' '{if ($8<0.00001) print $0}' | wc -l
### count peak number
for i in {1..3}
do
wc -l 'balanced_peaklist_for_compare.'$i'.pk.clean.bed'
done


~/group/software/ucsc/liftOver balanced_peaklist_for_compare.stable.submatch.highsig.bed ~/group/genome/mm9/mm9ToMm10.over.chain.gz stable.submatch.highsig.mm10.bed stable.submatch.highsig.mm10unmap.bed
~/group/software/ucsc/liftOver balanced_peaklist_for_compare.unstable.submatch.highsig.bed ~/group/genome/mm9/mm9ToMm10.over.chain.gz unstable.submatch.highsig.mm10.bed unstable.submatch.highsig.mm10unmap.bed
~/group/software/ucsc/liftOver balanced_peaklist_for_compare.unstable.submatch.lowsig.bed ~/group/genome/mm9/mm9ToMm10.over.chain.gz unstable.submatch.lowsig.mm10.bed unstable.submatch.lowsig.mm10unmap.bed

