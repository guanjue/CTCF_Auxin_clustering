cd ~/scratch/ctcf_auxin

### step0 get new blacklist
time bash blacklist/get_blacklist.sh

### step1 get ctcf pk
time bash ctcf_auxin.step1.getpk.sh

### step2 get raw rc bw
time bash bam_bw_bedgraph/qsub_get_raw_signal.1bp.allctrl.sh

### step3 get signal matrix
time bash get_signal_mat.sh

### step4 get bg substracted signal matrix
time Rscript get_bg_subtracted_sigmat.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.0kb.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt
#time Rscript get_bg_subtracted_sigmat.macs.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.1kb.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.5kb.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.10kb.sigmat.txt ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt

### step5 normalize to qPCR
time Rscript qPCR_LMnorm.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt bam_bw_bedgraph/file_list.all.timepoints.txt 20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.sig.txt 20180910_CTCF_ChIPqPCR_7pts.rmNA.qPCR.timepoints.txt ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.txt

### get normalized bw
time bash get_bw_norm.sh

### step6 get R2 before and after normalization
time Rscript get_replicates_R2.R

### step7 get Ward Statistics
time Rscript /Users/universe/Documents/2020_BG/ctcf_auxin/02_22_2020_pipeline/get_Wardstat.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.txt ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.Ward.txt
time Rscript /Users/universe/Documents/2020_BG/ctcf_auxin/02_22_2020_pipeline/get_Wardstat.R ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.txt ctcf.qPCR.randbg.blackrm.idsort.sigmat.bgsub.LMqPCRnorm.Ward.txt

### step8 clustering
mkdir all_pk
mkdir all_pk/LMqPCRnorm_decay_folder
time Rscript get_cluster_Decay_curve.R

### Step9 get same 0A peak sets (stable vs unstable)
time Rscript get_0A_same_sus.pk.R

### Step10 FIMO
time bash get_fimo_analysis.sh

### Step11 Mclust for validation
time Rscript get_cluster_Mclust.R

### Step12 plot exponential signal trend matching initial signal
time Rscript intensity_matching_JL_final.R
time Rscript intensity_matching_decay_JL_final.R



