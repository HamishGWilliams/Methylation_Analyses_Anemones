#!/usr/bin/env Rscript
#SBATCH --mem 24G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk	
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err



#Rscript /uoa/home/r02hw22/Equina_Methylation_Analysis/Scripts/runr_methylkit_CpG_experiment_2.sh

# Need to create this directory before using script
setwd("/uoa/home/r02hw22/Equina_Methylation_Analysis/outputs")

# Load in required packages
library(methylKit)
library(plyranges)
library(GenomicRanges)
library(qqman)
library(dplyr)
library(ggplot2)

# Input list of files for each sample to use:
file.list = list("../Data/Trimmed/Sample_6-6_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_7-7_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_8-8_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_10-10_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_14-14_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_20-20_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_21-21_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_29-11_redo_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_32-18_redo_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_36-2_redo_D/meth_CpG_cov_reads")


# Write methRead() function with selected options for Experimental Design                  
act_CpG = methRead(file.list,
                sample.id = list("6","7","8","10","14","20","21","11","18","2"),
                assembly="actinia",
                pipeline = "bismarkCoverage",
                treatment = c("B","B","C","B","C","C","C","B","B","C"),
                context = "CpG",
                dbtype="tabix")

# Print the methylation plot and save
png("methylation_plots_act_CpG_experiment_2.png")
lapply(act_CpG, getMethylationStats,  plot=T)
dev.off()

# Print the Coverage plot and save
png("coverage_plots_act_CpG_experiment_2.png")
lapply(act_CpG, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# Filter the methylation results by a coverage threshold
act_CpG_f = filterByCoverage(act_CpG, hi.count = 30, suffix = "CpG_f")

# Unite the 
act_CpG_fu = methylKit::unite(act_CpG_f, destrand=FALSE, min.per.group = 3L, suffix = "CpG_fu3")

png("cluster_plot_act_CpG_experiment_2.png")
clusterSamples(act_CpG_fu,
                        filterByQuantile = F,
                        sd.threshold = 0.1,
                        dist="correlation",
                        method="ward.D",
                        plot=T)
dev.off()

png("PCA_plot_act_CpG_multi_experiment_2.png")  #note this takes a very long time as it runs the PCA each time
PCASamples(act_CpG_fu, comp =c(1,2))
PCASamples(act_CpG_fu, comp =c(3,4))
PCASamples(act_CpG_fu, comp =c(5,6))
PCASamples(act_CpG_fu, comp =c(7,8))
PCASamples(act_CpG_fu, comp =c(9,10))
dev.off()


dmb_CpG = calculateDiffMeth(act_CpG_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CpGfu3odMNtestC")
                            

act_diff_CpG = getMethylDiff(dmb_CpG, difference = 5, qvalue=0.99, type="all")
#having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

save(act_diff_CpG, file = "Actinia_exp1_DMBs_d5_CpG_experiment_2.RData")

png("act_diff_chrom_CpG_experiment_2.png")
diffMethPerChr(act_diff_CpG,
                        plot = T,
                        qvalue.cutoff = 0.1,
                        meth.cutoff = 5,
                        exclude = F,
                        keep.empty.chrom = T)
dev.off()

act_diff_CpG_data = getData(act_diff_CpG)

act_diff_CpG_gr =  as(act_diff_CpG_data,"GRanges") %>% filter( ! (seqnames %in% c("HiC_scaffold_6" , "HiC_scaffold_16")))

genome = read_gff3("../../Genomes/DUM/trimmed/DUM_hifi_hic_scaffolded_trim.fa.out.gff3")     #replace with location of the genome

act_diff_CpG_gr_ann = join_overlap_left_directed(act_diff_CpG_gr, genome)

act_diff_CpG_df_ann = as(act_diff_CpG_gr_ann, "data.frame")

act_diff_CpG_df_ann$chrom_num = as.numeric(sapply(strsplit(as.character(act_diff_CpG_df_ann$seqnames), "_"), "[[",3))
act_diff_CpG_df_ann$p_fdr = p.adjust(act_diff_CpG_df_ann$pvalue, method = "fdr")
act_diff_CpG_df_ann$class = sapply(strsplit(act_diff_CpG_df_ann$Target, "/"), "[[",1)

DMBs_CpG = as.data.frame(filter(act_diff_CpG_df_ann, p_fdr < 0.05) %>% arrange(meth.diff))

write.table(DMBs_CpG, "Exp1 DMBs in CpG context.txt", sep="\t", row.names=FALSE) # have a look


png("act_diff_CpG_manhattan_experiment_2.png")
manhattan(act_diff_CpG_df_ann,
                           chr="chrom_num",
                           bp="start",
                           snp="Target",
                           p="pvalue",
                           col=c("dodgerblue","coral" ), 
                           suggestiveline = -log10(1e-05),
                           genomewideline = -log10(5e-08),,
                           annotateTop = TRUE,
                           main = "Manhattan plot of DMB p values CpG context")
dev.off()

png("act_diff_CpG_manmethdiff_experiment_2.png")
manhattan(act_diff_CpG_df_ann, p = "meth.diff",
                           chr="chrom_num",
                           bp="start",
                           snp="Target",
                           logp = FALSE,
                           ylab = "Difference in methylation percentage",
                           col=c("dodgerblue","coral" ),
                           genomewideline = FALSE,
                           suggestiveline = FALSE,
                           main = "Manhattan plot of methylation % differences CpG context")
dev.off()
