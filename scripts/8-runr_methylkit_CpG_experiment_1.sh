#!/usr/bin/env Rscript
#SBATCH --mem 48G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk	
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

#Rscript /uoa/home/r02hw22/Equina_Methylation_Analysis/Scripts/runr_methylkit_CpG_experiment_1.sh
�
# Need to create this directory before using script
setwd("/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results")

# Load in required packages
library(methylKit)
library(plyranges)
library(GenomicRanges)
library(qqman)

# Input list of files for each sample to use:
file.list = list("../../Data/Trimmed/Sample_3-3_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_9-9_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_13-13_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_15-15_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_16-16_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_17-17_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_22-22_D/meth_CpG_cov_reads",
		"../../Data/Trimmed/Sample_27-27_D/meth_CpG_cov_reads")


# Write methRead() function with selected options for Experimental Design                  
act_CpG = methRead(file.list,
                   sample.id = list("3","9","13","15","16","17","22","27"),
                   assembly="actinia",
                   pipeline = "bismarkCoverage",
                   treatment = c(1,0,0,1,0,0,1,1),
                   context = "CpG",
                   dbtype="tabix")

# Print the Raw methylation plot and save
png("raw_methylation_plots_act_CpG_experiment_1.png")
lapply(act_CpG, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save
png("raw_coverage_plots_act_CpG_experiment_1.png")
lapply(act_CpG, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

#2. Filter and Unite Methylation Results ----
# Filter the methylation results by a coverage threshold
act_CpG_f = filterByCoverage(act_CpG, lo.count = 10, hi.count = 32, suffix = "CpG_f")

# Unite the methylation results
act_CpG_fu = methylKit::unite(act_CpG_f, destrand=FALSE, min.per.group = 3L, suffix = "CpG_fu3")

# flatfile located at: /uoa/home/r02hw22/Equina_Methylation_Analysis/outputs/methylDB_2024-04-18_4SC/methylBase_CpG_fu3.txt.bgz
# Warning message:
# In .formatShortSlotValues(tmp) : NAs introduced by coercion

# When I changed the treatments to binary 0 or 1, this error was resolved

# Plot methylation profile (can choose specific if desired)
# png("meth_profile_act_CpG_experiment_1.png")
# plotMethylationProfile(act_CpG_fu)
# dev.off()

# for some reason this function isn't recognised, so I will skip it for now...

# Plot clustering of filtered data
png("cluster_plot_act_CpG_experiment_1.png")
clusterSamples(act_CpG_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

png("PCA_1-2_plot_act_CpG_multi_experiment_1.png")
PCASamples(act_CpG_fu, comp =c(1,2))
dev.off()

png("PCA_3-4_plot_act_CpG_multi_experiment_1.png")
PCASamples(act_CpG_fu, comp =c(3,4))
dev.off()

png("PCA_5-6_plot_act_CpG_multi_experiment_1.png")
PCASamples(act_CpG_fu, comp =c(5,6))
dev.off()

png("PCA_7-8_plot_act_CpG_multi_experiment_1.png")
PCASamples(act_CpG_fu, comp =c(7,8))
dev.off()

#png("PCA_9-10_plot_act_CpG_multi_experiment_1.png")
#PCASamples(act_CpG_fu, comp =c(9,10))
#dev.off()
# This last PCA doesn't seem to work, likely because there arent PC9 or 10 to calculate and plot

# Hierarchical clustering of samples
methMatrix <- getMethylationStats(act_CpG_fu, type="base")
distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
hc <- hclust(distMatrix)
plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))


# 3. Calculate Differential Methylation ----
dmb_CpG = calculateDiffMeth(act_CpG_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CpGfu3odMNtestC")

# Need to change the treatment labelling of the experiments to binary 0 and 1 for control and treatment     
# Doing this solved this issue, and a previous one. Also improved the images with addition of colour                       
                           
# Scatter plot of differential methylation
## Adjust the difference value to change the percentage of the difference
diffMethData <- getMethylDiff(dmb_CpG, difference=5, qvalue=0.1)
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")

# Error in as.double(y) :
  # cannot coerce type 'S4' to vector of type 'double'


# 4. Extract Methylation Difference Results ----
act_diff_CpG = getMethylDiff(dmb_CpG, difference = 5, qvalue=0.99, type="all")
## having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over
## Adjust the difference value to change the percentage of the difference

save(act_diff_CpG, file = "Actinia_exp1_DMBs_d5_CpG_experiment_1.RData")

png("act_diff_chrom_CpG_experiment_1.png")
diffMethPerChr(act_diff_CpG,
               plot = T,
               qvalue.cutoff = 0.99,
               meth.cutoff = 5,
               exclude = F,
               keep.empty.chrom = T)
dev.off()

# Warning message:
# In eval(quote(list(...)), env) : NAs introduced by coercion


# Does this work? do we have chromosone position data? I thought we didnt?

#5. Get Data from results ----
act_diff_CpG_data = getData(act_diff_CpG)
save(act_diff_CpG_data, file = "act_diff_CpG_data.RData")

# 6. Change results to a GRanges file
act_diff_CpG_gr =  as(act_diff_CpG_data,"GRanges")
save(act_diff_CpG_gr, file = "act_diff_CpG_gr.RData")

# 7. Overlap Granges positions with Genome annotated positions ----
genome = read_gff3("../../Data/A_Equina.gff3")

act_diff_CpG_gr_ann = join_overlap_left_directed(act_diff_CpG_gr, genome)

act_diff_CpG_gr_ann = join_overlap_left_directed(act_diff_CpG_gr, annotated_genome_gr)

act_diff_CpG_gr_ann2 = join_overlap_left_directed(annotated_genome_gr, act_diff_CpG_gr)

# When switching this to the seemingly more appropriate .ggf3 file, 
# now having a error message:
# Error: subscript contains NAs
# In addition: Warning message:
# In .merge_two_Seqinfo_objects(x, y) :
# The 2 combined objects have no sequence levels in common. (Use
# suppressWarnings() to suppress this warning.)

# I wonder if this has to do with a mismatch of names again?

act_diff_CpG_df_ann = as(act_diff_CpG_gr_ann, "data.frame")
act_diff_CpG_df_ann2 = as(act_diff_CpG_gr_ann2, "data.frame")

# 8. Create new variables and store key data ----
act_diff_CpG_df_ann$chrom_num = as.numeric(sapply(strsplit(as.character(act_diff_CpG_df_ann$seqnames), "_"), "[[",3))
# This probably won't work since we dont have chromosome data: ff_CpG_df_ann$seqnames), "_"), "[[",3))
# Error in FUN(X[[i]], ...) : subscript out of bounds

act_diff_CpG_df_ann$p_fdr = p.adjust(act_diff_CpG_df_ann$pvalue, method = "fdr")
act_diff_CpG_df_ann2$p_fdr = p.adjust(act_diff_CpG_df_ann2$pvalue, method = "fdr")
overlap_data$p_fdr <- p.adjust(overlap_data$pvalue, method = "fdr")

# This seemed to work just fine

act_diff_CpG_df_ann$class = sapply(strsplit(act_diff_CpG_df_ann$Target, "/"), "[[",1)
# This didnt work: Error in h(simpleError(msg, call)) :
# error in evaluating the argument 'X' in selecting a method for function 'sapply': non-character argument

DMBs_CpG_p01 = as.data.frame(filter(act_diff_CpG_df_ann, p_fdr < 0.1) %>% arrange(meth.diff))
DMBs_CpG_p005 = as.data.frame(filter(act_diff_CpG_df_ann, p_fdr < 0.05) %>% arrange(meth.diff))

DMBs2_CpG_p01 = as.data.frame(filter(act_diff_CpG_df_ann2, p_fdr < 0.1) %>% arrange(meth.diff))

overlap_DMBs_CpG_p01 <- as.data.frame(filter(overlap_data, p_fdr < 0.1) %>% arrange(meth_data))

write.table(DMBs_CpG, "Exp1 DMBs in CpG context.txt", sep="\t", row.names=FALSE) # have a look
write.table(overlap_DMBs_CpG_p01, "Exp1 overlap_DMBs in CpG context.txt", sep="\t", row.names=FALSE)
write.table(DMBs2_CpG_p01, "Exp1 DMBs2 in CpG context.txt", sep="\t", row.names=FALSE)
# There isn't any data here, something in the steps before didnt work, which is evident


# 9. Draw Manhattan plots
## Doesnt work for this data currently, need to work out a way to make the WHPX... groupings into arbitrary chromosome numbers currently.

'
png("act_diff_CpG_manhattan_experiment_1.png")
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

png("act_diff_CpG_manmethdiff_experiment_1.png")
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
'
