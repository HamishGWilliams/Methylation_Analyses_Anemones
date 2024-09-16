#!/usr/bin/env Rscript
#SBATCH --mem 24G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk	
#SBATCH --time=1-00:00:00

# Need to create this directory before using script
setwd("/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/outputs_promoters")

# Load in required packages
library(methylKit)
library(plyranges)
library(GenomicRanges)
library(qqman)
library(dplyr)
library(ggplot2)

## ------------------------ Experiment 1 ----------------------- ##

	  #  ---------- CpG ------------ #

# Input list of files for each sample to use:
file.list = list("../Data/Trimmed/Sample_3-3_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_9-9_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_13-13_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_15-15_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_16-16_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_17-17_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_22-22_D/meth_CpG_cov_reads",
		"../Data/Trimmed/Sample_27-27_D/meth_CpG_cov_reads")


# Write methRead() function with selected options for Experimental Design                  
act_CpG = methRead(file.list,
                   sample.id = list("3","9","13","15","16","17","22","27"),
                   assembly="actinia",
                   pipeline = "bismarkCoverage",
                   treatment = c(1,0,0,1,0,0,1,1),
                   context = "CpG",
                   dbtype="tabix")

# Print the Raw methylation plot and save - PLOT
png("raw_methylation_plots_act_CpG_experiment_2.png")
lapply(act_CpG, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save - PLOT
png("raw_coverage_plots_act_CpG_experiment_2.png")
lapply(act_CpG, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# 3. Filter the methylation results by a coverage threshold
act_CpG_f = filterByCoverage(act_CpG, lo.count = 10, hi.count = 32, suffix = "CpG_f")
				      # Change lo.count & hi.count based on the methylation & coverage plots

# 4. Unite the methylation results
act_CpG_fu = methylKit::unite(act_CpG_f, destrand=FALSE, min.per.group = 3L, suffix = "CpG_fu3")

# Plot clustering of filtered data
png("cluster_plot_act_CpG_experiment_1.png")
clusterSamples(act_CpG_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

# Plot PCAs
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

# Hierarchical clustering of samples
# methMatrix <- getMethylationStats(act_CpG_fu, type="base")
# distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
# hc <- hclust(distMatrix)
# plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))

# !! Currently not working, states:
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function �getMethylationStats� for signature �"methylBaseDB"�

# 5. Calculate Differential Methylation ----
dmb_CpG = calculateDiffMeth(act_CpG_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CpGfu3odMNtestC")

# Scatter plot of differential methylation
diffMethData <- getMethylDiff(dmb_CpG, qvalue=0.99)

png("scatter_plot_methylation_exp1_CpG.png")
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")
dev.off()

# !! Error !!
# Error in as.double(y) :
# cannot coerce type 'S4' to vector of type 'double'

# 6. Extract Methylation Difference Results ----
act_diff_CpG = getMethylDiff(dmb_CpG, qvalue=0.99, type="all")
# having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

# !! CHECKPOINT !! save work to RData file
save(act_diff_CpG, file = "Actinia_exp1_DMBs_d5_CpG_experiment_1.RData")

# 7. Get Data from results ----
act_diff_CpG_data = getData(act_diff_CpG)
# !! CHECKPOINT !! - save methylation data
save(act_diff_CpG_data, file = "act_diff_CpG_data_exp1.RData")

# 8. Change results to a GRanges file
act_diff_CpG_gr =  as(act_diff_CpG_data,"GRanges")
# !! CHECKPOINT !! - save methylation data as a GRanges file
save(act_diff_CpG_gr, file = "act_diff_CpG_gr_exp1.RData")

# 9. Load in the annotated genome file
methylation_matching <- read.table("methylation_matching_file.tsv", header = FALSE, sep = "\t", col.names = c("chr", "utg_name"))

# 10. Rename columns for clarity
colnames(act_diff_CpG_data)[1] <- "chr"
colnames(methylation_matching) <- c("chr", "utg_name")

# 11. Merge the data frames by the "chr" column
merged_data <- merge(act_diff_CpG_data, methylation_matching, by = "chr", all.x = TRUE)

# 12. Perform a left join to merge data frames
merged_data <- left_join(act_diff_CpG_data, methylation_matching, by = "chr")

# 13. Rename 'utg_name' to 'seqnames' and 'chr' to 'WHPX_names', then rearrange the columns
act_diff_CpG_data_matched <- merged_data %>%
  rename(chr = utg_name, WHPX_names = chr) %>%
  select(chr, everything(), WHPX_names)

# 14. remove the WHPX names
act_diff_CpG_data_match <- subset(act_diff_CpG_data_matched[,-2])

# 15. Transform to GRanges format
act_diff_CpG_gr <- as(act_diff_CpG_data_matched, "GRanges")

# 16. load in annotated genome
genome <- read_gff3("../Data/A_Equina_with_promoters.gff3")
# 17. Change to dataframe object
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
# 18. convert back to GRanges format
annotated_genome_gr <- as(genome, "GRanges")

# 19. Join_overlap_left_directed() on the methylation and genome files
act_diff_CpG_gr_ann = join_overlap_left_directed(act_diff_CpG_gr, annotated_genome_gr)

# 20. change to df format
act_diff_CpG_df_ann = as(act_diff_CpG_gr_ann, "data.frame")

# 21. Calculate p-adjusted values (bonferonni-hochberg)
act_diff_CpG_df_ann$p_fdr = p.adjust(act_diff_CpG_df_ann$pvalue, method = "fdr")

# Correct variable format
act_diff_CpG_df_ann$Parent <- as.character(act_diff_CpG_df_ann$Parent)

# 22. Write table of results - !! CHECKPOINT !!
write.table(act_diff_CpG_df_ann, "Exp1 DMBs in CpG context.txt", sep="\t", row.names=FALSE)

# 23. Volcano plot of differential methylation

data <- read.table("Exp1 DMBs in CpG context.txt", sep = "\t", header = T)

data$diffmeth <- 'NO'
data$diffmeth[data$meth.diff >= (5) & data$p_fdr < 0.1] <- 'UP'
data$diffmeth[data$meth.diff <= (-5) & data$p_fdr < 0.1] <- 'DOWN'
data$diffmeth <- as.factor(data$diffmeth)

volcano_plot <- ggplot(data = data,
       aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) + 
  geom_point() + 
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") + 
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") + 
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_color_manual(values = c("deepskyblue","gray", "brown1"),
                     labels = c("Down Methylated", "Not Significant", "Up Methylated")) + 
  labs(
    x = "Differential methyltion %", y = expression("-log"[10]*"p-adj")
  ) + 
  ggtitle("Differential Methylation of CpGs of Experiment 1")

ggsave("volcano_plot_exp2_cpg.png", volcano_plot, width = 18, height = 12, units = "cm")

	  #  ---------- CHG ------------ #

# Input list of files for each sample to use:
file.list = list("../Data/Trimmed/Sample_3-3_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_9-9_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_13-13_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_15-15_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_16-16_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_17-17_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_22-22_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_27-27_D/meth_CHG_cov_reads")


# Write methRead() function with selected options for Experimental Design                  
act_CHG = methRead(file.list,
                   sample.id = list("3","9","13","15","16","17","22","27"),
                   assembly="actinia",
                   pipeline = "bismarkCoverage",
                   treatment = c(1,0,0,1,0,0,1,1),
                   context = "CHG",
                   dbtype="tabix")

# Print the Raw methylation plot and save - PLOT
png("raw_methylation_plots_act_CHG_experiment_2.png")
lapply(act_CHG, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save - PLOT
png("raw_coverage_plots_act_CHG_experiment_2.png")
lapply(act_CHG, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# 3. Filter the methylation results by a coverage threshold
act_CHG_f = filterByCoverage(act_CHG, lo.count = 10, hi.count = 32, suffix = "CHG_f")
				      # Change lo.count & hi.count based on the methylation & coverage plots

# 4. Unite the methylation results
act_CHG_fu = methylKit::unite(act_CHG_f, destrand=FALSE, min.per.group = 3L, suffix = "CHG_fu3")

# Plot clustering of filtered data
png("cluster_plot_act_CHG_experiment_1.png")
clusterSamples(act_CHG_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

# Plot PCAs
png("PCA_1-2_plot_act_CHG_multi_experiment_1.png")
PCASamples(act_CHG_fu, comp =c(1,2))
dev.off()

png("PCA_3-4_plot_act_CHG_multi_experiment_1.png")
PCASamples(act_CHG_fu, comp =c(3,4))
dev.off()

png("PCA_5-6_plot_act_CHG_multi_experiment_1.png")
PCASamples(act_CHG_fu, comp =c(5,6))
dev.off()

png("PCA_7-8_plot_act_CHG_multi_experiment_1.png")
PCASamples(act_CHG_fu, comp =c(7,8))
dev.off()

# Hierarchical clustering of samples
# methMatrix <- getMethylationStats(act_CHG_fu, type="base")
# distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
# hc <- hclust(distMatrix)
# plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))

# !! Currently not working, states:
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function �getMethylationStats� for signature �"methylBaseDB"�

# 5. Calculate Differential Methylation ----
dmb_CHG = calculateDiffMeth(act_CHG_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CHGfu3odMNtestC")

# Scatter plot of differential methylation
diffMethData <- getMethylDiff(dmb_CHG, qvalue=0.99)

png("scatter_plot_methylation_exp1_CHG.png")
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")
dev.off()

# !! Error !!
# Error in as.double(y) :
# cannot coerce type 'S4' to vector of type 'double'

# 6. Extract Methylation Difference Results ----
act_diff_CHG = getMethylDiff(dmb_CHG, qvalue=0.99, type="all")
# having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

# !! CHECKPOINT !! save work to RData file
save(act_diff_CHG, file = "Actinia_exp1_DMBs_d5_CHG_experiment_1.RData")

# 7. Get Data from results ----
act_diff_CHG_data = getData(act_diff_CHG)
# !! CHECKPOINT !! - save methylation data
save(act_diff_CHG_data, file = "act_diff_CHG_data_exp1.RData")

# 8. Change results to a GRanges file
act_diff_CHG_gr =  as(act_diff_CHG_data,"GRanges")
# !! CHECKPOINT !! - save methylation data as a GRanges file
save(act_diff_CHG_gr, file = "act_diff_CHG_gr_exp1.RData")

# 9. Load in the annotated genome file
methylation_matching <- read.table("methylation_matching_file.tsv", header = FALSE, sep = "\t", col.names = c("chr", "utg_name"))

# 10. Rename columns for clarity
colnames(act_diff_CHG_data)[1] <- "chr"
colnames(methylation_matching) <- c("chr", "utg_name")

# 11. Merge the data frames by the "chr" column
merged_data <- merge(act_diff_CHG_data, methylation_matching, by = "chr", all.x = TRUE)

# 12. Perform a left join to merge data frames
merged_data <- left_join(act_diff_CHG_data, methylation_matching, by = "chr")

# 13. Rename 'utg_name' to 'seqnames' and 'chr' to 'WHPX_names', then rearrange the columns
act_diff_CHG_data_matched <- merged_data %>%
  rename(chr = utg_name, WHPX_names = chr) %>%
  select(chr, everything(), WHPX_names)

# 14. remove the WHPX names
act_diff_CHG_data_match <- subset(act_diff_CHG_data_matched[,-2])

# 15. Transform to GRanges format
act_diff_CHG_gr <- as(act_diff_CHG_data_matched, "GRanges")

# 16. load in annotated genome
genome <- read_gff3("../Data/A_Equina_with_promoters.gff3")
# 17. Change to dataframe object
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
# 18. convert back to GRanges format
annotated_genome_gr <- as(genome, "GRanges")

# 19. Join_overlap_left_directed() on the methylation and genome files
act_diff_CHG_gr_ann = join_overlap_left_directed(act_diff_CHG_gr, annotated_genome_gr)

# 20. change to df format
act_diff_CHG_df_ann = as(act_diff_CHG_gr_ann, "data.frame")

# 21. Calculate p-adjusted values (bonferonni-hochberg)
act_diff_CHG_df_ann$p_fdr = p.adjust(act_diff_CHG_df_ann$pvalue, method = "fdr")

# Correct variable format
act_diff_CHG_df_ann$Parent <- as.character(act_diff_CHG_df_ann$Parent)

# 22. Write table of results - !! CHECKPOINT !!
write.table(act_diff_CHG_df_ann, "Exp1 DMBs in CHG context.txt", sep="\t", row.names=FALSE)

# 23. Volcano plot of differential methylation

data <- read.table("Exp1 DMBs in CHG context.txt", sep = "\t", header = T)

data$diffmeth <- 'NO'
data$diffmeth[data$meth.diff >= (5) & data$p_fdr < 0.1] <- 'UP'
data$diffmeth[data$meth.diff <= (-5) & data$p_fdr < 0.1] <- 'DOWN'
data$diffmeth <- as.factor(data$diffmeth)

volcano_plot <- ggplot(data = data,
       aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) + 
  geom_point() + 
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") + 
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") + 
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_color_manual(values = c("deepskyblue","gray", "brown1"),
                     labels = c("Down Methylated", "Not Significant", "Up Methylated")) + 
  labs(
    x = "Differential methyltion %", y = expression("-log"[10]*"p-adj")
  ) + 
  ggtitle("Differential Methylation of CHGs of Experiment 1")

ggsave("volcano_plot_exp2_CHG.png", volcano_plot, width = 18, height = 12, units = "cm")

	  #  ---------- CHH ------------ #

# Input list of files for each sample to use:
file.list = list("../Data/Trimmed/Sample_3-3_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_9-9_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_13-13_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_15-15_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_16-16_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_17-17_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_22-22_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_27-27_D/meth_CHH_cov_reads")


# Write methRead() function with selected options for Experimental Design                  
act_CHH = methRead(file.list,
                   sample.id = list("3","9","13","15","16","17","22","27"),
                   assembly="actinia",
                   pipeline = "bismarkCoverage",
                   treatment = c(1,0,0,1,0,0,1,1),
                   context = "CHH",
                   dbtype="tabix")

# Print the Raw methylation plot and save - PLOT
png("raw_methylation_plots_act_CHH_experiment_2.png")
lapply(act_CHH, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save - PLOT
png("raw_coverage_plots_act_CHH_experiment_2.png")
lapply(act_CHH, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# 3. Filter the methylation results by a coverage threshold
act_CHH_f = filterByCoverage(act_CHH, lo.count = 10, hi.count = 32, suffix = "CHH_f")
				      # Change lo.count & hi.count based on the methylation & coverage plots

# 4. Unite the methylation results
act_CHH_fu = methylKit::unite(act_CHH_f, destrand=FALSE, min.per.group = 3L, suffix = "CHH_fu3")

# Plot clustering of filtered data
png("cluster_plot_act_CHH_experiment_1.png")
clusterSamples(act_CHH_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

# Plot PCAs
png("PCA_1-2_plot_act_CHH_multi_experiment_1.png")
PCASamples(act_CHH_fu, comp =c(1,2))
dev.off()

png("PCA_3-4_plot_act_CHH_multi_experiment_1.png")
PCASamples(act_CHH_fu, comp =c(3,4))
dev.off()

png("PCA_5-6_plot_act_CHH_multi_experiment_1.png")
PCASamples(act_CHH_fu, comp =c(5,6))
dev.off()

png("PCA_7-8_plot_act_CHH_multi_experiment_1.png")
PCASamples(act_CHH_fu, comp =c(7,8))
dev.off()

# Hierarchical clustering of samples
# methMatrix <- getMethylationStats(act_CHH_fu, type="base")
# distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
# hc <- hclust(distMatrix)
# plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))

# !! Currently not working, states:
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function �getMethylationStats� for signature �"methylBaseDB"�

# 5. Calculate Differential Methylation ----
dmb_CHH = calculateDiffMeth(act_CHH_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CHHfu3odMNtestC")

# Scatter plot of differential methylation
diffMethData <- getMethylDiff(dmb_CHH, qvalue=0.99)

png("scatter_plot_methylation_exp1_CHH.png")
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")
dev.off()

# !! Error !!
# Error in as.double(y) :
# cannot coerce type 'S4' to vector of type 'double'

# 6. Extract Methylation Difference Results ----
act_diff_CHH = getMethylDiff(dmb_CHH, qvalue=0.99, type="all")
# having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

# !! CHECKPOINT !! save work to RData file
save(act_diff_CHH, file = "Actinia_exp1_DMBs_d5_CHH_experiment_1.RData")

# 7. Get Data from results ----
act_diff_CHH_data = getData(act_diff_CHH)
# !! CHECKPOINT !! - save methylation data
save(act_diff_CHH_data, file = "act_diff_CHH_data_exp1.RData")

# 8. Change results to a GRanges file
act_diff_CHH_gr =  as(act_diff_CHH_data,"GRanges")
# !! CHECKPOINT !! - save methylation data as a GRanges file
save(act_diff_CHH_gr, file = "act_diff_CHH_gr_exp1.RData")

# 9. Load in the annotated genome file
methylation_matching <- read.table("methylation_matching_file.tsv", header = FALSE, sep = "\t", col.names = c("chr", "utg_name"))

# 10. Rename columns for clarity
colnames(act_diff_CHH_data)[1] <- "chr"
colnames(methylation_matching) <- c("chr", "utg_name")

# 11. Merge the data frames by the "chr" column
merged_data <- merge(act_diff_CHH_data, methylation_matching, by = "chr", all.x = TRUE)

# 12. Perform a left join to merge data frames
merged_data <- left_join(act_diff_CHH_data, methylation_matching, by = "chr")

# 13. Rename 'utg_name' to 'seqnames' and 'chr' to 'WHPX_names', then rearrange the columns
act_diff_CHH_data_matched <- merged_data %>%
  rename(chr = utg_name, WHPX_names = chr) %>%
  select(chr, everything(), WHPX_names)

# 14. remove the WHPX names
act_diff_CHH_data_match <- subset(act_diff_CHH_data_matched[,-2])

# 15. Transform to GRanges format
act_diff_CHH_gr <- as(act_diff_CHH_data_matched, "GRanges")

# 16. load in annotated genome
genome <- read_gff3("../Data/A_Equina_with_promoters.gff3")
# 17. Change to dataframe object
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
# 18. convert back to GRanges format
annotated_genome_gr <- as(genome, "GRanges")

# 19. Join_overlap_left_directed() on the methylation and genome files
act_diff_CHH_gr_ann = join_overlap_left_directed(act_diff_CHH_gr, annotated_genome_gr)

# 20. change to df format
act_diff_CHH_df_ann = as(act_diff_CHH_gr_ann, "data.frame")

# 21. Calculate p-adjusted values (bonferonni-hochberg)
act_diff_CHH_df_ann$p_fdr = p.adjust(act_diff_CHH_df_ann$pvalue, method = "fdr")

# Correct variable format
act_diff_CHH_df_ann$Parent <- as.character(act_diff_CHH_df_ann$Parent)

# 22. Write table of results - !! CHECKPOINT !!
write.table(act_diff_CHH_df_ann, "Exp1 DMBs in CHH context.txt", sep="\t", row.names=FALSE)

# 23. Volcano plot of differential methylation

data <- read.table("Exp1 DMBs in CHH context.txt", sep = "\t", header = T)

data$diffmeth <- 'NO'
data$diffmeth[data$meth.diff >= (5) & data$p_fdr < 0.1] <- 'UP'
data$diffmeth[data$meth.diff <= (-5) & data$p_fdr < 0.1] <- 'DOWN'
data$diffmeth <- as.factor(data$diffmeth)

volcano_plot <- ggplot(data = data,
       aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) + 
  geom_point() + 
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") + 
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") + 
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_color_manual(values = c("deepskyblue","gray", "brown1"),
                     labels = c("Down Methylated", "Not Significant", "Up Methylated")) + 
  labs(
    x = "Differential methyltion %", y = expression("-log"[10]*"p-adj")
  ) + 
  ggtitle("Differential Methylation of CHHs of Experiment 1")

ggsave("volcano_plot_exp2_CHH.png", volcano_plot, width = 18, height = 12, units = "cm")


## ------------------------ Experiment 2 ----------------------- ##

	  #  ---------- CpG ------------ #


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
                treatment = c(0,0,1,0,1,1,1,0,0,1),
                context = "CpG",
                dbtype="tabix")

# Print the Raw methylation plot and save - PLOT
png("raw_methylation_plots_act_CpG_experiment_2.png")
lapply(act_CpG, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save - PLOT
png("raw_coverage_plots_act_CpG_experiment_2.png")
lapply(act_CpG, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# 3. Filter the methylation results by a coverage threshold
act_CpG_f = filterByCoverage(act_CpG, lo.count = 10, hi.count = 32, suffix = "CpG_f")
				      # Change lo.count & hi.count based on the methylation & coverage plots

# 4. Unite the methylation results
act_CpG_fu = methylKit::unite(act_CpG_f, destrand=FALSE, min.per.group = 3L, suffix = "CpG_fu3")

# Plot clustering of filtered data
png("cluster_plot_act_CpG_experiment_2.png")
clusterSamples(act_CpG_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

# Plot PCAs
png("PCA_1-2_plot_act_CpG_multi_experiment_2.png")
PCASamples(act_CpG_fu, comp =c(1,2))
dev.off()

png("PCA_3-4_plot_act_CpG_multi_experiment_2.png")
PCASamples(act_CpG_fu, comp =c(3,4))
dev.off()

png("PCA_5-6_plot_act_CpG_multi_experiment_2.png")
PCASamples(act_CpG_fu, comp =c(5,6))
dev.off()

png("PCA_7-8_plot_act_CpG_multi_experiment_2.png")
PCASamples(act_CpG_fu, comp =c(7,8))
dev.off()

# Hierarchical clustering of samples
# methMatrix <- getMethylationStats(act_CpG_fu, type="base")
# distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
# hc <- hclust(distMatrix)
# plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))

# !! Currently not working, states:
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function �getMethylationStats� for signature �"methylBaseDB"�

# 5. Calculate Differential Methylation ----
dmb_CpG = calculateDiffMeth(act_CpG_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CpGfu3odMNtestC")

# Scatter plot of differential methylation
diffMethData <- getMethylDiff(dmb_CpG, qvalue=0.99)

png("scatter_plot_methylation_exp1_CpG.png")
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")
dev.off()

# !! Error !!
# Error in as.double(y) :
# cannot coerce type 'S4' to vector of type 'double'

# 6. Extract Methylation Difference Results ----
act_diff_CpG = getMethylDiff(dmb_CpG, qvalue=0.99, type="all")
# having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

# !! CHECKPOINT !! save work to RData file
save(act_diff_CpG, file = "Actinia_exp1_DMBs_d5_CpG_experiment_2.RData")

# 7. Get Data from results ----
act_diff_CpG_data = getData(act_diff_CpG)
# !! CHECKPOINT !! - save methylation data
save(act_diff_CpG_data, file = "act_diff_CpG_data_exp2.RData")

# 8. Change results to a GRanges file
act_diff_CpG_gr =  as(act_diff_CpG_data,"GRanges")
# !! CHECKPOINT !! - save methylation data as a GRanges file
save(act_diff_CpG_gr, file = "act_diff_CpG_gr_exp2.RData")

# 9. Load in the annotated genome file
methylation_matching <- read.table("methylation_matching_file.tsv", header = FALSE, sep = "\t", col.names = c("chr", "utg_name"))

# 10. Rename columns for clarity
colnames(act_diff_CpG_data)[1] <- "chr"
colnames(methylation_matching) <- c("chr", "utg_name")

# 11. Merge the data frames by the "chr" column
merged_data <- merge(act_diff_CpG_data, methylation_matching, by = "chr", all.x = TRUE)

# 12. Perform a left join to merge data frames
merged_data <- left_join(act_diff_CpG_data, methylation_matching, by = "chr")

# 13. Rename 'utg_name' to 'seqnames' and 'chr' to 'WHPX_names', then rearrange the columns
act_diff_CpG_data_matched <- merged_data %>%
  rename(chr = utg_name, WHPX_names = chr) %>%
  select(chr, everything(), WHPX_names)

# 14. remove the WHPX names
act_diff_CpG_data_match <- subset(act_diff_CpG_data_matched[,-2])

# 15. Transform to GRanges format
act_diff_CpG_gr <- as(act_diff_CpG_data_matched, "GRanges")

# 16. load in annotated genome
genome <- read_gff3("../Data/A_Equina_with_promoters.gff3")
# 17. Change to dataframe object
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
# 18. convert back to GRanges format
annotated_genome_gr <- as(genome, "GRanges")

# 19. Join_overlap_left_directed() on the methylation and genome files
act_diff_CpG_gr_ann = join_overlap_left_directed(act_diff_CpG_gr, annotated_genome_gr)

# 20. change to df format
act_diff_CpG_df_ann = as(act_diff_CpG_gr_ann, "data.frame")

# 21. Calculate p-adjusted values (bonferonni-hochberg)
act_diff_CpG_df_ann$p_fdr = p.adjust(act_diff_CpG_df_ann$pvalue, method = "fdr")

# Correct variable format
act_diff_CpG_df_ann$Parent <- as.character(act_diff_CpG_df_ann$Parent)

# 22. Write table of results - !! CHECKPOINT !!
write.table(act_diff_CpG_df_ann, "Exp2 DMBs in CpG context.txt", sep="\t", row.names=FALSE)

# 23. Volcano plot of differential methylation

data <- read.table("Exp2 DMBs in CpG context.txt", sep = "\t", header = T)

data$diffmeth <- 'NO'
data$diffmeth[data$meth.diff >= (5) & data$p_fdr < 0.1] <- 'UP'
data$diffmeth[data$meth.diff <= (-5) & data$p_fdr < 0.1] <- 'DOWN'
data$diffmeth <- as.factor(data$diffmeth)

volcano_plot <- ggplot(data = data,
       aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) + 
  geom_point() + 
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") + 
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") + 
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_color_manual(values = c("deepskyblue","gray", "brown1"),
                     labels = c("Down Methylated", "Not Significant", "Up Methylated")) + 
  labs(
    x = "Differential methyltion %", y = expression("-log"[10]*"p-adj")
  ) + 
  ggtitle("Differential Methylation of CpGs of Experiment 2")

ggsave("volcano_plot_exp2_CpG.png", volcano_plot, width = 18, height = 12, units = "cm")

	  #  ---------- CHG ------------ #


# Input list of files for each sample to use:
file.list = list("../Data/Trimmed/Sample_6-6_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_7-7_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_8-8_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_10-10_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_14-14_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_20-20_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_21-21_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_29-11_redo_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_32-18_redo_D/meth_CHG_cov_reads",
		"../Data/Trimmed/Sample_36-2_redo_D/meth_CHG_cov_reads")

# Write methRead() function with selected options for Experimental Design                  
act_CHG = methRead(file.list,
                sample.id = list("6","7","8","10","14","20","21","11","18","2"),
                assembly="actinia",
                pipeline = "bismarkCoverage",
                treatment = c(0,0,1,0,1,1,1,0,0,1),
                context = "CHG",
                dbtype="tabix")

# Print the Raw methylation plot and save - PLOT
png("raw_methylation_plots_act_CHG_experiment_2.png")
lapply(act_CHG, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save - PLOT
png("raw_coverage_plots_act_CHG_experiment_2.png")
lapply(act_CHG, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# 3. Filter the methylation results by a coverage threshold
act_CHG_f = filterByCoverage(act_CHG, lo.count = 10, hi.count = 32, suffix = "CHG_f")
				      # Change lo.count & hi.count based on the methylation & coverage plots

# 4. Unite the methylation results
act_CHG_fu = methylKit::unite(act_CHG_f, destrand=FALSE, min.per.group = 3L, suffix = "CHG_fu3")

# Plot clustering of filtered data
png("cluster_plot_act_CHG_experiment_2.png")
clusterSamples(act_CHG_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

# Plot PCAs
png("PCA_1-2_plot_act_CHG_multi_experiment_2.png")
PCASamples(act_CHG_fu, comp =c(1,2))
dev.off()

png("PCA_3-4_plot_act_CHG_multi_experiment_2.png")
PCASamples(act_CHG_fu, comp =c(3,4))
dev.off()

png("PCA_5-6_plot_act_CHG_multi_experiment_2.png")
PCASamples(act_CHG_fu, comp =c(5,6))
dev.off()

png("PCA_7-8_plot_act_CHG_multi_experiment_2.png")
PCASamples(act_CHG_fu, comp =c(7,8))
dev.off()

# Hierarchical clustering of samples
# methMatrix <- getMethylationStats(act_CHG_fu, type="base")
# distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
# hc <- hclust(distMatrix)
# plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))

# !! Currently not working, states:
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function �getMethylationStats� for signature �"methylBaseDB"�

# 5. Calculate Differential Methylation ----
dmb_CHG = calculateDiffMeth(act_CHG_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CHGfu3odMNtestC")

# Scatter plot of differential methylation
diffMethData <- getMethylDiff(dmb_CHG, qvalue=0.99)

png("scatter_plot_methylation_exp1_CHG.png")
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")
dev.off()

# !! Error !!
# Error in as.double(y) :
# cannot coerce type 'S4' to vector of type 'double'

# 6. Extract Methylation Difference Results ----
act_diff_CHG = getMethylDiff(dmb_CHG, qvalue=0.99, type="all")
# having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

# !! CHECKPOINT !! save work to RData file
save(act_diff_CHG, file = "Actinia_exp1_DMBs_d5_CHG_experiment_2.RData")

# 7. Get Data from results ----
act_diff_CHG_data = getData(act_diff_CHG)
# !! CHECKPOINT !! - save methylation data
save(act_diff_CHG_data, file = "act_diff_CHG_data_exp2.RData")

# 8. Change results to a GRanges file
act_diff_CHG_gr =  as(act_diff_CHG_data,"GRanges")
# !! CHECKPOINT !! - save methylation data as a GRanges file
save(act_diff_CHG_gr, file = "act_diff_CHG_gr_exp2.RData")

# 9. Load in the annotated genome file
methylation_matching <- read.table("methylation_matching_file.tsv", header = FALSE, sep = "\t", col.names = c("chr", "utg_name"))

# 10. Rename columns for clarity
colnames(act_diff_CHG_data)[1] <- "chr"
colnames(methylation_matching) <- c("chr", "utg_name")

# 11. Merge the data frames by the "chr" column
merged_data <- merge(act_diff_CHG_data, methylation_matching, by = "chr", all.x = TRUE)

# 12. Perform a left join to merge data frames
merged_data <- left_join(act_diff_CHG_data, methylation_matching, by = "chr")

# 13. Rename 'utg_name' to 'seqnames' and 'chr' to 'WHPX_names', then rearrange the columns
act_diff_CHG_data_matched <- merged_data %>%
  rename(chr = utg_name, WHPX_names = chr) %>%
  select(chr, everything(), WHPX_names)

# 14. remove the WHPX names
act_diff_CHG_data_match <- subset(act_diff_CHG_data_matched[,-2])

# 15. Transform to GRanges format
act_diff_CHG_gr <- as(act_diff_CHG_data_matched, "GRanges")

# 16. load in annotated genome
genome <- read_gff3("../Data/A_Equina_with_promoters.gff3")
# 17. Change to dataframe object
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
# 18. convert back to GRanges format
annotated_genome_gr <- as(genome, "GRanges")

# 19. Join_overlap_left_directed() on the methylation and genome files
act_diff_CHG_gr_ann = join_overlap_left_directed(act_diff_CHG_gr, annotated_genome_gr)

# 20. change to df format
act_diff_CHG_df_ann = as(act_diff_CHG_gr_ann, "data.frame")

# 21. Calculate p-adjusted values (bonferonni-hochberg)
act_diff_CHG_df_ann$p_fdr = p.adjust(act_diff_CHG_df_ann$pvalue, method = "fdr")

# Correct variable format
act_diff_CHG_df_ann$Parent <- as.character(act_diff_CHG_df_ann$Parent)

# 22. Write table of results - !! CHECKPOINT !!
write.table(act_diff_CHG_df_ann, "Exp2 DMBs in CHG context.txt", sep="\t", row.names=FALSE)

# 23. Volcano plot of differential methylation

data <- read.table("Exp2 DMBs in CHG context.txt", sep = "\t", header = T)

data$diffmeth <- 'NO'
data$diffmeth[data$meth.diff >= (5) & data$p_fdr < 0.1] <- 'UP'
data$diffmeth[data$meth.diff <= (-5) & data$p_fdr < 0.1] <- 'DOWN'
data$diffmeth <- as.factor(data$diffmeth)

volcano_plot <- ggplot(data = data,
       aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) + 
  geom_point() + 
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") + 
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") + 
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_color_manual(values = c("deepskyblue","gray", "brown1"),
                     labels = c("Down Methylated", "Not Significant", "Up Methylated")) + 
  labs(
    x = "Differential methyltion %", y = expression("-log"[10]*"p-adj")
  ) + 
  ggtitle("Differential Methylation of CHGs of Experiment 2")

ggsave("volcano_plot_exp2_CHG.png", volcano_plot, width = 18, height = 12, units = "cm")

	  #  ---------- CHH ------------ #


# Input list of files for each sample to use:
file.list = list("../Data/Trimmed/Sample_6-6_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_7-7_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_8-8_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_10-10_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_14-14_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_20-20_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_21-21_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_29-11_redo_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_32-18_redo_D/meth_CHH_cov_reads",
		"../Data/Trimmed/Sample_36-2_redo_D/meth_CHH_cov_reads")

# Write methRead() function with selected options for Experimental Design                  
act_CHH = methRead(file.list,
                sample.id = list("6","7","8","10","14","20","21","11","18","2"),
                assembly="actinia",
                pipeline = "bismarkCoverage",
                treatment = c(0,0,1,0,1,1,1,0,0,1),
                context = "CHH",
                dbtype="tabix")

# Print the Raw methylation plot and save - PLOT
png("raw_methylation_plots_act_CHH_experiment_2.png")
lapply(act_CHH, getMethylationStats,  plot=T)
dev.off()

# Print the Raw Coverage plot and save - PLOT
png("raw_coverage_plots_act_CHH_experiment_2.png")
lapply(act_CHH, getCoverageStats,  plot=T, xlim=c(log10(5),2))
dev.off()

# 3. Filter the methylation results by a coverage threshold
act_CHH_f = filterByCoverage(act_CHH, lo.count = 10, hi.count = 32, suffix = "CHH_f")
				      # Change lo.count & hi.count based on the methylation & coverage plots

# 4. Unite the methylation results
act_CHH_fu = methylKit::unite(act_CHH_f, destrand=FALSE, min.per.group = 3L, suffix = "CHH_fu3")

# Plot clustering of filtered data
png("cluster_plot_act_CHH_experiment_2.png")
clusterSamples(act_CHH_fu,
               filterByQuantile = F,
               sd.threshold = 0.1,
               dist="correlation",
               method="ward.D",
               plot=T)
dev.off()

# Plot PCAs
png("PCA_1-2_plot_act_CHH_multi_experiment_2.png")
PCASamples(act_CHH_fu, comp =c(1,2))
dev.off()

png("PCA_3-4_plot_act_CHH_multi_experiment_2.png")
PCASamples(act_CHH_fu, comp =c(3,4))
dev.off()

png("PCA_5-6_plot_act_CHH_multi_experiment_2.png")
PCASamples(act_CHH_fu, comp =c(5,6))
dev.off()

png("PCA_7-8_plot_act_CHH_multi_experiment_2.png")
PCASamples(act_CHH_fu, comp =c(7,8))
dev.off()

# Hierarchical clustering of samples
# methMatrix <- getMethylationStats(act_CHH_fu, type="base")
# distMatrix <- dist(t(methMatrix))  # transpose to cluster by sample
# hc <- hclust(distMatrix)
# plot(hc)

# Using heatmap to visualize clustering and methylation levels
heatmap(as.matrix(methMatrix), Colv = NA, scale="row", margins=c(10, 10))

# !! Currently not working, states:
# Error in (function (classes, fdef, mtable)  :
# unable to find an inherited method for function �getMethylationStats� for signature �"methylBaseDB"�

# 5. Calculate Differential Methylation ----
dmb_CHH = calculateDiffMeth(act_CHH_fu,
                            covariates = NULL,
                            overdispersion = "MN",
                            #test = "F",
                            mc.cores = 32,
                            suffix = "CHHfu3odMNtestC")

# Scatter plot of differential methylation
diffMethData <- getMethylDiff(dmb_CHH, qvalue=0.99)

png("scatter_plot_methylation_exp1_CHH.png")
plot(diffMethData, xlab="Methylation Difference (%)", ylab="Adjusted P-Value", main="Differential Methylation Scatterplot")
dev.off()

# !! Error !!
# Error in as.double(y) :
# cannot coerce type 'S4' to vector of type 'double'

# 6. Extract Methylation Difference Results ----
act_diff_CHH = getMethylDiff(dmb_CHH, qvalue=0.99, type="all")
# having diff = 5 reduces size of object which allows us to handle it, otherwise the system always falls over

# !! CHECKPOINT !! save work to RData file
save(act_diff_CHH, file = "Actinia_exp1_DMBs_d5_CHH_experiment_2.RData")

# 7. Get Data from results ----
act_diff_CHH_data = getData(act_diff_CHH)
# !! CHECKPOINT !! - save methylation data
save(act_diff_CHH_data, file = "act_diff_CHH_data_exp2.RData")

# 8. Change results to a GRanges file
act_diff_CHH_gr =  as(act_diff_CHH_data,"GRanges")
# !! CHECKPOINT !! - save methylation data as a GRanges file
save(act_diff_CHH_gr, file = "act_diff_CHH_gr_exp2.RData")

# 9. Load in the annotated genome file
methylation_matching <- read.table("methylation_matching_file.tsv", header = FALSE, sep = "\t", col.names = c("chr", "utg_name"))

# 10. Rename columns for clarity
colnames(act_diff_CHH_data)[1] <- "chr"
colnames(methylation_matching) <- c("chr", "utg_name")

# 11. Merge the data frames by the "chr" column
merged_data <- merge(act_diff_CHH_data, methylation_matching, by = "chr", all.x = TRUE)

# 12. Perform a left join to merge data frames
merged_data <- left_join(act_diff_CHH_data, methylation_matching, by = "chr")

# 13. Rename 'utg_name' to 'seqnames' and 'chr' to 'WHPX_names', then rearrange the columns
act_diff_CHH_data_matched <- merged_data %>%
  rename(chr = utg_name, WHPX_names = chr) %>%
  select(chr, everything(), WHPX_names)

# 14. remove the WHPX names
act_diff_CHH_data_match <- subset(act_diff_CHH_data_matched[,-2])

# 15. Transform to GRanges format
act_diff_CHH_gr <- as(act_diff_CHH_data_matched, "GRanges")

# 16. load in annotated genome
genome <- read_gff3("../Data/A_Equina_with_promoters.gff3")
# 17. Change to dataframe object
genome <- as.data.frame(genome)
# replace NA data with "NA"
genome[] <- lapply(genome, function(x) {
    # Check if the column is a factor and handle it accordingly
    if(is.factor(x)) {
        # Convert factor to character to avoid levels issues
        x <- as.character(x)
    }
    # Replace NA with "NA" (as character string)
    x[is.na(x)] <- "NA"
    return(x)
})
# 18. convert back to GRanges format
annotated_genome_gr <- as(genome, "GRanges")

# 19. Join_overlap_left_directed() on the methylation and genome files
act_diff_CHH_gr_ann = join_overlap_left_directed(act_diff_CHH_gr, annotated_genome_gr)

# 20. change to df format
act_diff_CHH_df_ann = as(act_diff_CHH_gr_ann, "data.frame")

# 21. Calculate p-adjusted values (bonferonni-hochberg)
act_diff_CHH_df_ann$p_fdr = p.adjust(act_diff_CHH_df_ann$pvalue, method = "fdr")

# Correct variable format
act_diff_CHH_df_ann$Parent <- as.character(act_diff_CHH_df_ann$Parent)

# 22. Write table of results - !! CHECKPOINT !!
write.table(act_diff_CHH_df_ann, "Exp2 DMBs in CHH context.txt", sep="\t", row.names=FALSE)

# 23. Volcano plot of differential methylation

data <- read.table("Exp2 DMBs in CHH context.txt", sep = "\t", header = T)

data$diffmeth <- 'NO'
data$diffmeth[data$meth.diff >= (5) & data$p_fdr < 0.1] <- 'UP'
data$diffmeth[data$meth.diff <= (-5) & data$p_fdr < 0.1] <- 'DOWN'
data$diffmeth <- as.factor(data$diffmeth)

volcano_plot <- ggplot(data = data,
       aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) + 
  geom_point() + 
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") + 
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") + 
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") + 
  theme_classic() +
  theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = "black"),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_color_manual(values = c("deepskyblue","gray", "brown1"),
                     labels = c("Down Methylated", "Not Significant", "Up Methylated")) + 
  labs(
    x = "Differential methyltion %", y = expression("-log"[10]*"p-adj")
  ) + 
  ggtitle("Differential Methylation of CHHs of Experiment 2")

ggsave("volcano_plot_exp2_CHH.png", volcano_plot, width = 18, height = 12, units = "cm")
