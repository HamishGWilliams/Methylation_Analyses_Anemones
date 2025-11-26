# Packages
library(BiocManager)
library(methylKit)
library(plyranges)
library(GenomicRanges)
library(qqman)
library(ggplot2)

# Load Data
getwd()
setwd("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/Methylation_Analyses/All contexts DMBs results")
getwd()

exp1_CpG <- read.table("act_diff_df_ann_exp1_CpG_removing_outliers.txt", sep = '\t', header = T)
exp1_CpG$context <- "CpG"
exp1_CHG <- read.table("act_diff_df_ann_exp1_CHG_removing_outliers.txt", sep = '\t', header = T)
exp1_CHG$context <- "CHG"
exp1_CHH <- read.table("act_diff_df_ann_exp1_CHH_removing_outliers.txt", sep = '\t', header = T)
exp1_CHH$context <- "CHH"

exp2_CpG <- read.table("act_diff_df_ann_exp2_CpG_removing_outliers.txt", sep = '\t', header = T)
exp2_CpG$context <- "CpG"
exp2_CHG <- read.table("act_diff_df_ann_exp2_CHG_removing_outliers.txt", sep = '\t', header = T)
exp2_CHG$context <- "CHG"
exp2_CHH <- read.table("act_diff_df_ann_exp2_CHH_removing_outliers.txt", sep = '\t', header = T)
exp2_CHH$context <- "CHH"

# Merge DFs by experiment

# exp1
exp1_data <- rbind(exp1_CpG, exp1_CHG, exp1_CHH)

# exp2
exp2_data <- rbind(exp2_CpG, exp2_CHG, exp2_CHH)

# # Create volcano plot - Exp 1----
data_exp1 <- exp1_data

data_exp1$diffmeth <- "NS"  # Default is grey

data_exp1$diffmeth[data_exp1$meth.diff < -5 & -log10(data_exp1$p_fdr) > 1] <- "Left"   # blue
data_exp1$diffmeth[data_exp1$meth.diff > 5 & -log10(data_exp1$p_fdr) > 1]  <- "Right"  # red

data_exp1$diffmeth <- factor(data_exp1$diffmeth, levels = c("NS", "Left", "Right"))


volcano_plot_exp1 <- ggplot(data = data_exp1, aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth, shape = context)) +
                geom_point(size = 3, alpha = 1) +
                geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") +
                geom_vline(xintercept = c(5), col = "red", linetype = "dashed") +
                geom_hline(yintercept = c(1), col = "black", linetype = "dashed") +
                theme_classic() +
                theme(
                  axis.title.y = element_text(
                    face = "bold",
                    margin = margin(0, 20, 0, 0),
                    size = rel(1.1),
                    color = "black"
                  ),
                  axis.title.x = element_text(
                    hjust = 0.5,
                    face = "bold",
                    margin = margin(20, 0, 0, 0),
                    size = rel(1.1),
                    color = "black"
                  ),
                  plot.title = element_text(hjust = 0.5)
                ) +
                  scale_color_manual(values = c("NS" = "grey", "Left" = "deepskyblue2", "Right" = "coral2")) +
                labs(
                  x = "Differential methylation %",
                  y = expression("-log"[10] * "p-adj")
                )  +
  guides(color = "none", shape = guide_legend(title = "Context")) +
                # ggtitle(paste("Differential Methylation of", context, "in", experiment, "removed outliers")) + 
  scale_x_continuous(limits = c(-65, 65), breaks = c(-60,-50,-40,-30,-20,-10, 0, 10, 20,30,40,50,60)) + 
                scale_y_continuous(limits = c(0, 3), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0,2.2,2.4,2.6,2.8,3.0))

volcano_plot_exp1

ggsave("Exp1_Base_level_DMA_volcano_plot.png", volcano_plot_exp1, width = 18, height = 12, units = "cm")

# exp2 plot ----

data <- exp2_data

data$diffmeth <- "NS"  # Default is grey

data$diffmeth[data$meth.diff < -5 & -log10(data$p_fdr) > 1] <- "Left"   # blue
data$diffmeth[data$meth.diff > 5 & -log10(data$p_fdr) > 1]  <- "Right"  # red

data$diffmeth <- factor(data$diffmeth, levels = c("NS", "Left", "Right"))


volcano_plot_exp2 <- ggplot(data = data, aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth, shape = context)) +
  geom_point(size = 3, alpha = 1) +
  geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") +
  geom_vline(xintercept = c(5), col = "red", linetype = "dashed") +
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    axis.title.y = element_text(
      face = "bold",
      margin = margin(0, 20, 0, 0),
      size = rel(1.1),
      color = "black"
    ),
    axis.title.x = element_text(
      hjust = 0.5,
      face = "bold",
      margin = margin(20, 0, 0, 0),
      size = rel(1.1),
      color = "black"
    ),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c("NS" = "grey", "Left" = "deepskyblue2", "Right" = "coral2")) +
  labs(
    x = "Differential methylation %",
    y = expression("-log"[10] * "p-adj")
  )  +
  guides(color = "none", shape = guide_legend(title = "Context")) +
  # ggtitle(paste("Differential Methylation of", context, "in", experiment, "removed outliers")) + 
  scale_x_continuous(limits = c(-65, 65), breaks = c(-60,-50,-40,-30,-20,-10, 0, 10, 20,30,40,50,60)) + 
  scale_y_continuous(limits = c(0, 3), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0,2.2,2.4,2.6,2.8,3.0))

volcano_plot_exp2

ggsave("Exp2_Base_level_DMA_volcano_plot.png", volcano_plot_exp2, width = 18, height = 12, units = "cm")
getwd()

