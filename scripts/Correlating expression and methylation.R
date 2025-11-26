# ---- Experiment 1 ----

# Load packages
library(ggplot2)

# Install necessary packages
if (!requireNamespace("energy", quietly = TRUE)) {
  install.packages("energy")
}
if (!requireNamespace("Hmisc", quietly = TRUE)) {
  install.packages("Hmisc")
}
if (!requireNamespace("minerva", quietly = TRUE)) {
  install.packages("minerva")
}

# Load necessary packages
library(energy)
library(Hmisc)
library(minerva)


# Load in methylation data
CpG_exp1 <- read.table("../Exp1 DMBs in CpG context.txt", sep = "\t", header = T)
CHG_exp1 <- read.table("../Exp1 DMBs in CHG context.txt", sep = "\t", header = T)
CHH_exp1 <- read.table("../Exp1 DMBs in CHH context.txt", sep = "\t", header = T)

# load in expression data

data <- read.csv("../DEresults_exp1_outlier.csv")
expression_data <- subset(data[,c("X","log2FoldChange","padj")])
colnames(expression_data) <- c("ID","LFC","expression.p.adj")

# CpG exp1
meth_data <- ""
data <- subset(CpG_exp1[,c("meth.diff","ID","p_fdr")])
colnames(data) <- c("CpG.meth.diff","ID","CpG.p.fdr")
meth_data <- data

combined_data <- ""
combined_data <- merge(meth_data, expression_data, by = "ID")
plot(combined_data$LFC, combined_data$CpG.meth.diff)

CpG_plot_exp1 <- ggplot(data = combined_data, aes(x = LFC, y = CpG.meth.diff)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Corrlation between Expression and Methylation changes in CpGs - Experiment 1") + 
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(limits = c(-50,50), breaks = seq(-50,50,5)) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black")

print(CpG_plot_exp1)
ggsave(CpG_plot_exp1, file = "CpG_Expression_Exp1_Corr_plot.png")

shapiro.test(combined_data$CpG.meth.diff)
shapiro.test(combined_data$LFC)

# compute distance correlation
dcor.test(combined_data$LFC, combined_data$CpG.meth.diff)
# compute Hoeffdings D statistic
hoeffd(combined_data$LFC, combined_data$CpG.meth.diff)
# compute maximal information coeffidient
mine(combined_data$LFC, combined_data$CpG.meth.diff)

## looks like we cant calculate an accurate p value for correlation using these alternative corr
## methods (alternative because neither variable conforms to normal distribution, monotonic or 
## a linear relationship (although one might expect one))

## Might want to consider keeping the raw data without any threshold restrictions
## to refine the dataset into a smaller and more tnagible form to interpret results from
## we are missing a whole range of values within the > -5% & < 5% methylation difference

## Also thinking back to when I used the lo.count and hi.count thresholds to restrict the data
## for the methylation results, I could probably remove that, since we already have trimmed
## the data based on a methylation bias plot, and then re-analyse the results to see if 
## we can retain some of the methylated cytosine areas which can be correlated to the 
## expression data

## Another point to add here is that there is a significantly smaller number of 
## observations in the combined (merged) dataset than there are in either the 
## methylation data or the expression data. There seems to be a large number of areas
## of methylated genes which are not related to named genes, and therefore cannot be
## matched and merged together appropriatley, leading to a list of 2608 (CpG_exp1)
## cytosines and 46,331 aligned expressed genes being reduced to a list of only
## 486 points (~18% of methylated points, and ~1% of the aligned expression genes)

  ## We expect there to be much lower number of genes from the expression file
## becuase methylation does not act on RNA and transcripts all the time, and can
## act on areas around the coding regions of genes (and therefore not attributed
## to a specific gene name).

CpG_exp1$type <- as.factor(CpG_exp1$type)
levels(CpG_exp1$type) <- c(levels(CpG_exp1$type), "Missing")
CpG_exp1$type[is.na(CpG_exp1$type)] <- "Missing"
summary(CpG_exp1$type)
plot(CpG_exp1$type)

## PLEASE NOTE: the same areas are often labelled multiple times, therefore creating
## inflated number of data points, especially within coding regions of genomic sequences
  
# CHG exp1
meth_data <- ""
data <- subset(CHG_exp1[,c("meth.diff","ID","p_fdr")])
colnames(data) <- c("CHG.meth.diff","ID","CHG.p.fdr")
meth_data <- data

combined_data <- ""
combined_data <- merge(meth_data, expression_data, by = "ID")
plot(combined_data$LFC, combined_data$CHG.meth.diff)

CHG_plot_exp1 <- ggplot(data = combined_data, aes(x = LFC, y = CHG.meth.diff)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Corrlation between Expression and Methylation changes in CHGs - Experiment 1") + 
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(limits = c(-50,50), breaks = seq(-50,50,5)) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black")

print(CHG_plot_exp1)
ggsave(CHG_plot_exp1, file = "CHG_Expression_Exp1_Corr_plot.png")

CHG_exp1$type <- as.factor(CHG_exp1$type)
levels(CHG_exp1$type) <- c(levels(CHG_exp1$type), "Missing")
CHG_exp1$type[is.na(CHG_exp1$type)] <- "Missing"
summary(CHG_exp1$type)
plot(CHG_exp1$type, xlab = "type", ylab = "Count", main = "CHG Experiment 1 Types")


# CHH exp1
data <- ""
meth_data <- ""
data <- subset(CHH_exp1[,c("meth.diff","ID","p_fdr")])
colnames(data) <- c("CHH.meth.diff","ID","CHH.p.fdr")
meth_data <- data

combined_data <- ""
combined_data <- merge(meth_data, expression_data, by = "ID")
plot(combined_data$LFC, combined_data$CHH.meth.diff)

CHH_plot_exp1 <- ggplot(data = combined_data, aes(x = LFC, y = CHH.meth.diff)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Corrlation between Expression and Methylation changes in CHHs - Experiment 1") + 
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(limits = c(-50,50), breaks = seq(-50,50,5)) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black")

print(CHH_plot_exp1)
ggsave(CHH_plot_exp1, file = "CHH_Expression_Exp1_Corr_plot.png")

CHH_exp1$type <- as.factor(CHH_exp1$type)
levels(CHH_exp1$type) <- c(levels(CHH_exp1$type), "Missing")
CHH_exp1$type[is.na(CHH_exp1$type)] <- "Missing"
summary(CHH_exp1$type)
plot(CHH_exp1$type, xlab = "type", ylab = "Count", main = "CHH Experiment 1 Types")

# ---- Experiment 2 ----

# Load in data
CpG_exp2 <- read.table("../Exp2 DMBs in CpG context.txt", sep = "\t", header = T)
CHG_exp2 <- read.table("../Exp2 DMBs in CHG context.txt", sep = "\t", header = T)
CHH_exp2 <- read.table("../Exp2 DMBs in CHH context.txt", sep = "\t", header = T)

data <- read.csv("../DEresults_exp2.csv")
expression_data <- subset(data[,c("X","log2FoldChange","padj")])
colnames(expression_data) <- c("ID","LFC","expression.p.adj")

# CpG exp2
meth_data <- ""
data <- ""
data <- subset(CpG_exp2[,c("meth.diff","ID","p_fdr")])
colnames(data) <- c("meth.diff","ID","p.fdr")
meth_data <- data

combined_data <- ""
combined_data <- merge(meth_data, expression_data, by = "ID")
combined_data$ID <- as.factor(combined_data$ID)
plot(combined_data$LFC, combined_data$meth.diff)

CpG_plot_exp2 <- ggplot(data = combined_data, aes(x = LFC, y = meth.diff)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Corrlation between Expression and Methylation changes in CpGs - Experiment 2") + 
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(limits = c(-50,50), breaks = seq(-50,50,5)) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black")

print(CpG_plot_exp2)
ggsave(CpG_plot_exp2, file = "CpG_Expression_Exp2_Corr_plot.png")

combined_data <- combined_data[order(combined_data$meth.diff, decreasing = T),]
plot(combined_data$ID, combined_data$meth.diff, col = "grey",
     ylab = "methylation difference (%)", xlab = "unique genes",
     main = "methylation difference % across and within genes")

CpG_exp2$type <- as.factor(CpG_exp2$type)
levels(CpG_exp2$type) <- c(levels(CpG_exp2$type), "Missing")
CpG_exp2$type[is.na(CpG_exp2$type)] <- "Missing"
summary(CpG_exp2$type)
plot(CpG_exp2$type, xlab = "type", ylab = "Count", main = "CpG Experiment 2 Types")

# CHG exp 2
meth_data <- ""
data <- ""
data <- subset(CHG_exp2[,c("meth.diff","ID","p_fdr")])
colnames(data) <- c("meth.diff","ID","p.fdr")
meth_data <- data

combined_data <- ""
combined_data <- merge(meth_data, expression_data, by = "ID")
plot(combined_data$LFC, combined_data$meth.diff)

CHG_plot_exp2 <- ggplot(data = combined_data, aes(x = LFC, y = meth.diff)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Corrlation between Expression and Methylation changes in CHGs - Experiment 2") + 
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(limits = c(-50,50), breaks = seq(-50,50,5)) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black")

print(CHG_plot_exp2)
ggsave(CHG_plot_exp2, file = "CHG_Expression_Exp2_Corr_plot.png")

CHG_exp2$type <- as.factor(CHG_exp2$type)
levels(CHG_exp2$type) <- c(levels(CHG_exp2$type), "Missing")
CHG_exp2$type[is.na(CHG_exp2$type)] <- "Missing"
summary(CHG_exp2$type)
plot(CHG_exp2$type, xlab = "type", ylab = "Count", main = "CHG Experiment 2 Types")

# CHH exp 2
meth_data <- ""
data <- ""
data <- subset(CHH_exp2[,c("meth.diff","ID","p_fdr")])
colnames(data) <- c("meth.diff","ID","p.fdr")
meth_data <- data

combined_data <- ""
combined_data <- merge(meth_data, expression_data, by = "ID")
plot(combined_data$LFC, combined_data$meth.diff)

CHH_plot_exp2 <- ggplot(data = combined_data, aes(x = LFC, y = meth.diff)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Corrlation between Expression and Methylation changes in CHHs - Experiment 2") + 
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(limits = c(-50,50), breaks = seq(-50,50,5)) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black")

print(CHH_plot_exp2)
ggsave(CHH_plot_exp2, file = "CHH_Expression_Exp2_Corr_plot.png")

CHH_exp2$type <- as.factor(CHH_exp2$type)
levels(CHH_exp2$type) <- c(levels(CHH_exp2$type), "Missing")
CHH_exp2$type[is.na(CHH_exp2$type)] <- "Missing"
summary(CHH_exp2$type)
plot(CHH_exp2$type, xlab = "type", ylab = "Count", main = "CHH Experiment 2 Types")

par(mfrow = c(3,2))
plot(CpG_exp1$type, xlab = "type", ylab = "Count", main = "CpG Experiment 1 Types")
plot(CpG_exp2$type, xlab = "type", ylab = "Count", main = "CpG Experiment 2 Types")
plot(CHG_exp1$type, xlab = "type", ylab = "Count", main = "CHG Experiment 1 Types")
plot(CHG_exp2$type, xlab = "type", ylab = "Count", main = "CHG Experiment 2 Types")
plot(CHH_exp1$type, xlab = "type", ylab = "Count", main = "CHH Experiment 1 Types")
plot(CHH_exp2$type, xlab = "type", ylab = "Count", main = "CHH Experiment 2 Types")
par(mfrow = c(1,1))
