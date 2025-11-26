# Load packages
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)
library(tidyr)

# Import Data
getwd()
setwd("C:/Users/r02hw22/OneDrive/Documents/PhD/Chapter 2 - Differential Methylation Analyses/Methylation_Analyses")
meth_cov <- read.csv("Data/Bismark Results/Bismark_Results_All.csv")

# Convert percentage to numeric
meth_cov$CPG_meth <- as.numeric(sub("%", "", meth_cov$CPG_meth))
meth_cov$CHG_meth <- as.numeric(sub("%", "", meth_cov$CHG_meth))
meth_cov$CHH_meth <- as.numeric(sub("%", "", meth_cov$CHH_meth))

# Reshape the data to long format
meth_cov_long <- melt(meth_cov, id.vars = c("Sample", "Group", "Clone", "TC_Analysed", "TM_CPG", "TM_CHG", "TM_CHH", "C2T_CPG", "C2T_CHG", "C2T_CHH"), measure.vars = c("CPG_meth", "CHG_meth", "CHH_meth"))

# Rename the columns to more descriptive names
colnames(meth_cov_long) <- c("Sample", "Group", "Clone", "TC_Analysed", "TM_CPG", "TM_CHG", "TM_CHH", "C2T_CPG", "C2T_CHG", "C2T_CHH", "Methylation_Type", "Percentage")

# Modify clone format
meth_cov_long$Clone <- as.factor(meth_cov_long$Clone)
meth_cov$Clone <- as.factor(meth_cov$Clone)

# Factor the 'Group' variable to specify the order of the groups
meth_cov_long$Group <- factor(meth_cov_long$Group, levels = c("X", "A", "B", "C"))
meth_cov$Group <- factor(meth_cov$Group, levels = c("X", "A", "B", "C"))

# Plot the boxplots on the same plot


ggplot(meth_cov_long, aes(x = Methylation_Type, y = Percentage, fill = Group)) + 
  geom_boxplot() + 
  facet_wrap(~Methylation_Type, scales = "free") +
  labs(title = "Methylation Percentages by Type and Group") + 
  theme_classic()

meth_exp1 <- subset(meth_cov[meth_cov$Group == "A" | meth_cov$Group == "X",])
meth_exp2 <- subset(meth_cov[meth_cov$Group == "B" | meth_cov$Group == "C",])

t.test(meth_exp1$CPG_meth~meth_exp1$Group)
t.test(meth_exp2$CPG_meth~meth_exp2$Group)

t.test(meth_exp1$CHG_meth~meth_exp1$Group)
t.test(meth_exp2$CHG_meth~meth_exp2$Group)

t.test(meth_exp1$CHH_meth~meth_exp1$Group)
t.test(meth_exp2$CHH_meth~meth_exp2$Group)

# removing the outlier clone pairs
meth_cov_long_rm_outliers <- subset(meth_cov_long, !Sample %in% c("3", "17", "8", "18"))

ggplot(meth_cov_long_rm_outliers, aes(x = Methylation_Type, y = Percentage, fill = Group)) + 
  geom_boxplot() + 
  facet_wrap(~Methylation_Type, scales = "free") +
  labs(title = "Methylation Percentages by Type and Group, no outliers") + 
  theme_classic()

library(gghalves)
library(ggplot2)  
library(ggside)

ggplot(meth_cov_long_rm_outliers, aes(x = Methylation_Type, y = Percentage, fill = Group)) + 
  geom_half_boxplot(side = "l") +
  geom_half_point(side = "r", aes(color = Group)) +
  facet_wrap(~Methylation_Type, scales = "free") +
  labs(title = "Methylation Percentages by Type and Group, no outliers") + 
  theme_classic()

ggplot(meth_cov_long_rm_outliers, aes(x = Methylation_Type, y = Percentage, color = Group)) + 
  geom_point() + 
  facet_wrap(~Methylation_Type, scales = "free") +
  labs(title = "Methylation Percentages by Type and Group, no outliers") + 
  theme_classic()

meth_cov_long_rm_outliers <- as.data.frame(meth_cov_long_rm_outliers)

#-----------------------------------------------------------------------------------------------------------

# 2D directional plot ----
meth_cov <- read.csv("./Data/Bismark Results/Bismark_Results_All.csv")
meth_cov$Clone <- as.factor(meth_cov$Clone)

# Convert percentage to numeric
meth_cov$CPG_meth <- as.numeric(sub("%", "", meth_cov$CPG_meth))
meth_cov$CHG_meth <- as.numeric(sub("%", "", meth_cov$CHG_meth))
meth_cov$CHH_meth <- as.numeric(sub("%", "", meth_cov$CHH_meth))

# Reshape data for plotting
df <- meth_cov[,c(1:3,11:13)]
df$experiment <- c(1,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,1,1)

# re-order
df$Clone <- factor(df$Clone, levels = c("2","5","9","15","6","8","12","19","20"))

# Split data into group A and X, and B & C
data_group_A <- df %>% filter(Group == 'A')
data_group_X <- df %>% filter(Group == 'X')
data_group_B <- df %>% filter(Group == 'B')
data_group_C <- df %>% filter(Group == 'C')

# Join on Clone to get start and end points for arrows
arrow_data_1 <- inner_join(data_group_X, data_group_A, by = "Clone", suffix = c("_X", "_A"))
arrow_data_2 <- inner_join(data_group_C, data_group_B, by = "Clone", suffix = c("_C", "_B"))

# CPG Methylation percentages
ggplot(data = df, aes(x = Clone, y = CPG_meth)) +
  geom_point(aes(color = Group), size = 8) +
  geom_segment(data = arrow_data_1, aes(x = Clone, 
                                      y = CPG_meth_X, 
                                      xend = Clone, 
                                      yend = CPG_meth_A,), 
               arrow = arrow(type = "closed", 
                             length = unit(0.20, "inches")), 
               color = "black") +
  geom_segment(data = arrow_data_2, aes(x = Clone, 
                                      y = CPG_meth_B, 
                                      xend = Clone, 
                                      yend = CPG_meth_C,), 
               arrow = arrow(type = "closed", 
                             length = unit(0.20, "inches")), 
               color = "black") +
  scale_color_manual(values = c("A" = "red", "X" = "blue", "B" = "green3", "C" = "goldenrod2")) +
  scale_fill_manual(values = c("A" = "red", "X" = "blue", "B" = "green3", "C" = "goldenrod2")) +
  theme_classic() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(
    x = paste("Clone pair ID"),
    y = paste("CpG Methylation %"), 
    title = "CpG Methylation Changes in Experiment") # + facet_wrap(~experiment)

# Looks to me like there is quite an inconsistent change of metylation in experiemnt 1,
# but in expriment 2, there is a much more consistent change in methylation.

# CHG methylation percentages
ggplot(data = df, aes(x = Clone, y = CHG_meth)) +
  geom_point(aes(color = Group), size = 8) +
  geom_segment(data = arrow_data_1, aes(x = Clone, 
                                        y = CHG_meth_X, 
                                        xend = Clone, 
                                        yend = CHG_meth_A,), 
               arrow = arrow(type = "closed", 
                             length = unit(0.20, "inches")), 
               color = "black") +
  geom_segment(data = arrow_data_2, aes(x = Clone, 
                                        y = CHG_meth_B, 
                                        xend = Clone, 
                                        yend = CHG_meth_C,), 
               arrow = arrow(type = "closed", 
                             length = unit(0.20, "inches")), 
               color = "black") +
  scale_color_manual(values = c("A" = "red", "X" = "blue", "B" = "green3", "C" = "goldenrod2")) +
  scale_fill_manual(values = c("A" = "red", "X" = "blue", "B" = "green3", "C" = "goldenrod2")) +
  theme_classic() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(
    x = paste("Clone pair ID"),
    y = paste("CHG Methylation %"), 
    title = "CHG Methylation Changes in Experiment") # + facet_wrap(~experiment)

# No obvious patterns, but more changes observed in experiment 2

# CHH methylation percentages
ggplot(data = df, aes(x = Clone, y = CHH_meth)) +
  geom_point(aes(color = Group), size = 8) +
  geom_segment(data = arrow_data_1, aes(x = Clone, 
                                        y = CHH_meth_X, 
                                        xend = Clone, 
                                        yend = CHH_meth_A,), 
               arrow = arrow(type = "closed", 
                             length = unit(0.20, "inches")), 
               color = "black") +
  geom_segment(data = arrow_data_2, aes(x = Clone, 
                                        y = CHH_meth_B, 
                                        xend = Clone, 
                                        yend = CHH_meth_C,), 
               arrow = arrow(type = "closed", 
                             length = unit(0.20, "inches")), 
               color = "black") +
  scale_color_manual(values = c("A" = "red", "X" = "blue", "B" = "green3", "C" = "goldenrod2")) +
  scale_fill_manual(values = c("A" = "red", "X" = "blue", "B" = "green3", "C" = "goldenrod2")) +
  theme_classic() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(
    x = paste("Clone pair ID"),
    y = paste("CHH Methylation %"), 
    title = "CHH Methylation Changes in Experiment") # + facet_wrap(~experiment)

# More conistent changes to an increase of methylation, but clone pair 20 and 9 are against this majority trend



# -----------------------------------------------------------------------------------------------------------
