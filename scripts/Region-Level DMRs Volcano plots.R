# load packages:
{
  library(dplyr) # for data manipulation  
  library(ggplot2) # for plotting  
  library(data.table) # for data table manipulation
}

# set worrking wd
getwd()
setwd("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/Methylation_Analyses/All context DMR results")
getwd()

# load data
exp1_DMRs <- read.table("exp1_DMRs.txt", header = T, sep = "\t")
exp2_DMRs <- read.table("exp2_DMRs.txt", header = T, sep = "\t")

# manipulate variables namez
str(exp1_DMRs)
exp1_DMRs$region <- as.factor(exp1_DMRs$region)
exp1_DMRs$Context <- as.factor(exp1_DMRs$context)
summary(exp1_DMRs$region)

str(exp2_DMRs)
exp2_DMRs$region <- as.factor(exp2_DMRs$region)
exp2_DMRs$Context <- as.factor(exp2_DMRs$context)
summary(exp2_DMRs$region)

# exp1 plot ----

volcano_plot_DMRs_exp1 <- ggplot(data = exp1_DMRs, 
                       aes(x = meth.diff, 
                           y = -log10(p_fdr), 
                           col = region,
                           shape = Context)) +
  geom_point(size = 2, alpha = 0.4) +
  #geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") +
  #geom_vline(xintercept = c(5), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed") +
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") +
  theme_bw() +
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
  scale_color_manual(
    name = "Region",
    values = c("#E69F00", "#56B4E9", "#009E73", "#a1d61a", "#0072B2", "#D55E00", "#CC79A7", "red"),
    labels = c("Repeat", "Downstream", "Exon", "5' UTR", "Gene","ncRNA", "Promoter", "3' UTR")
  ) +
  labs(
    x = "Differential methylation %",
    y = expression("-log"[10] * "p-adj")
  ) +
  #ggtitle(paste("DMRs of experiment 1")) + 
  #scale_x_continuous(limits = c(-40, 40), breaks = c(-40,-20,0,20,40)) + 
  scale_y_continuous(limits = c(0, 6), breaks = c(0,1, 2, 3, 4, 5,6)) + 
  facet_wrap(~region, nrow = 2, scales = "free") + 
  guides(color = "none")

volcano_plot_DMRs_exp1

ggsave("Exp1_DMRs_volcano_plot_all_contexts.png", volcano_plot_DMRs_exp1, width = 18, height = 12, units = "cm")

# exp2 plot ----

volcano_plot_DMRs_exp2 <- ggplot(data = exp2_DMRs, 
                                 aes(x = meth.diff, 
                                     y = -log10(p_fdr), 
                                     col = region,
                                     shape = Context)) +
  geom_point(size = 2, alpha = 0.4) +
  #geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") +
  #geom_vline(xintercept = c(5), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed") +
  geom_hline(yintercept = c(1), col = "black", linetype = "dashed") +
  theme_bw() +
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
  scale_color_manual(
    name = "Region",
    values = c("#E69F00", "#56B4E9", "#009E73", "#a1d61a", "#0072B2", "#D55E00", "#CC79A7", "red"),
    labels = c("Repeat", "Downstream", "Exon", "5' UTR", "Gene","ncRNA", "Promoter", "3' UTR")
  ) +
  labs(
    x = "Differential methylation %",
    y = expression("-log"[10] * "p-adj")
  ) +
  #ggtitle(paste("DMRs of experiment 1")) + 
  #scale_x_continuous(limits = c(-40, 40), breaks = c(-40,-20,0,20,40)) + 
  scale_y_continuous(limits = c(0, 7), breaks = c(0,1, 2, 3, 4, 5,6,7)) + 
  facet_wrap(~region, nrow = 2, scales = "free") + 
  guides(color = "none")

volcano_plot_DMRs_exp2

ggsave("Exp2_DMRs_volcano_plot_all_contexts.png", volcano_plot_DMRs_exp2, width = 18, height = 12, units = "cm")

## Summarising the data for the chapter results ----

# exp 1 ----

# plot for context
volcano_plot_DMRs_exp1

# only padj < 0.1
Padj_exp1_DMRs <- subset(exp1_DMRs[exp1_DMRs$p_fdr <= 0.1,])
summary(Padj_exp1_DMRs)

Padj_up_methylated_exp1_DMRs <- subset(Padj_exp1_DMRs[Padj_exp1_DMRs$meth.diff >= 0,])
summary(Padj_up_methylated_exp1_DMRs)

# exp 2 ----

# plot for context
volcano_plot_DMRs_exp2

# only padj < 0.1
Padj_exp2_DMRs <- subset(exp2_DMRs[exp2_DMRs$p_fdr <= 0.1,])
summary(Padj_exp2_DMRs)
summary(Padj_exp2_DMRs$region)

Padj_up_methylated_exp2_DMRs <- subset(Padj_exp2_DMRs[Padj_exp2_DMRs$meth.diff >= 0,])
summary(Padj_up_methylated_exp2_DMRs)


# Extracting data for table

Padj_exp1_DMRs$direction <- ifelse(Padj_exp1_DMRs$meth.diff > 0, "hyper", "hypo")

summary_table_exp1 <- Padj_exp1_DMRs %>%
  group_by(region, Context, direction) %>%
  summarise(count = n()) %>%
  ungroup()

library(tidyr)

summary_wide_exp1 <- summary_table_exp1 %>%
  pivot_wider(
    names_from = c(Context, direction),
    values_from = count,
    values_fill = 0
  )


# exp2 

Padj_exp2_DMRs$direction <- ifelse(Padj_exp2_DMRs$meth.diff > 0, "hyper", "hypo")

summary_table_exp2 <- Padj_exp2_DMRs %>%
  group_by(region, Context, direction) %>%
  summarise(count = n()) %>%
  ungroup()

library(tidyr)

summary_wide_exp2 <- summary_table_exp2 %>%
  pivot_wider(
    names_from = c(Context, direction),
    values_from = count,
    values_fill = 0
  )


# Combined data

DMR_table_data <- rbind(Padj_exp1_DMRs, Padj_exp2_DMRs)
DMR_table_data$experiment <- as.factor(DMR_table_data$experiment)

library(dplyr)

summary_table <- DMR_table_data %>%
  group_by(experiment, region, context, direction) %>%
  summarise(count = n(), .groups = "drop")

direction_totals <- summary_table %>%
  group_by(experiment, region, direction) %>%
  summarise(dir_total = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = dir_total, values_fill = 0) %>%
  mutate(total = hyper + hypo)

full_summary <- summary_table %>%
  left_join(direction_totals, by = c("experiment", "region")) %>%
  arrange(experiment, region, context)


# Count DMRs by experiment, region, context, direction
summary <- DMR_table_data %>%
  group_by(experiment, region, context, direction) %>%
  summarise(count = n(), .groups = "drop")

# Pivot to wide format for context-direction pairs
summary_wide <- summary %>%
  unite("context_dir", context, direction) %>%
  pivot_wider(
    names_from = context_dir,
    values_from = count,
    values_fill = 0
  )

# Add total hyper, hypo, and overall total columns
summary_wide <- summary_wide %>%
  mutate(
    Total_Hyper = CHG_hyper + CHH_hyper + CpG_hyper,
    Total_Hypo  = CHG_hypo + CHH_hypo + CpG_hypo,
    Total       = Total_Hyper + Total_Hypo
  ) %>%
  relocate(starts_with("CHG"), starts_with("CHH"), starts_with("CpG"),
           Total_Hyper, Total_Hypo, Total)

write.csv(summary_wide, "DMR_summary_table_data.csv")

# install.packages("flextable")
library(flextable)
#install.packages("officer")
library(officer)
library(dplyr)

DMR_summary_table_data <- read.csv("DMR_summary_table_data.csv", sep = ',', header = T)

# Replace "-" with "–"
DMR_summary_table_data[DMR_summary_table_data == "0"] <- "–"

# Find where the experiment changes (e.g., row index before second experiment starts)
exp_change_row <- 8

# Set global defaults once (outside the ft pipeline)
set_flextable_defaults(background.color = "transparent")

# Format the table
ft <- flextable(DMR_summary_table_data) %>%
  # Relabel headers
  set_header_labels(
    Experiment = "Experiment", Region = "Region",
    CpG.hypo = "CpG Hypo", CpG.hyper = "CpG Hyper",
    CHG.hypo = "CHG Hypo", CHG.hyper = "CHG Hyper",
    CHH.hypo = "CHH Hypo", CHH.hyper = "CHH Hyper",
    Total.Hypo = "Total Hypo", Total.Hyper = "Total Hyper", Total = "Total"
  ) %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  theme_booktabs() %>%
  
  # Bold the header row
  bold(part = "header") %>%
  
  # Bold the 'Total Hypo', 'Total Hyper', and 'Total' columns
  bold(j = c("Total.Hypo", "Total.Hyper", "Total"), part = "body") %>%
  
  # Merge identical values in the 'experiment' column
  merge_v(j = "Experiment") %>%
  
  # Center vertically merged cells for cleaner look
  valign(j = "Experiment", valign = "center") %>%
  
  # Add horizontal line below last row of first experiment
  hline(i = exp_change_row, border = fp_border(color = "black", width = 1)) %>%
  border(i = 1, border.top = fp_border(color = "black", width = 1.5)) %>%
  border(i = nrow(DMR_summary_table_data), border.bottom = fp_border(color = "black", width = 1.5)) %>%
  bg(i = seq(2, nrow(DMR_summary_table_data), by = 2), bg = "#F9F9F9", part = "body") %>%
  fontsize(j = NULL, size = 11, part = "body") %>%
  fontsize(part = "header", size = 12) %>%
  bg(j = c("Total.Hypo", "Total.Hyper", "Total"), bg = "#EFEFEF", part = "body")

ft <- set_caption(ft, caption = "Table 1. Counts of differentially methylated regions (DMRs) by experiment, genomic feature, and cytosine context.")

ft <- ft %>%
  # Shrink column widths to fit ~6.3 inches total
  width(j = "Experiment", width = 0.6) %>%
  width(j = "Region", width = 1.0) %>%
  width(j = c("CpG.hypo", "CpG.hyper",
              "CHG.hypo", "CHG.hyper",
              "CHH.hypo", "CHH.hyper"), width = 0.5) %>%
  width(j = c("Total.Hypo", "Total.Hyper", "Total"), width = 0.55) %>%
  
  # Font size and padding for compactness
  fontsize(part = "header", size = 9) %>%
  fontsize(part = "body", size = 8) %>%
  padding(part = "all", padding.top = 1, padding.bottom = 1, padding.left = 1, padding.right = 1)

ft

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft, align = "left")
print(doc, target = "DMR_summary_table.docx")


