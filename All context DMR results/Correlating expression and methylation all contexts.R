# load packages:
{
  library(dplyr) # for data manipulation  
  library(ggplot2) # for plotting  
  library(data.table) # for data table manipulation
  library(plyranges)     # For genomic ranges manipulation
  library(GenomicRanges) # For handling GRanges objects
  library(genomation)    # For annotateWithFeatures and related tools
  library(ggpubr)
  library(tidyr)
  library(flextable)
  library(officer)
  library(ggh4x)
}


# load data
# set worrking wd
getwd()
setwd("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/Methylation_Analyses/All context DMR results")
getwd()

# Data ----
# Load data (created using vlookup on methyulation and expression data, capturing all relevant data)
data <- read.table("Corelating Methylation and Expression All context data.txt", header = T, sep = '\t')
data$region <- factor(data$region)
data$context <- factor(data$context)
data$experiment <- factor(data$experiment)

levels(data$region)[levels(data$region) == "dispersed_repeat"] <- "Repeat"
levels(data$region)[levels(data$region) == "downstream_region"] <- "Downstream"
levels(data$region)[levels(data$region) == "exon"] <- "Exon"
levels(data$region)[levels(data$region) == "five_prime_UTR"] <- "5` UTR"
levels(data$region)[levels(data$region) == "gene"] <- "Gene"
levels(data$region)[levels(data$region) == "ncRNA_gene"] <- "ncRNA"
levels(data$region)[levels(data$region) == "promoter"] <- "Promoter"
levels(data$region)[levels(data$region) == "three_prime_UTR"] <- "3` UTR"


## exp1 correlation plot ----
exp1_corr_plot <- ggplot(data = data[data$experiment == "exp1",],
                         aes(y = exp_l2fc, 
                             x = meth.diff, 
                             #shape = context, 
                             col = "black",
                             fill = context)) + 
  
  # Points (no colour)
  geom_point(#color = "#585858", 
    alpha = 0.7, 
    size = 1.5,
    shape = 21) + 
  
  # Correlation lines (coloured by context)
  geom_smooth(aes(color = context),
              method = "lm", 
              se = FALSE, 
              linetype = "solid") +
  
  # Custom line colours for context
  scale_color_manual(
    name = "Context",
    values = c(
      "CpG" = "#855c01",  # orange
      "CHG" = "blue",  # blue
      "CHH" = "#d61c2e"   # green
    )
  ) +
  
  scale_fill_manual(
    name = "Context",
    values = c(
      "CpG" = "#f0c565",
      "CHG" = "#83d0fc",
      "CHH" = "#eb6e7b"
    )
  ) +
  
  # Theme and reference lines
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Facet by region
  facet_wrap(~context~region, 
             scales = "free_x", 
             nrow = 3
  ) +
  
  # Label Axes better
  labs(
    x = "Differential methylation %",
    y = expression("Log"[2] * " Fold Change")
  ) + 
  
  # panel spacing
  theme(
    panel.spacing = unit(1, "lines"), # increase to add more vertical and horizontal space
    strip.text = element_text(size = 7, face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_text(size = 7)  # or smaller if needed
  )



exp1_corr_plot

ggsave("exp1_Correlation_Plot.png", exp1_corr_plot, width = 24, height = 14, units = "cm")
 
# Get correlation results ----
cor_results_pearson <- data %>%
  filter(experiment == "exp1") %>%
  group_by(region, context) %>%
  summarise(
    r = cor(meth.diff, exp_l2fc, 
            method = "pearson", 
            use = "complete.obs"),
    p = cor.test(meth.diff, exp_l2fc, 
                 method = "pearson")$p.value,
    .groups = "drop"
  )

cor_results <- cor_results_pearson %>%
  mutate(
    r = round(r, 3),
    p = signif(p, 3),
    sig = case_when(
      p < 0.001 ~ "****",
      p < 0.01  ~ " *** ",
      p < 0.05  ~ "  **",
      p < 0.1 ~ "  *",
      TRUE      ~ ""
    )
  )

# cor_results_exp2 <- cor_results %>%
#   mutate(
#     label = sig,
#     x = 0,         # center (you can adjust based on your data range)
#     y = Inf        # place at top of each panel
#   )

cor_results_exp1 <- cor_results %>%
  mutate(
    label_x = 0,
    label_y = 2.8
  )


print(cor_results_exp1, n = 24)

# Optional: write to CSV
write.csv(cor_results, "exp1_region_correlations.csv", row.names = FALSE)

# Add sigs to plot
# Adding astericks of sig to plot
exp1_corr_plot_sig <- exp1_corr_plot + 
  
  geom_text(data = cor_results_exp1,
            aes(x = label_x,
                y = label_y,
                label = sig),
            inherit.aes = F,
            hjust = 1, 
            vjust = 1,
            size = 7)

exp1_corr_plot_sig

ggsave("exp1_Correlation_Plot_sig.png", exp1_corr_plot_sig, width = 24, height = 14, units = "cm")


## exp2 correlation plot ----
exp2_corr_plot <- ggplot(data = data[data$experiment == "exp2",],
                         aes(y = exp_l2fc, 
                             x = meth.diff, 
                             #shape = context, 
                             col = "black",
                             fill = context)) + 
  
  # Points (no colour)
  geom_point(#color = "#585858", 
    alpha = 0.7, 
    size = 1.5,
    shape = 21) + 
  
  # Correlation lines (coloured by context)
  geom_smooth(aes(color = context),
              method = "lm", 
              se = FALSE, 
              linetype = "solid") +
  
  # Custom line colours for context
  scale_color_manual(
    name = "Context",
    values = c(
      "CpG" = "#855c01",  # orange
      "CHG" = "blue",  # blue
      "CHH" = "#d61c2e"   # green
    )
  ) +
  
  scale_fill_manual(
    name = "Context",
    values = c(
      "CpG" = "#f0c565",
      "CHG" = "#83d0fc",
      "CHH" = "#eb6e7b"
    )
  ) +
  
  # Theme and reference lines
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Facet by region
  facet_wrap(~context~region, 
             scales = "free_x", 
             nrow = 3
  ) +
  
  # Label Axes better
  labs(
    x = "Differential methylation %",
    y = expression("Log"[2] * " Fold Change")
  ) + 
  
  
  # panel spacing
  theme(
    panel.spacing = unit(1, "lines"), # increase to add more vertical and horizontal space
    strip.text = element_text(size = 7, face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_text(size = 7)  # or smaller if needed
  )

exp2_corr_plot

ggsave("exp2_Correlation_Plot.png", exp2_corr_plot, width = 24, height = 14, units = "cm")

# Get correlation results ----
cor_results_pearson <- data %>%
  filter(experiment == "exp2") %>%
  group_by(region, context) %>%
  summarise(
    r = cor(meth.diff, exp_l2fc, 
            method = "pearson", 
            use = "complete.obs"),
    p = cor.test(meth.diff, exp_l2fc, 
                 method = "pearson")$p.value,
    .groups = "drop"
  )

cor_results <- cor_results_pearson %>%
  mutate(
    r = round(r, 3),
    p = signif(p, 3),
    sig = case_when(
      p < 0.001 ~ "****",
      p < 0.01  ~ " *** ",
      p < 0.05  ~ " **",
      p < 0.1 ~ "     *",
      TRUE      ~ ""
    )
  )

# cor_results_exp2 <- cor_results %>%
#   mutate(
#     label = sig,
#     x = 0,         # center (you can adjust based on your data range)
#     y = Inf        # place at top of each panel
#   )

cor_results_exp2 <- cor_results %>%
  mutate(
    label_x = 0,
    label_y = 2.8
  )


print(cor_results_exp2, n = 24)

# Optional: write to CSV
write.csv(cor_results, "exp2_region_correlations.csv", row.names = FALSE)

# Adding astericks of sig to plot
exp2_corr_plot_sig <- exp2_corr_plot + 
  
  geom_text(data = cor_results_exp2,
            aes(x = label_x,
                y = label_y,
                label = sig),
            inherit.aes = F,
            hjust = 1, 
            vjust = 1,
            size = 7)

exp2_corr_plot_sig

ggsave("exp2_Correlation_Plot_sig.png", exp2_corr_plot_sig, width = 24, height = 14, units = "cm")

# Making a Table for the supplementary materials: ----

# Add experiment data to each cor_result df
cor_results_exp1$experiment <- "Acute"
cor_results_exp2$experiment <- "Primed"
  # join
  cor_results_combined <- rbind(cor_results_exp1, cor_results_exp2)
  


# remove unnecessary variables & rearrange
print(cor_results_combined, n = 48)
cor_results_combined_refined <- cor_results_combined[-c(6:8)]
cor_results_combined_refined <- cor_results_combined_refined[c(6,1,2,3,4,5)]
print(cor_results_combined_refined)



# Find where the experiment changes (e.g., row index before second experiment starts)
exp_change_row <- 24

# Set global defaults once (outside the ft pipeline)
set_flextable_defaults(background.color = "transparent")

# Format the table
ft <- flextable(cor_results_combined_refined) %>%
  # Relabel headers
  set_header_labels(
    experiment = "Experiment", 
    region = "Region",
    context = "Context",
    r = "R",
    p = "P-value",
    sig = "Sig"
  ) %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  theme_booktabs()  %>%

  
  # Bold the header row
  bold(part = "header") %>%
  

  # Bold the 'Total Hypo', 'Total Hyper', and 'Total' columns
  #bold(j = c("Total.Hypo", "Total.Hyper", "Total"), part = "body") %>%
  
  # Merge identical values in the 'experiment' column
  merge_v(j = "experiment") %>%
  
  # Center vertically merged cells for cleaner look
  valign(j = "experiment", valign = "top") %>%
  
  # merge identical value in the 'region' column
  merge_v(j="region") %>%
  valign(j = "region", valign = "top") %>%

  # Add horizontal line below last row of first experiment
  hline(i = exp_change_row, border = fp_border(color = "black", width = 1)) %>%
  border(i = 1, border.top = fp_border(color = "black", width = 1.5)) %>%
  border(i = nrow(cor_results_combined_refined), border.bottom = fp_border(color = "black", width = 1.5)) %>%
  bg(i = seq(2, nrow(cor_results_combined_refined), by = 2), bg = "#F9F9F9", part = "body") %>%
  fontsize(j = NULL, size = 11, part = "body") %>%
  fontsize(part = "header", size = 12)
  #bg(j = c("Total.Hypo", "Total.Hyper", "Total"), bg = "#EFEFEF", part = "body")

ft <- set_caption(ft, caption = "Table 2. Correlation of regionalized methylation and associated gene expression.")

ft <- ft %>%
  # Shrink column widths to fit ~6.3 inches total
# width(j = "experiment", width = 0.6) %>%
# width(j = "region", width = 1.0) %>%
# width(j = c("CpG.hypo", "CpG.hyper",
#             "CHG.hypo", "CHG.hyper",
#             "CHH.hypo", "CHH.hyper"), width = 0.5) %>%
# width(j = c("Total.Hypo", "Total.Hyper", "Total"), width = 0.55) %>%
  
  # Font size and padding for compactness
  fontsize(part = "header", size = 9) %>%
  fontsize(part = "body", size = 8) %>%
  padding(part = "all", padding.top = 1, padding.bottom = 1, padding.left = 1, padding.right = 1)

ft

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft, align = "left")
print(doc, target = "Meth-Expr_table.docx")
