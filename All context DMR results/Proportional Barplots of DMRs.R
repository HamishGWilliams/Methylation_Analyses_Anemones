# load packages:
{
  library(dplyr) # for data manipulation  
  library(ggplot2) # for plotting  
  library(data.table) # for data table manipulation
  library(plyranges)     # For genomic ranges manipulation
  library(GenomicRanges) # For handling GRanges objects
  library(genomation)    # For annotateWithFeatures and related tools
}

# set worrking wd
getwd()
setwd("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/Methylation_Analyses/All context DMR results")
getwd()

# Load in Data
exp1_DMRs <- read.table("exp1_DMRs.txt", header = T, sep = "\t")
exp2_DMRs <- read.table("exp2_DMRs.txt", header = T, sep = "\t")

data <- read_gff3("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/methylation_summarisation/combined_annotations.gff3")
data <- as(data, 'data.frame')
data$region <- data$type
summary(data$region)
levels(data$region)[levels(data$region) == "dispersed_repeat"] <- "Repeat"
levels(data$region)[levels(data$region) == "downstream_region"] <- "Downstream"
levels(data$region)[levels(data$region) == "exon"] <- "Exon"
levels(data$region)[levels(data$region) == "five_prime_UTR"] <- "5` UTR"
levels(data$region)[levels(data$region) == "gene"] <- "Gene"
levels(data$region)[levels(data$region) == "ncRNA_gene"] <- "ncRNA"
levels(data$region)[levels(data$region) == "promoter"] <- "Promoter"
levels(data$region)[levels(data$region) == "three_prime_UTR"] <- "3` UTR"

genome <- subset(data[data$region == c("Repeat",
                                    "Downstream",
                                    "Exon",
                                    "5` UTR",
                                    "Gene",
                                    "ncRNA",
                                    "Promoter",
                                    "3` UTR"),])

summary(genome$region)

# change order of all dataframes regions to match:
desired_order <- c("Repeat","Downstream","Exon","5` UTR", "Gene", "ncRNA","Promoter","3` UTR")

exp1_DMRs$region <- factor(exp1_DMRs$region, levels = desired_order)
exp2_DMRs$region <- factor(exp2_DMRs$region, levels = desired_order)
genome$region <- factor(genome$region, levels = desired_order)

# Making barplots ----

# exp1 ----
# subset to only DMRs
exp1_barplot_data <- subset(exp1_DMRs[exp1_DMRs$p_fdr <= 0.1,]) 
exp1_barplot_data$experiment <- as.factor(exp1_barplot_data$experiment)
# Count occurrences by region and context
exp1_barplot_data <- exp1_barplot_data %>%
  count(region, context, experiment)
# Calculate proportions relative to full dataset
exp1_barplot_data <- exp1_barplot_data %>%
  mutate(proportion = (n / sum(n)))
# geom_text on top of bars
bar_tops_exp1 <- exp1_barplot_data %>%
  group_by(region) %>%
  summarise(total_prop = sum(proportion),
            total_n = sum(n))
# make plot
exp1_barplot <- ggplot(exp1_barplot_data, aes(x = region, y = proportion, fill = context)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(
    x = "Region",
    y = "Proportion"# , title = "Counts of Genomic Regions - Experiment 1"
  ) +
  # Add total count at the top of each bar
  geom_text(data = bar_tops_exp1, aes(x = region, y = total_prop + 0.04, label = total_n),
            vjust = -0.3, inherit.aes = FALSE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1)) 

# exp2 ----
# subset to only DMRs
exp2_barplot_data <- subset(exp2_DMRs[exp2_DMRs$p_fdr <= 0.1,]) 
exp2_barplot_data$experiment <- as.factor(exp2_barplot_data$experiment)
# Count occurrences by region and context
exp2_barplot_data <- exp2_barplot_data %>%
  count(region, context, experiment)
# Calculate proportions relative to full dataset
exp2_barplot_data <- exp2_barplot_data %>%
  mutate(proportion = (n / sum(n)))
# geom_text on top of bars
bar_tops_exp2 <- exp2_barplot_data %>%
  group_by(region) %>%
  summarise(total_prop = sum(proportion),
            total_n = sum(n))
# make plot
exp2_barplot <- ggplot(exp2_barplot_data, aes(x = region, y = proportion, fill = context)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(
    x = "Region",
    y = "Proportion"# , title = "Counts of Genomic Regions - Experiment 2"
  ) +
  # Add total count at the top of each bar
  geom_text(data = bar_tops_exp2, aes(x = region, y = total_prop + 0.04, label = total_n),
            vjust = -0.3, inherit.aes = FALSE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1)) 

# genome ----
# Count occurrences by region and context
genome_barplot_data <- genome %>%
  count(region)
# Calculate proportions relative to full dataset
genome_barplot_data <- genome_barplot_data %>%
  mutate(proportion = (n / sum(n)))
# geom_text on top of bars
bar_tops_genome <- genome_barplot_data %>%
  group_by(region) %>%
  summarise(total_prop = sum(proportion),
            total_n = sum(n))

# make plot
genome_barplot <- ggplot(genome_barplot_data, aes(x = region, y = proportion)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(
    x = "Region",
    y = "Proportion" #, title = "Counts of Genomic Regions - genome"
  ) +
  # Add total count at the top of each bar
  geom_text(data = bar_tops_genome, aes(x = region, y = total_prop + 0.04, label = total_n),
            vjust = -0.3, inherit.aes = FALSE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.2,0.4,0.6,0.8,1)) 

# all plots
exp1_bar_plot <- exp1_barplot + scale_fill_manual(values = c("CpG" = "#E69F00", "CHG" = "#56B4E9", "CHH" = "#009E73")) 
exp2_bar_plot <- exp2_barplot + scale_fill_manual(values = c("CpG" = "#E69F00", "CHG" = "#56B4E9", "CHH" = "#009E73")) 
genome_barplot

ggsave("Exp1_DMRs_Barplot.png", exp1_bar_plot, width = 10, height = 12, units = "cm")
ggsave("Exp2_DMRs_Barplot.png", exp2_bar_plot, width = 10, height = 12, units = "cm")
ggsave("Genome_Barplot.png", genome_barplot, width = 10.5, height = 12, units = "cm")



# Proportion tests ()chi-squared:----

## exp 1: ----

sum(genome_barplot_data$n)
sum(exp1_barplot_data$n)

  # Repeats
prop.test(x = c(122277, 67), n = c(195009,145))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(122277, 67) out of c(195009, 147)
  # X-squared = 17.693, df = 1, p-value = 5.821e-05
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #    0.08033586 0.24959141
  # sample estimates:
  #   prop 1    prop 2 
  # 0.6270326 0.4620690  

  # Downstream
prop.test(x = c(5959, 1), n = c(195009,145))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(5959, 1) out of c(195009, 147)
  # X-squared = 2.0548, df = 1, p-value = 0.1574
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   0.007041318 0.040468368
  # sample estimates:
  #   prop 1      prop 2 
  # 0.030557564 0.006896552  
  # 
  # Warning message:
  #   In prop.test(x = c(5959, 1), n = c(195009, 147)) :
  #   Chi-squared approximation may be incorrect

  # Exon
prop.test(x = c(47918, 10), n = c(195009,145))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(47918, 10) out of c(195009, 147)
  # X-squared = 24.084, df = 1, p-value = 1.257e-06
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   0.1320172 0.2214958
  # sample estimates:
  #   prop 1     prop 2 
  # 0.24572199 0.06896552  

  # 5` UTR`
    # NO 5` UTRS in exp1 DMRs

  # Gene
prop.test(x = c(6244, 60), n = c(195009,145))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(6244, 60) out of c(195009, 147)
  # X-squared = 652.87, df = 1, p-value < 2.2e-16
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   -0.4590045 -0.2932840
  # sample estimates:
  #   prop 1     prop 2 
  # 0.03201904 0.41379310 
  # 
  # Warning message:
  #   In prop.test(x = c(6244, 60), n = c(195009, 147)) :
  #   Chi-squared approximation may be incorrect

  # Promoter
prop.test(x = c(5855, 7), n = c(195009,145))
  # data:  c(5855, 7) out of c(195009, 147)
  # X-squared = 1.0153, df = 1, p-value = 0.2966
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #    -0.05659936  0.02009614
  # sample estimates:
  #   prop 1     prop 2 
  # 0.03002426 0.04827586  
  # 
  # Warning message:
  #   In prop.test(x = c(5855, 7), n = c(195009, 147)) :
  #   Chi-squared approximation may be incorrect

  # 3` UTR`
prop.test(x = c(2805, 2), n = c(195009,145))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(2805, 2) out of c(195009, 147)
  # X-squared = 8.7598e-29, df = 1, p-value = 1
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   -0.01873458  0.02029160
  # sample estimates:
  #   prop 1     prop 2 
  # 0.01438395 0.01360544 
  # 
  # Warning message:
  #   In prop.test(x = c(2805, 2), n = c(195009, 147)) :
  #   Chi-squared approximation may be incorrect


# exp2:
  # repeats
prop.test(x = c(122277, 105), n = c(195009,230))

  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(122277, 105) out of c(195009, 230)
  # X-squared = 27.83, df = 1, p-value = 1.325e-07
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   0.1039253 0.2370964
  # sample estimates:
  #   prop 1    prop 2 
  # 0.6270326 0.4565217 

  # downstream
prop.test(x = c(5959, 24), n = c(195009,230))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(5959, 24) out of c(195009, 230)
  # X-squared = 39.662, df = 1, p-value = 3.02e-10
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   -0.11548308 -0.03209744
  # sample estimates:
  #   prop 1     prop 2 
  # 0.03055756 0.10434783 

  # exons
prop.test(x = c(47918, 18), n = c(195009,230))

  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(47918, 18) out of c(195009, 230)
  # X-squared = 33.88, df = 1, p-value = 5.862e-09
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   0.1305216 0.2044006
  # sample estimates:
  #   prop 1     prop 2 
  # 0.24572199 0.07826087 

  # genes
prop.test(x = c(6244, 67), n = c(195009,230))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(6244, 67) out of c(195009, 230)
  # X-squared = 485.5, df = 1, p-value < 2.2e-16
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   -0.3201872 -0.1983834
  # sample estimates:
  #   prop 1     prop 2 
  # 0.03201904 0.29130435 

  # ncRNA
prop.test(x = c(1000, 3), n = c(195009,230))

  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(1000, 3) out of c(195009, 230)
  # X-squared = 1.4805, df = 1, p-value = 0.2237
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   -0.024758652  0.008927632
  # sample estimates:
  #   prop 1      prop 2 
  # 0.005127968 0.013043478 
  # 
  # Warning message:
  #   In prop.test(x = c(1000, 3), n = c(195009, 230)) :
  #   Chi-squared approximation may be incorrect

  # promoter
prop.test(x = c(5855, 12), n = c(195009,230))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(5855, 12) out of c(195009, 230)
  # X-squared = 3.1442, df = 1, p-value = 0.0762
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   -0.053075349  0.008776034
  # sample estimates:
  #   prop 1     prop 2 
  # 0.03002426 0.05217391 

  # 5` UTR`
prop.test(x = c(2951, 1), n = c(195009,230))
  # 2-sample test for equality of proportions with continuity correction
  # 
  # data:  c(2951, 1) out of c(195009, 230)
  # X-squared = 1.1432, df = 1, p-value = 0.285
  # alternative hypothesis: two.sided
  # 95 percent confidence interval:
  #   8.804856e-05 2.148157e-02
  # sample estimates:
  #   prop 1      prop 2 
  # 0.015132635 0.004347826 
  # 
  # Warning message:
  #   In prop.test(x = c(2951, 1), n = c(195009, 230)) :
  #   Chi-squared approximation may be incorrect


