# Load required packages
library(plyranges)     # For genomic ranges manipulation
library(GenomicRanges) # For handling GRanges objects
library(qqman)         # For Q-Q and Manhattan plots (not directly used, but loaded in original code)
library(dplyr)         # For data manipulation
library(ggplot2)       # For plotting
library(data.table)    # For fast I/O on data tables
library(genomation)    # For annotateWithFeatures and related tools

#--------------------------
# Data Input
#--------------------------
# Read in raw coverage data; these files contain CpG methylation info
raw_cov <- fread("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/methylation_summarisation/meth_CpG_cov_reads")
colnames(raw_cov) <- c("chrom", "start", "end", "pc_meth", "count_meth", "count_un")

# Read genome annotation (GFF3) which contains various feature types (genes, exons, promoters, etc.)
genome <- read_gff3("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/methylation_summarisation/combined_annotations.gff3")

#--------------------------
# Methylation Analysis
#--------------------------
# Perform binomial tests to determine which CpGs are significantly methylated
cov_data <- raw_cov %>%
  mutate(coverage = count_meth + count_un) %>%
  rowwise() %>%
  mutate(pv_meth = binom.test(count_meth, coverage, 0.01, alternative = "greater")$p.value) %>%
  ungroup() %>%
  mutate(
    pv_meth_c = p.adjust(pv_meth, method = "fdr"),
    is_meth   = ifelse(pv_meth_c < 0.05, 1, 0)
  )

# Calculate overall mean methylation rate for the sample
mean_meth_rate <- sum(cov_data$is_meth) / nrow(cov_data)
write.table(mean_meth_rate, "mean_cpg_meth.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Convert cov_data to a GRanges object for annotation
cov_gr <- GRanges(
  seqnames = cov_data$chrom,
  ranges   = IRanges(start = cov_data$start, end = cov_data$end),
  strand   = "*",
  is_meth  = cov_data$is_meth
)

#--------------------------
# Feature Annotation Setup
#--------------------------
# Extract features of interest from the genome. 
# Here, we pick out specific feature types.
genes             <- genome[genome$type == "gene"]
exons             <- genome[genome$type == "exon"]
promoters         <- genome[genome$type == "promoter"]
dispersed_repeats <- genome[genome$type == "dispersed_repeat"]
downstream_regions<- genome[genome$type == "downstream_region"]

# Combine selected features into a GRangesList
featuresList <- GRangesList(
  gene              = genes,
  exon              = exons,
  promoter          = promoters,
  dispersed_repeat  = dispersed_repeats,
  downstream_region = downstream_regions
)

# Annotate the target CpGs (cov_gr) with the selected features
annotation_result <- annotateWithFeatures(target = cov_gr, features = featuresList, intersect.chr = TRUE)

#--------------------------
# Summarizing Annotation Results
#--------------------------
# Extract annotation stats (percentage of CpGs in each feature)
annotation_stats <- getTargetAnnotationStats(annotation_result, percentage = TRUE)
annotation_df    <- data.frame(feature_type = names(annotation_stats),
                               percentage   = as.numeric(annotation_stats))

# Plot the annotation as a pie chart with ggplot2 for flexibility
summarising_piechart <- ggplot(annotation_df, aes(x = "", y = percentage, fill = feature_type)) +
  geom_col(width = 1, color = "black") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(percentage, "%")), position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "Distribution of CpGs - Sample 3-3", fill = "Feature Type") +
  theme_minimal() +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

summarising_piechart

ggsave("Sample_3-3_D_Distribution_of_CpGs.png", summarising_piechart)

# Optionally, you can still save a base R plot if needed:
plotTargetAnnotation(annotation_result, precedence = TRUE, main = "Distribution of CpGs - Sample 3-3")
png("TargetAnnotationPlot.png", width = 1200, height = 800, res = 150)
plotTargetAnnotation(annotation_result, precedence = TRUE, main = "Distribution of CpGs - Sample 3-3")
dev.off()

#--------------------------
# Mean Methylation by Feature Type
#--------------------------
# Assign each CpG a single feature type (first match found)
feature_assignment <- rep(NA_character_, length(cov_gr))
for (f in names(featuresList)) {
  ov <- findOverlaps(cov_gr, featuresList[[f]])
  q_hits <- queryHits(ov)
  unassigned <- is.na(feature_assignment[q_hits])
  feature_assignment[q_hits[unassigned]] <- f
}

# Add feature_type to raw_cov
raw_cov$feature_type <- feature_assignment

# Calculate mean methylation per feature type, with coverage filter
mean_meth_rate_feature <- raw_cov %>%
  mutate(coverage = count_meth + count_un) %>%
  filter(coverage %in% 5:29) %>%
  group_by(feature_type) %>%
  summarise(
    n          = n(),
    total_meth = sum(count_meth),
    mean       = total_meth / n
  )

write.table(mean_meth_rate_feature, "mean_cpg_meth_by_feature.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Plot the mean methylation by feature type as a bar chart
mean_methylation_by_feature <- ggplot(mean_meth_rate_feature, aes(x = feature_type, y = mean, fill = feature_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Mean Methylation by Feature Type",
    x = "Feature Type",
    y = "Mean Methylation"
  ) +
  theme(
    legend.position = "none", # Hide legend if feature_type names are clear enough as x-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if many feature types
  )

mean_methylation_by_feature

#--------------------------
# Proportion of Methylated Cytosines by Feature Type
#--------------------------
# Convert cov_gr to a data frame with assigned features
mcols(cov_gr)$feature_type <- feature_assignment
cov_df <- as.data.frame(cov_gr)

# Calculate the proportion of methylated CpGs per feature
prop_meth_by_feature <- cov_df %>%
  filter(!is.na(feature_type)) %>%
  group_by(feature_type) %>%
  summarise(
    total_cpgs   = n(),
    meth_cpgs    = sum(is_meth, na.rm = TRUE),
    proportion_meth = meth_cpgs / total_cpgs
  )

# Plot the proportion of methylated CpGs as a bar chart
proportion_methylated <- ggplot(prop_meth_by_feature, aes(x = feature_type, y = proportion_meth, fill = feature_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of Methylated Cytosines by Feature Type",
       x = "Feature Type",
       y = "Proportion of Methylated CpGs") +
  coord_flip()

ggsave("proportion_methylated_Sample_3-3_D.png", proportion_methylated)
