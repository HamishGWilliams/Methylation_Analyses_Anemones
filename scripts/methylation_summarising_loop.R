# ------------------------------
# Load in all necessary packages
# ------------------------------
library(plyranges)     # For genomic ranges manipulation
library(GenomicRanges) # For handling GRanges objects
library(qqman)         # For Q-Q and Manhattan plots (not directly used, but loaded in original code)
library(dplyr)         # For data manipulation
library(ggplot2)       # For plotting
library(data.table)    # For fast I/O on data tables
library(genomation)    # For annotateWithFeatures and related tools

# ------------------------------
# Initialise the loop with parameters
# ------------------------------
base_path <- "/uoa/scratch/users/r02hw22/Methylation_Analyses/Data2/Trimmed" # path to data
setwd(base_path) # change path to this directory to access samples

# Inside this directory is the data which is at the sample level. 
# Therefore for each "experiment" we need to manually separate out the 
# different samples to that each analyses is completed seperately for 
# collating each experiment's data together.

# Experiment 1 samples
exp1_samples <- c("Sample_3-3_D",
                  "Sample_9-9_D",
                  "Sample_13-13_D",
                  "Sample_15-15_D",
                  "Sample_16-16_D",
                  "Sample_17-17_D",
                  "Sample_22-22_D",
                  "Sample_27-27_D")

# Experiment 2 samples
exp2_samples <- c("Sample_6-6_D",
                  "Sample_7-7_D",
                  "Sample_8-8_D",
                  "Sample_10-10_D",
                  "Sample_14-14_D",
                  "Sample_20-20_D",
                  "Sample_21-21_D",
                  "Sample_29-11_redo_D",
                  "Sample_32-18_redo_D",
                  "Sample_36-2_redo_D")

# Make an object to allow you to direct through each experiment separately
experiments = c("exp1","exp2")
# And another for each specific context
contexts = c("CpG","CHG","CHH")
# Load in genome
genome <- read_gff3("/uoa/scratch/users/r02hw22/Methylation_Analyses/Data2/combined_annotations.gff3")

# ----
# Start Loop

for (experiment in experiments) {
  if (experiment == "exp1") {
    file.list = exp1_samples
    sample.id = list("3","9","13","15","16","17","22","27")
    treatment = c(1, 0, 0, 1, 0, 0, 1, 1)
  } else if (experiment == "exp2") {
    file.list = exp2_samples
    sample.id = list("6", "7", "8", "10", "14", "20", "21", "29_11", "32_18", "36_2")
    treatment = c(0, 0, 1, 0, 1, 1, 1, 0, 0, 1)
  }
  
  for (context in contexts) {
    context = context                                                           # make an object to house the specific context for loop
      for (sample in file.list) {                                               # start a loop for each sample in experiment
        fpath <- paste0(base_path, "/", sample)                                 # make dynamic path to specific file
        setwd(fpath)                                                            # change to sample directory
        reads_file <- paste0("meth_", context, "_cov_reads")                    # dynamically named file to import
        sample_name <- gsub("meth_","_cov_reads", reads_file)                   # sample_name object to dynamically name items later
        raw_cov <- fread(reads_file)                                            # import cov_read file
        colnames(raw_cov) <- c("chrom", "start", "end", "pc_meth", "count_meth",# change column names
                               "count_un")
        
        # Perform methylation analysis (binomial test)
        cov_data <- raw_cov %>%
          mutate(coverage = count_meth + count_un) %>%
          rowwise() %>%
          mutate(pv_meth = binom.test(count_meth, coverage, 0.01, 
                                      alternative = "greater")$p.value) %>%
          ungroup() %>%
          mutate(
            pv_meth_c = p.adjust(pv_meth, method = "fdr"),
            is_meth   = ifelse(pv_meth_c < 0.05, 1, 0)
          )
        
        # Convert to GRanges
        cov_gr <- GRanges(
          seqnames = cov_data$chrom,
          ranges   = IRanges(start = cov_data$start, end = cov_data$end),
          strand   = "*",
          is_meth  = cov_data$is_meth
        )
        
        # Annotate with features
        annotation_result <- annotateWithFeatures(target = cov_gr, features = featuresList, intersect.chr = TRUE)
        
        # Assign feature type to each CpG
        feature_assignment <- rep(NA_character_, length(cov_gr))
        for (ft in names(featuresList)) {
          ov <- findOverlaps(cov_gr, featuresList[[ft]])
          q_hits <- queryHits(ov)
          unassigned <- is.na(feature_assignment[q_hits])
          feature_assignment[q_hits[unassigned]] <- ft
        }
        
        # Add feature_type to raw_cov
        raw_cov$feature_type <- feature_assignment
        
        # Calculate mean methylation per feature type (with coverage filter)
        mean_meth_rate_feature <- raw_cov %>%
          mutate(coverage = count_meth + count_un) %>%
          filter(coverage %in% 5:29) %>%
          group_by(feature_type) %>%
          summarise(
            n          = n(),
            total_meth = sum(count_meth),
            mean       = total_meth / n
          ) %>%
          # Add a column to indicate which sample these results come from
          mutate(sample = sample_name)
        
        # Append results to the combined list
        all_sample_results[[sample_name]] <- mean_meth_rate_feature
      } # now back at the context level
    
    # Combine all sample results into a single data frame
    combined_results <- bind_rows(all_sample_results)
    
    #--------------------------
    # Visualization: Box Plots
    #--------------------------
    # Now combined_results has columns: feature_type, n, total_meth, mean, and sample
    # We can create a box plot to show the distribution of average methylation per feature type across samples
    
    ggplot(combined_results, aes(x = feature_type, y = mean, fill = feature_type)) +
      geom_boxplot() +
      theme_minimal() +
      labs(
        title = "Distribution of Average Methylation per Feature Type Across Samples",
        x = "Feature Type",
        y = "Average Methylation"
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      coord_flip()
    
    # Save the box plot if you want
    boxplot_name <- paste0("Average_Methylation_per_feature_boxplot-", experiment, context)
    ggsave(boxplot_name, width = 8, height = 6, dpi = 300)
  }
}


