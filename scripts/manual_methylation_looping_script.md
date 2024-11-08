```
---- # Load Packages ----

library(methylKit) # for methylation analyses  
library(plyranges) # ...  
library(GenomicRanges) # for converting anf using files as GRanges  
library(qqman) # ...  
library(dplyr) # provides access to functions to manipulate data  
library(ggplot2) # for plotting purposes
library(data.table) # formatting data
library(Rsamtools)  # For bgzip compression
library(data.table) # manipulating dfs

# Set working directory
setwd("/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results")

# ---- Write a loop to generate file.list objects for each experiment and context ----

# Define the sample numbers for each experiment
exp1_samples <- c("3-3", "9-9", "13-13", "15-15", "16-16", "17-17", "22-22", "27-27")
exp2_samples <- c("6-6", "7-7", "8-8", "10-10", "14-14", "20-20", "21-21", "29-11_redo", "32-18_redo", "36-2_redo")

# Define the contexts
contexts <- c("CpG", "CHG", "CHH")

# Initialize empty lists to store the file paths
file.list.exp1 <- list(CpG = list(), CHG = list(), CHH = list())
file.list.exp2 <- list(CpG = list(), CHG = list(), CHH = list())

# Base path for the files
base_path <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Data2/Trimmed/"

# Function to generate file paths for a given sample list and experiment
generate_file_paths <- function(samples, contexts) {
  file_list <- list()
  
  # Loop over each context
  for (context in contexts) {
    # Initialize a list for the current context
    file_list[[context]] <- list()
    
    # Loop over each sample
    for (sample in samples) {
      # Construct the file path
      file_path <- paste0(base_path, "Sample_", sample, "_D/meth_", context, "_cov_reads")
      
      # Append the file path to the current context list
      file_list[[context]] <- append(file_list[[context]], file_path)
    }
  }
  
  return(file_list)
}

# Generate file paths for Experiment 1 and Experiment 2
file.list.exp1 <- generate_file_paths(exp1_samples, contexts)
file.list.exp2 <- generate_file_paths(exp2_samples, contexts)

# Now you have lists with file paths for Experiment 1 and 2
# Access them like file.list.exp1$CpG, file.list.exp1$CHG, file.list.exp1$CHH
# Same goes for file.list.exp2

# Example of how to print the generated lists
print(file.list.exp1$CpG)
print(file.list.exp2$CHG)

# ---- methRead all data into objects ----

# Sample IDs for each experiment
exp1_sample_ids <- list("3", "9", "13", "15", "16", "17", "22", "27")
exp2_sample_ids <- list("6", "7", "8", "10", "14", "20", "21", "29_11", "32_18", "36_2")

# Treatments for each experiment
exp1_treatment <- c(1, 0, 0, 1, 0, 0, 1, 1)
exp2_treatment <- c(0, 0, 1, 0, 1, 1, 1, 0, 0, 1)

# Define the base path for storing results
results_base_path <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results"

# List of contexts
contexts <- c("CpG", "CHG", "CHH")

# Create a list of sample details for each experiment
experiments <- list(
  exp1 = list(
    file.list = file.list.exp1,
    sample.id = exp1_sample_ids,
    treatment = exp1_treatment
  ),
  exp2 = list(
    file.list = file.list.exp2,
    sample.id = exp2_sample_ids,
    treatment = exp2_treatment
  )
)

# Loop through each experiment and context
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Create a directory name based on the experiment and context
    result_dir <- file.path(results_base_path, paste0(exp_name, "_", context))
    
    # Create the directory if it does not exist
    if (!dir.exists(result_dir)) {
      dir.create(result_dir, recursive = TRUE)
      cat("Created directory:", result_dir, "\n")
    }
    
    # Run the methRead function with the correct arguments
    meth_obj <- methRead(
      location = exp_data$file.list[[context]],  # Use 'location' instead of 'file.list'
      sample.id = exp_data$sample.id,
      assembly = "actinia",
      pipeline = "bismarkCoverage",
      treatment = exp_data$treatment,
      context = context,
      dbtype = "tabix",
      dbdir = result_dir
    )
    
    # Store the resulting methRead object in a dynamically named variable
    assign(paste0("methRead_", exp_name, "_", context), meth_obj)
    
    # Output the status
    cat("Processed:", exp_name, "with context:", context, "\n")
    cat("Results stored in:", result_dir, "\n")
  }
}
                   
# Loop through each experiment and context to generate raw methylation and coverage plots ----
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Get the directory name based on the experiment and context
    result_dir <- file.path(results_base_path, paste0(exp_name, "_", context))
    
    # Construct the dynamically named methRead object
    meth_obj_name <- paste0("methRead_", exp_name, "_", context)
    meth_obj <- get(meth_obj_name)
    
    # Generate the Raw Methylation plot and save to the appropriate directory
    png(file.path(result_dir, paste0("raw_meth_", exp_name, "_", context, ".png")), 
        width = 12, height = 8, units = "in", res = 300)
    lapply(meth_obj, getMethylationStats, plot = TRUE)
    dev.off()
    
    # Generate the Raw Coverage plot and save to the appropriate directory
    png(file.path(result_dir, paste0("raw_coverage_", exp_name, "_", context, ".png")), 
        width = 12, height = 8, units = "in", res = 300)
    lapply(meth_obj, getCoverageStats, plot = TRUE)
    dev.off()
    
    # Output the status
    cat("Generated plots for:", exp_name, "with context:", context, "\n")
  }
}

# Filter the data, standardize by median coverage, and then unite the results ----
# Loop through each experiment and context to filter, normalize, and unite methylation data
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Get the directory name based on the experiment and context
    result_dir <- file.path(results_base_path, paste0(exp_name, "_", context))
    
    # Construct the dynamically named methRead object
    meth_obj_name <- paste0("methRead_", exp_name, "_", context)
    meth_obj <- get(meth_obj_name)
    
    # Step 1: Filter by coverage and save results to the appropriate directory
    filtered_obj <- filterByCoverage(meth_obj, 
                                     lo.count = 5, 
                                     hi.count = 32, 
                                     suffix = paste0(context, "_f"), 
                                     dbdir = result_dir)
    
    # Step 2: Normalize coverage and save results to the appropriate directory
    normalized_obj <- normalizeCoverage(filtered_obj, 
                                        method = "median", 
                                        dbdir = result_dir)
    
    # Step 3: Unite the methylation results and save to the appropriate directory
    united_obj <- methylKit::unite(normalized_obj, 
                                   destrand = FALSE, 
                                   min.per.group = 3L, 
                                   suffix = paste0(context, "_fu3"), 
                                   dbdir = result_dir)
    
    # Save the resulting united object dynamically
    assign(paste0("united_", exp_name, "_", context), united_obj)
    
    # Output the status
    cat("Processed filtering, normalization, and uniting for:", exp_name, "with context:", context, "\n")
    cat("Results stored in:", result_dir, "\n")
  }
}

## plotting clustering and PCAs ----
# Loop through each experiment and context to generate clustering and PCA plots
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Get the directory name based on the experiment and context
    result_dir <- file.path(results_base_path, paste0(exp_name, "_", context))
    
    # Construct the dynamically named united object
    united_obj_name <- paste0("united_", exp_name, "_", context)
    united_obj <- get(united_obj_name)
    
    # Generate the Clustering Plot and save to the appropriate directory
    png(file.path(result_dir, paste0("cluster_plot_", exp_name, "_", context, ".png")),
        width = 12, height = 8, units = "in", res = 300)
    clusterSamples(united_obj,
                   filterByQuantile = FALSE,
                   sd.threshold = 0.1,
                   dist = "correlation",
                   method = "ward.D",
                   plot = TRUE)
    dev.off()
    
    # Generate PCA Plots for different components and save to the appropriate directory
    for (comp_pair in list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))) {
      comp_str <- paste0(comp_pair[1], "-", comp_pair[2])
      png(file.path(result_dir, paste0("PCA_", comp_str, "_plot_", exp_name, "_", context, ".png")),
          width = 12, height = 8, units = "in", res = 300)
      PCASamples(united_obj, comp = comp_pair)
      dev.off()
    }
    
    # Output the status
    cat("Generated clustering and PCA plots for:", exp_name, "with context:", context, "\n")
  }
}

# Calculate Diff in Methylation ----
# Loop through each experiment and context to calculate differential methylation
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Get the directory name based on the experiment and context
    result_dir <- file.path(results_base_path, paste0(exp_name, "_", context))
    
    # Construct the dynamically named united object
    united_obj_name <- paste0("united_", exp_name, "_", context)
    united_obj <- get(united_obj_name)
    
    # Define a context-specific suffix for the differential methylation analysis
    diff_meth_suffix <- paste0(context, "_", exp_name, "_DiffMeth")
    
    # Calculate Differential Methylation and save to the appropriate directory
    diff_meth_obj <- calculateDiffMeth(
      united_obj,
      overdispersion = "MN",
      mc.cores = 32,
      suffix = diff_meth_suffix,
      dbdir = result_dir
    )
    
    # Save the resulting differential methylation object dynamically
    assign(paste0("DiffMeth_", exp_name, "_", context), diff_meth_obj)
    
    # Output the status
    cat("Calculated differential methylation for:", exp_name, "with context:", context, "\n")
    cat("Results stored in:", result_dir, "\n")
  }
}

# Filter results for diffMeth = 5 ----
# Loop through each experiment and context to retrieve differential methylation data
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Construct the dynamically named differential methylation object
    diff_meth_obj_name <- paste0("DiffMeth_", exp_name, "_", context)
    diff_meth_obj <- get(diff_meth_obj_name)
    
    # Retrieve the differential methylation data
    diff_meth_data <- getMethylDiff(
      diff_meth_obj,
      qvalue = 0.99,
      diff = 5,
      type = "all"
    )
    
    # Save the differential methylation data to a dynamically named object
    assign(paste0("diffMethData_", exp_name, "_", context), diff_meth_data)
    
    # Output the status
    cat("Retrieved differential methylation data for:", exp_name, "with context:", context, "\n")
  }
}

# Save various formats of the data, and join annotation data to methylation data ----

# Import genome file, and remove NAs
genome <- read_gff3("/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Data2/genome_plus_upstream_and_downstream.gff3")
    # CHANGE THIS TO FINAL ANNOTATED GENOME FILE ONCE IT HAS FINISHED RUNNING!
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
genome <- as(genome, "GRanges")

# Loop through each experiment and context to save data and perform annotation
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Construct the dynamically named differential methylation object
    diff_meth_obj_name <- paste0("diffMethData_", exp_name, "_", context)
    diff_meth_obj <- get(diff_meth_obj_name)
    
    # Convert to data frame and save in the environment
    diff_meth_df <- getData(diff_meth_obj) # changed to df
    diff_meth_df_name <- paste0("DiffMeth_", exp_name, "_", context, "_df") # write a name for the file
    assign(diff_meth_df_name, diff_meth_df)
    save(diff_meth_df, file = file.path(results_base_path, paste0("Actinia_DMBs_d5_", exp_name, "_", context, "_df.RData")))
    
    # Convert to GRanges object and save in the environment
    diff_meth_gr <- as(diff_meth_df, "GRanges")
    diff_meth_gr_name <- paste0("DiffMeth_", exp_name, "_", context, "_gr")
    assign(diff_meth_gr_name, diff_meth_gr)
    save(diff_meth_gr, file = file.path(results_base_path, paste0("Actinia_DMBs_d5_", exp_name, "_", context, "_gr.RData")))
    
    # Join with the genome annotation
    diff_meth_gr_ann <- join_overlap_left_directed(diff_meth_gr, genome)
    
    # Convert the joined annotated object to a data frame and save in the environment
    diff_meth_df_ann <- as(diff_meth_gr_ann, "data.frame") # change to df
    diff_meth_df_ann$p_fdr <- p.adjust(diff_meth_df_ann$pvalue, method = "fdr") # p-adjust
    diff_meth_df_ann$Parent <- as.character(diff_meth_df_ann$Parent) # fix class of variable
    diff_meth_df_ann_name <- paste0("DiffMeth_", exp_name, "_", context, "_df_ann")
    assign(diff_meth_df_ann_name, diff_meth_df_ann)
    
    # Save the final annotated data frame to a file
    save(diff_meth_df_ann, file = file.path(results_base_path, paste0("Actinia_DMBs_d5_", exp_name, "_", context, "_df_ann.RData")))
    
    # Output the status
    cat("Processed and saved data for:", exp_name, "with context:", context, "\n")
  }
}

# Make Volcano plots for each experiment and context ----
# Loop through each experiment and context to create and save volcano plots
for (exp_name in names(experiments)) {
  exp_data <- experiments[[exp_name]]
  
  for (context in contexts) {
    
    # Get the directory name based on the experiment and context
    result_dir <- file.path(results_base_path, paste0(exp_name, "_", context))
    
    # Construct the dynamically named differential methylation data frame object
    diff_meth_data_name <- paste0("DiffMeth_", exp_name, "_", context, "_df_ann")
    data <- get(diff_meth_data_name)
    
    # Prepare data for plotting
    data$diffmeth <- ''
    data$diffmeth[data$meth.diff >= 0] <- 'UP'
    data$diffmeth[data$meth.diff <= 0] <- 'DOWN'
    data$diffmeth <- as.factor(data$diffmeth)
    
    # Create the volcano plot using ggplot
    volcano_plot <- ggplot(data = data,
                           aes(x = meth.diff, y = -log10(p_fdr), col = diffmeth)) +
      geom_point() +
      geom_vline(xintercept = c(-5), col = "blue", linetype = "dashed") +
      geom_vline(xintercept = c(5), col = "red", linetype = "dashed") +
      geom_hline(yintercept = c(1), col = "black", linetype = "dashed") +
      theme_classic() +
      theme(
        axis.title.y = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1.1), color = "black"),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1.1), color = "black")
      ) +
      scale_color_manual(values = c("deepskyblue", "brown1"),
                         labels = c("Down Methylated", "Up Methylated")) +
      labs(
        x = "Differential methylation %", y = expression("-log"[10]*"p-adj")
      )
    
    # Save the volcano plot to the appropriate directory
    plot_filename <- paste0("volcano_plot_", exp_name, "_", context, ".png")
    ggsave(file.path(result_dir, plot_filename), volcano_plot, width = 18, height = 12, units = "cm")
    
    # Output the status
    cat("Generated and saved volcano plot for:", exp_name, "with context:", context, "\n")
  }
}


## Regionlization code ----

# Set export directory
export_dir <- "/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/results"  # Replace with your desired path

# Get the unique types in the 'genome' data
types <- unique(genome$type)

# Loop over each unique 'type'
for (current_type in types) {

  # reset starting position to results directory
  setwd(export_dir)
  
  # Create directory only if it doesn't already exist
  if (!dir.exists(export_dir, current_type)) {
    dir.create(export_dir, current_type)
  
  # Change to this new directory
  new_dir <- file.path(export_dir, current_type)
  setwd(new_dir)
 
}
  # Subset 'genome' to regions of the current type
  genome_type <- genome[genome$type == current_type]
  
  # Check if there are any regions of this type
  if (length(genome_type) == 0) {
    cat("No regions found for type:", current_type, "\n", file = log_file, append = TRUE)
    next  # Skip to the next type
  }
  
  cat("Processing regions of type:", current_type, "\n", file = log_file, append = TRUE)
  
  # Loop over each sample
  for (i in seq_along(act_CpG_raw_DB)) {
    # Get the methylRawDB object for the current sample
    methyl_raw_db <- act_CpG_raw_DB[[i]]
    sample_id <- methyl_raw_db@sample.id
    
    # Log sample processing start
    cat("Processing sample:", sample_id, "for type:", current_type, "\n", file = log_file, append = TRUE)
    
    # Retrieve data using getData()
    sample_data <- getData(methyl_raw_db)
    
    # Convert methylation data to GRanges
    methyl_gr <- GRanges(
      seqnames = sample_data$chr,
      ranges = IRanges(start = sample_data$start, end = sample_data$end),
      strand = sample_data$strand,
      coverage = sample_data$coverage,
      numCs = sample_data$numCs,
      numTs = sample_data$numTs
    )
    
    # Ensure chromosome names match, pruning sequences not in genome_type
    seqlevels(methyl_gr, pruning.mode = "coarse") <- seqlevels(genome_type)
    
    # Find overlaps between methylation sites and regions of current_type
    overlaps <- findOverlaps(
      query = genome_type,
      subject = methyl_gr,
      ignore.strand = TRUE  # Set to FALSE if strand-specific
    )
    
    # Log number of overlaps found
    num_overlaps <- length(overlaps)
    cat("Number of overlaps found for sample", sample_id, "and type", current_type, ":", num_overlaps, "\n", file = log_file, append = TRUE)
    
    # Save overlaps information to a file
    overlaps_info_file <- file.path(export_dir, paste0(sample_id, "_", current_type, "_overlaps_info.tsv"))
    
    overlaps_info <- data.frame(
      Region = as.character(seqnames(genome_type[queryHits(overlaps)])),
      Region_Start = start(genome_type[queryHits(overlaps)]),
      Region_End = end(genome_type[queryHits(overlaps)]),
      Methylation_Site = as.character(seqnames(methyl_gr[subjectHits(overlaps)])),
      Site_Position = start(methyl_gr[subjectHits(overlaps)]),
      Coverage = mcols(methyl_gr)$coverage[subjectHits(overlaps)]
    )
    
    write.table(
      overlaps_info,
      file = overlaps_info_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    cat("Overlaps information for sample", sample_id, "and type", current_type, "saved to", overlaps_info_file, "\n", file = log_file, append = TRUE)
    
    # Initialize data.table to store results for this sample and type
    results_dt <- data.table(
      chr = as.character(seqnames(genome_type)),
      start = start(genome_type),
      end = end(genome_type),
      strand = as.character(strand(genome_type)),
      coverage = 0,
      numCs = 0,
      numTs = 0
    )
    
    # If overlaps are found, aggregate counts
    if (num_overlaps > 0) {
      # Extract overlapping indices
      region_idx <- queryHits(overlaps)
      methyl_idx <- subjectHits(overlaps)
      
      # Create a data.table with methylation data and corresponding region IDs
      overlap_dt <- data.table(
        region_id = region_idx,
        coverage = mcols(methyl_gr)$coverage[methyl_idx],
        numCs = mcols(methyl_gr)$numCs[methyl_idx],
        numTs = mcols(methyl_gr)$numTs[methyl_idx]
      )
      
      # Aggregate counts over regions
      agg_dt <- overlap_dt[, .(
        coverage = sum(coverage),
        numCs = sum(numCs),
        numTs = sum(numTs)
      ), by = region_id]
      
      # Update results_dt with aggregated counts
      results_dt[agg_dt$region_id, coverage := agg_dt$coverage]
      results_dt[agg_dt$region_id, numCs := agg_dt$numCs]
      results_dt[agg_dt$region_id, numTs := agg_dt$numTs]
    } else {
      cat("No overlaps found for sample", sample_id, "and type", current_type, ". All counts set to zero.\n", file = log_file, append = TRUE)
    }
    
    # Step 6: Save the Results to a File
    output_file <- file.path(export_dir, paste0(sample_id, "_", current_type, "_aggregated_methylation_counts.txt"))
    write.table(
      results_dt,
      file = output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    cat("Aggregated methylation counts for sample", sample_id, "and type", current_type, "have been saved to", output_file, "\n", file = log_file, append = TRUE)
  }
}

cat("Processing completed.\n", file = log_file, append = TRUE)

```