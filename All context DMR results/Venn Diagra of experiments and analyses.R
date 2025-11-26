# Import dataframes from analyses ----
# Build a function which will automatically detect installed packages and install required packages
check_packages_ClusterProfileR <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      print(paste(pkg, "is not installed, installing package..."))
      install.packages(pkg) # installs regular R packages
      BiocManager::install(pkg) # installs packages requiring BiocManager
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
      print(paste(pkg, "is installed and loaded"))
    }
  }
}

# Generate a dataframe with a list of packages required

pkg_list_ClusterProfileR<-c("BiocManager", # Have this first to ensure BiocManager packages will be installed
                            "tidyverse", 
                            "tibble", 
                            "stats",
                            "ggplot2", 
                            "dplyr", 
                            "tidyr", 
                            "ggfortify", 
                            "knitr",
                            "clusterProfiler",
                            "AnnotationForge",
                            "enrichplot",
                            "data.table",
                            "stats",
                            "plotly",
                            "ggVennDiagram",
                            "VennDiagram",
                            "limma",
                            "ggridges",
                            "pathview",
                            "KEGGREST",
                            "stringr",
                            "viridis",
                            "gghalves",
                            "ggside",
                            "flextable")

# Run the code
# NOTE: this code should be run twice, double checking is always good!

# Run Function::
check_packages_ClusterProfileR(pkg_list_ClusterProfileR)

# Setwd()
setwd("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/Methylation_Analyses/Methylation_Analyses/All context DMR results")

#--- Load methylation analyses ---
ORA_results_meth <- read.csv("ORA_results_combined_QUICKGO.csv", header = T, sep = ',')
GSEA_results_meth <- read.csv("GSEA_combined_data_named.csv", header = T, sep = ",")

#--- Load expression analyses ---
ORA_expr_exp1 <- read.csv("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/GSEA/Results/Iteration 2/ORA_exp1_results_QuickGO.csv") %>%
  select(-Description2) %>%
  mutate(Experiment = "Acute")

ORA_expr_exp2 <- read.csv("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/GSEA/Results/Iteration 2/ORA_exp2_results.csv") %>%
  mutate(Experiment = "Primed")

ORA_results_expr <- bind_rows(ORA_expr_exp1, ORA_expr_exp2)

GSEA_expr_exp1 <- read.csv("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/GSEA/Results/Iteration 2/GSEA_exp1_results_QUICKGO.csv") %>%
  select(-Description_QuickGO) %>%
  mutate(Experiment = "Acute")

GSEA_expr_exp2 <- read.csv("C:/Users/r02hw22/OneDrive - University of Aberdeen/Documents/GitHub/GSEA/Results/Iteration 2/GSEA_exp2_results_QUICKGO.csv") %>%
  select(-Description_QuickGo) %>%
  mutate(Experiment = "Primed")

GSEA_results_expr <- bind_rows(GSEA_expr_exp1, GSEA_expr_exp2)

#--- Process ORA Methylation ---
ORA_meth_venn_genes <- ORA_results_meth %>%
  select(Experiment = experiment, Gene_ID = geneID) %>%
  mutate(
    Experiment = recode(Experiment,
                        "Acute Experiment" = "Acute",
                        "Primed Experiment" = "Primed"),
    Analysis = "Methylation"
  ) %>%
  distinct()

#--- Process GSEA Methylation ---
GSEA_meth_venn_genes <- GSEA_results_meth %>%
  separate_rows(core_enrichment, sep = "/") %>%
  select(Experiment = exp, Gene_ID = core_enrichment) %>%
  mutate(
    Experiment = recode(Experiment,
                        "exp1" = "Acute",
                        "exp2" = "Primed"),
    Analysis = "Methylation"
  ) %>%
  distinct()

#--- Process ORA Expression ---
ORA_expr_venn_genes <- ORA_results_expr %>%
  select(Experiment, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(Gene_ID = geneID) %>%
  mutate(Analysis = "Expression") %>%
  distinct()

#--- Process GSEA Expression ---
GSEA_expr_venn_genes <- GSEA_results_expr %>%
  separate_rows(core_enrichment, sep = "/") %>%
  select(Experiment, Gene_ID = core_enrichment) %>%
  mutate(Analysis = "Expression") %>%
  distinct()

#--- Combine all sets together ---
venn_genes_df <- bind_rows(
  ORA_meth_venn_genes,
  GSEA_meth_venn_genes,
  ORA_expr_venn_genes,
  GSEA_expr_venn_genes
) %>% distinct()


gene_sets <- venn_genes_df %>%
  mutate(Group = paste(Experiment, Analysis, sep = " ")) %>%
  group_by(Group) %>%
  summarise(Genes = list(unique(Gene_ID))) %>%
  deframe()

str(gene_sets)

library(purrr)
library(dplyr)
library(tibble)

# Generate all pairwise combinations
pairwise_comparisons <- combn(names(gene_sets), 2, simplify = FALSE)

# Build a tibble showing group names and the genes they share
intersections_named <- map_df(pairwise_comparisons, function(pair) {
  common_genes <- intersect(gene_sets[[pair[1]]], gene_sets[[pair[2]]])
  
  if (length(common_genes) > 0) {
    tibble(
      Groups = paste(pair, collapse = " & "),
      Shared_Genes = list(common_genes),
      n_shared = length(common_genes)
    )
  } else {
    NULL
  }
})

# View results
intersections_named[[2]]


#install.packages("ggvenn")
library(ggvenn)

ggvenn(
  gene_sets,
  show_elements = FALSE,
  stroke_size = 0.5,
  fill_alpha = 0.4,
  set_name_size = 4
)


library(VennDiagram)

venn.plot <- venn.diagram(
  gene_sets,
  category.names = names(gene_sets),
  filename = NULL,
  fill = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
  alpha = 0.4,
  cex = 1.2,
  cat.cex = 1.1,
  margin = 0.05
)

grid::grid.draw(venn.plot)

# Open PNG device
png("venn_gene_plot.png", width = 2000, height = 1200, res = 300)

# Draw the Venn plot
grid::grid.draw(venn.plot)

# Close the device
dev.off()


# Now with the GO IDs
#--- Process ORA Methylation ---
ORA_meth_venn_go <- ORA_results_meth %>%
  select(Experiment = experiment, GO_ID = ID) %>%
  mutate(
    Experiment = recode(Experiment,
                        "Acute Experiment" = "Acute",
                        "Primed Experiment" = "Primed"),
    Analysis = "Methylation"
  ) %>%
  distinct()

#--- Process GSEA Methylation ---
GSEA_meth_venn_go <- GSEA_results_meth %>%
  separate_rows(core_enrichment, sep = "/") %>%
  select(Experiment = exp, GO_ID = ID) %>%
  mutate(
    Experiment = recode(Experiment,
                        "exp1" = "Acute",
                        "exp2" = "Primed"),
    Analysis = "Methylation"
  ) %>%
  distinct()

#--- Process ORA Expression ---
ORA_expr_venn_go <- ORA_results_expr %>%
  select(Experiment, ID) %>%
  separate_rows(ID, sep = "/") %>%
  rename(GO_ID = ID) %>%
  mutate(Analysis = "Expression") %>%
  distinct()

#--- Process GSEA Expression ---
GSEA_expr_venn_go <- GSEA_results_expr %>%
  separate_rows(core_enrichment, sep = "/") %>%
  select(Experiment, GO_ID = ID) %>%
  mutate(Analysis = "Expression") %>%
  distinct()

#--- Combine all sets together ---
venn_go_df <- bind_rows(
  ORA_meth_venn_go,
  GSEA_meth_venn_go,
  ORA_expr_venn_go,
  GSEA_expr_venn_go
) %>% distinct()

# plot together ----
go_sets <- venn_go_df %>%
  mutate(Group = paste(Experiment, Analysis, sep = " ")) %>%
  group_by(Group) %>%
  summarise(Genes = list(unique(GO_ID))) %>%
  deframe()


install.packages("combinat")
library(combinat)

# Get all group names
group_names <- names(go_sets)

# Create all combinations of group names, for sizes 2 to N
all_combinations <- map(2:length(group_names), function(k) {
  combn(group_names, k, simplify = FALSE)
}) %>% flatten()

# Calculate intersection for each combination
all_intersections <- map_df(all_combinations, function(combo) {
  intersected_genes <- purrr::reduce(go_sets[combo], intersect)
  
  tibble(
    Groups = paste(combo, collapse = " & "),
    Group_List = list(combo),
    Shared_Genes = list(intersected_genes),
    n_shared = length(intersected_genes)
  )
}) %>%
  filter(n_shared > 0)  # Optional: remove empty intersections

# View the results
all_intersections[[3]]

ggvenn(
  go_sets,
  show_elements = FALSE,
  stroke_size = 0.5,
  fill_alpha = 0.4,
  set_name_size = 4
)

venn.plot <- venn.diagram(
  go_sets,
  category.names = names(go_sets),
  filename = NULL,
  fill = c("red", "blue", "green", "yellow"),
  alpha = 0.4,
  cex = 1.2,
  cat.cex = 1.1,
  margin = 0.05
)



# Open PNG device
png("venn_gos_plot.png", width = 2000, height = 1200, res = 300)

# Draw the Venn plot
grid::grid.draw(venn.plot)

# Close the device
dev.off()

library(purrr)
# Making a table of the results
# Step 1: Create named GO term sets
go_sets <- venn_go_df %>%
  mutate(Group = paste(Experiment, Analysis, sep = " ")) %>%
  group_by(Group) %>%
  summarise(GO_Terms = list(unique(GO_ID)), .groups = "drop") %>%
  deframe()

# Step 2: All pairwise combinations (excluding self-pairs and duplicates)
set_names <- names(go_sets)
pairwise_combos <- expand.grid(SetA = set_names, SetB = set_names, stringsAsFactors = FALSE) %>%
  filter(SetA < SetB)

# Step 3: Add shared GO term details
intersection_summary <- pairwise_combos %>%
  rowwise() %>%
  mutate(
    Shared_GOs = list(intersect(go_sets[[SetA]], go_sets[[SetB]])),
    Shared_Count = length(Shared_GOs),
    Shared_GO_IDs = paste(Shared_GOs, collapse = ", ")
  ) %>%
  ungroup() %>%
  select(SetA, SetB, Shared_Count, Shared_GO_IDs)

# Step 4: Format with flextable
ft <- flextable(intersection_summary) %>%
  set_header_labels(
    SetA = "Set A",
    SetB = "Set B",
    Shared_Count = "Shared GO Terms",
    Shared_GO_IDs = "GO IDs"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  fontsize(part = "header", size = 10) %>%
  fontsize(part = "body", size = 8) %>%
  padding(part = "all", padding.top = 2, padding.bottom = 2) %>%
  set_caption("Table: Shared GO terms between experimental and analytical groups")

# View the table
ft
