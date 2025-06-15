# =============================================================================
# PARAMETERS - All hard-coded parameters at the beginning
# =============================================================================

rm(list=ls())
UTILS_DIR <- "utils/"
DATA_DIR <- "../de_analysis/nullpower/datasets/"
FDR <- 0.05  # False Discovery Rate threshold
GENE_FILTER_THRESHOLD <- 0.05  # Minimum mean expression for gene filtering
EXPRESSION_QUANTILES <- c(0.333, 0.666)  # Quantiles for expression groups
RANDOM_SEED <- 12345  # Seed for reproducibility
N_PSEUDO_DONORS <- 8  # Number of pseudo-donors for pseudo-replicates experiment
N_RANDOMIZATIONS <- 2  # Number of times to repeat each experiment for each dataset

# =============================================================================
# SETUP
# =============================================================================

# Source utility functions
source(file.path(UTILS_DIR, "utils.R"))
source(file.path(UTILS_DIR, "edgeR.R"))
source(file.path(UTILS_DIR, "limma.R"))
source(file.path(UTILS_DIR, "glmGamPoi.R"))
source(file.path(UTILS_DIR, "nebula.R"))
source(file.path(UTILS_DIR, "devil.R"))

# Load libraries
library(ggplot2)
library(muscData)
library(scRNAseq)

# Define analysis methods
list.func <- list(
  glmgp.mult,
  edger.mult,
  limma.mult,
  glmgp.cell.mult,
  nebula.mult,
  devil.base,
  devil.mixed
)

method_names <- c('glmGamPoi (Pb)', 'edgeR (Pb)', 'limma (Pb)', 
                  'glmGamPoi (cell)', 'Nebula', 'Devil (base)', 'Devil (mixed)')

# Define datasets to analyze
dataset_names <- c("Crowell19", "BacherTCellData", "ReprocessedFluidigmData", "HuCortexData")
dataset_names = c("Crowell19", "BacherTCellData", "HuCortexData", "LedergorMyelomaData", "GrunHSCData")
# =============================================================================
# EXPERIMENT FUNCTIONS
# =============================================================================

get_dataset = function(dataset_name) {
  print(paste("Loading dataset:", dataset_name))
  
  if (dataset_name == "Crowell19") {
    data <- muscData::Crowell19_4vs4(metadata = TRUE)
    sce <- data[["EH3297"]]
    sce$donor_id <- sce$sample_id
    sce.sub <- sce[, sce$group_id == "Vehicle"]
    cnt_mat = sce.sub@assays@data$counts
    meta = sce.sub@colData
  } else if (dataset_name == "BacherTCellData") {
    sce = scRNAseq::BacherTCellData()
    sce$donor_id <- sce$sample_id
    sce.sub <- sce[, sce$diagnosis == "Healthy"]
    cnt_mat = sce.sub@assays@data$counts
    meta = sce.sub@colData
    all(colnames(cnt_mat) == rownames(meta))
  } else if (dataset_name == "HuCortexData") {
    sce = scRNAseq::HuCortexData(mode = c("ctx"))
    sce$donor_id <- sce$Sample
    sce.sub <- sce[, sce$Mode == "ctx"]
    cnt_mat = sce.sub@assays@data$counts
    meta = sce.sub@colData
    all(colnames(cnt_mat) == rownames(meta))
  } else if (dataset_name == "LedergorMyelomaData") {
    # LedergorMyelomaData
    sce = scRNAseq::LedergorMyelomaData()
    sce$donor_id <- sce$Subject_ID
    sce.sub <- sce[, sce$Condition == "Control"]
    cnt_mat = sce.sub@assays@data$counts
    meta = sce.sub@colData
    all(colnames(cnt_mat) == rownames(meta))
  } else if (dataset_name == "GrunHSCData") {
    # GrunHSCData
    sce = scRNAseq::GrunHSCData()
    sce$donor_id <- sce$sample
    sce.sub <- sce[, sce$protocol == "sorted hematopoietic stem cells"]
    cnt_mat = sce.sub@assays@data$counts
    meta = sce.sub@colData
    all(colnames(cnt_mat) == rownames(meta))
  } else {
    stop(paste("Unknown dataset:", dataset_name))
  } 
  print(paste("Dataset loaded -", nrow(cnt_mat), "genes x", ncol(cnt_mat), "cells"))
  list(cnt_mat = as.matrix(cnt_mat), meta = meta)
}

apply_qc_filtering <- function(cnt_matrix, col_data) {
  print("Applying QC filtering...")
  
  # Calculate QC metrics
  library_sizes <- colSums(cnt_matrix)
  genes_per_cell <- colSums(cnt_matrix > 0)
  cells_per_gene <- rowSums(cnt_matrix > 0)
  
  # Calculate mitochondrial gene percentage
  mito_genes <- grep("^MT-|^mt-", rownames(cnt_matrix), value = TRUE)
  if (length(mito_genes) == 0) {
    mito_genes <- grep("^MTRNR|^mtrnr", rownames(cnt_matrix), value = TRUE)
  }
  
  if (length(mito_genes) > 0) {
    mito_counts <- colSums(cnt_matrix[mito_genes, , drop = FALSE])
    mito_percent <- (mito_counts / library_sizes) * 100
    print(paste("Found", length(mito_genes), "mitochondrial genes"))
  } else {
    mito_percent <- rep(0, ncol(cnt_matrix))
    print("Warning: No mitochondrial genes found with standard naming conventions")
  }
  
  # Apply QC thresholds
  cell_lib_filter <- library_sizes >= 1000
  cell_mito_filter <- mito_percent <= 15
  cell_gene_filter <- genes_per_cell >= 500
  
  # Updated gene filter: expressed in >= 5% of cells
  min_cells <- ceiling(0.05 * ncol(cnt_matrix))  # 5% of total cells
  gene_cell_filter <- cells_per_gene >= min_cells
  
  # Combine filters
  cells_to_keep <- cell_lib_filter & cell_mito_filter & cell_gene_filter
  genes_to_keep <- gene_cell_filter
  
  print("=== QC FILTERING SUMMARY ===")
  print(paste("Library size filter (>=1000):", sum(cell_lib_filter), "/", length(cell_lib_filter), "cells pass"))
  print(paste("Mitochondrial filter (<=15%):", sum(cell_mito_filter), "/", length(cell_mito_filter), "cells pass"))
  print(paste("Genes per cell filter (>=500):", sum(cell_gene_filter), "/", length(cell_gene_filter), "cells pass"))
  print(paste("Combined cell filters:", sum(cells_to_keep), "/", length(cells_to_keep), "cells pass"))
  print(paste("Gene filter (>=5% of cells):", sum(genes_to_keep), "/", length(genes_to_keep), "genes pass"))
  
  # Apply filters
  cnt_matrix_filtered <- cnt_matrix[genes_to_keep, cells_to_keep]
  col_data_filtered <- col_data[cells_to_keep, ]
  
  print(paste("After QC filtering:", nrow(cnt_matrix_filtered), "genes x", ncol(cnt_matrix_filtered), "cells"))
  
  return(list(cnt_matrix = cnt_matrix_filtered, col_data = col_data_filtered))
}

# Function for biological replicates experiment
run_bio_replicates_experiment <- function(cnt_matrix, col_data, dataset_name, randomization_id) {
  # Get unique donors/samples
  unique_donors <- unique(col_data$donor_id)
  
  if (length(unique_donors) < 2) {
    print("Warning: Less than 2 donors available for biological replicates")
    return(NULL)
  }
  
  # Randomly assign donors to groups (use dataset and randomization for unique seed)
  set.seed(RANDOM_SEED + match(dataset_name, dataset_names) * 100 + randomization_id)
  donor_assignment <- data.frame(
    donor_id = unique_donors,
    group = sample(rep(c(0, 1), length.out = length(unique_donors)))
  )
  
  print(paste("Biological replicates - Randomization", randomization_id, ":"))
  print(paste("  Group 0:", sum(donor_assignment$group == 0), "donors"))
  print(paste("  Group 1:", sum(donor_assignment$group == 1), "donors"))
  
  # Create modified metadata
  modified_coldata <- col_data
  modified_coldata$tx_cell <- donor_assignment$group[
    match(modified_coldata$donor_id, donor_assignment$donor_id)
  ]
  modified_coldata$id <- modified_coldata$donor_id
  
  # Order cells by donor to group them together (same as pseudo-replicates function)
  cell_order <- order(modified_coldata$donor_id)
  
  # Reorder count matrix columns and metadata rows
  cnt_matrix_ordered <- cnt_matrix[, cell_order]
  modified_coldata_ordered <- modified_coldata[cell_order, ]
  
  return(list(cnt_matrix = cnt_matrix_ordered, col_data = modified_coldata_ordered))
}

# Function for pseudo-replicates experiment
run_pseudo_replicates_experiment <- function(cnt_matrix, col_data, dataset_name, randomization_id) {
  # Set seed for reproducibility (use dataset and randomization for unique seed)
  set.seed(RANDOM_SEED + match(dataset_name, dataset_names) * 1000 + randomization_id)
  
  # Randomly assign cells to pseudo-donors
  n_cells <- ncol(cnt_matrix)
  pseudo_donor_ids <- paste0("pseudo_donor_", 1:N_PSEUDO_DONORS)
  
  # Randomly assign cells to pseudo-donors
  cell_to_pseudo_donor <- sample(pseudo_donor_ids, n_cells, replace = TRUE)
  
  # Randomly assign pseudo-donors to groups
  pseudo_donor_assignment <- data.frame(
    donor_id = pseudo_donor_ids,
    group = sample(rep(c(0, 1), length.out = N_PSEUDO_DONORS))
  )
  
  print(paste("Pseudo-replicates - Randomization", randomization_id, ":"))
  print(paste("  Group 0:", sum(pseudo_donor_assignment$group == 0), "pseudo-donors"))
  print(paste("  Group 1:", sum(pseudo_donor_assignment$group == 1), "pseudo-donors"))
  
  # Create modified metadata
  modified_coldata <- col_data
  modified_coldata$sample_id <- cell_to_pseudo_donor
  modified_coldata$donor_id <- cell_to_pseudo_donor
  modified_coldata$id <- cell_to_pseudo_donor
  
  # Assign group based on pseudo-donor assignment
  modified_coldata$tx_cell <- pseudo_donor_assignment$group[
    match(modified_coldata$donor_id, pseudo_donor_assignment$donor_id)
  ]
  
  # Order cells by pseudo-donor to group them together
  cell_order <- order(modified_coldata$donor_id)
  
  # Reorder count matrix columns and metadata rows
  cnt_matrix_ordered <- cnt_matrix[, cell_order]
  modified_coldata_ordered <- modified_coldata[cell_order, ]
  
  return(list(cnt_matrix = cnt_matrix_ordered, col_data = modified_coldata_ordered))
}

# Function to run analysis methods
run_analysis_methods <- function(cnt_matrix, col_data) {
  list.result.method <- list()
  timings <- c()
  
  for (int.test in 1:length(list.func)) {
    print(paste("Running method", int.test, "-", method_names[int.test]))
    tryCatch({
      result <- list.func[[int.test]](cnt_matrix, col_data)[, 4:5]
      timings <- c(timings, unique(result[, 2]))
      list.result.method[[int.test]] <- result[, 1]
    }, error = function(e) {
      print(paste("Error in method", int.test, ":", e$message))
      list.result.method[[int.test]] <- rep(NA, nrow(cnt_matrix))
    })
  }
  
  df.result <- do.call(rbind, list.result.method)
  df.result <- df.result[1:length(method_names), ]
  rownames(df.result) <- method_names
  
  return(df.result)
}

# Function to process results
process_results <- function(df_result, dataset_name, experiment_type, randomization_id, mean_exp) {
  # Calculate expression quantiles for this specific dataset
  expression_quantiles <- quantile(mean_exp, EXPRESSION_QUANTILES, na.rm = TRUE)
  
  df_res <- lapply(1:nrow(df_result), function(j) {
    method <- rownames(df_result)[j]
    p_values <- df_result[j, ]
    
    # Handle NA values
    valid_p_values <- !is.na(p_values)
    if(sum(valid_p_values) == 0) {
      return(data.frame(
        dataset = dataset_name,
        experiment_type = experiment_type,
        randomization_id = randomization_id,
        method = method,
        expression_group = c("low", "medium", "high"),
        n_genes = c(0, 0, 0),
        n_deg = c(0, 0, 0),
        prop_deg = c(0, 0, 0)
      ))
    }
    
    adjusted_p <- rep(NA, length(p_values))
    adjusted_p[valid_p_values] <- p.adjust(p_values[valid_p_values], method = "BH")
    
    # Create expression groups
    low_exp_genes <- mean_exp <= expression_quantiles[1]
    high_exp_genes <- mean_exp >= expression_quantiles[2]
    medium_exp_genes <- !low_exp_genes & !high_exp_genes
    
    # Count DEGs in each group
    n_deg_low <- sum(adjusted_p[low_exp_genes] <= FDR, na.rm = TRUE)
    n_deg_medium <- sum(adjusted_p[medium_exp_genes] <= FDR, na.rm = TRUE)
    n_deg_high <- sum(adjusted_p[high_exp_genes] <= FDR, na.rm = TRUE)
    
    # Create result data frame
    data.frame(
      dataset = dataset_name,
      experiment_type = experiment_type,
      randomization_id = randomization_id,
      method = method,
      expression_group = c("low", "medium", "high"),
      n_genes = c(sum(low_exp_genes), sum(medium_exp_genes), sum(high_exp_genes)),
      n_deg = c(n_deg_low, n_deg_medium, n_deg_high),
      prop_deg = c(n_deg_low, n_deg_medium, n_deg_high) / c(sum(low_exp_genes), sum(medium_exp_genes), sum(high_exp_genes))
    )
  })
  
  return(do.call("rbind", df_res))
}

# =============================================================================
# MAIN ANALYSIS LOOP
# =============================================================================
print("Starting multi-dataset analysis...")
start_time <- Sys.time()

# Initialize results list
all_results <- list()
result_counter <- 1

# Run analysis for each dataset
for (dataset_name in dataset_names) {
  print(paste("=== Processing dataset:", dataset_name, "==="))
  
  tryCatch({
    # Load dataset
    dataset <- get_dataset(dataset_name)
    cnt_matrix <- dataset$cnt_mat
    col_data <- dataset$meta
    
    # Apply QC filtering
    filtered_data <- apply_qc_filtering(cnt_matrix, col_data)
    cnt_matrix <- filtered_data$cnt_matrix
    col_data <- filtered_data$col_data
    
    # Filter genes by mean expression
    mean_exp <- rowMeans(cnt_matrix)
    gene_filter <- mean_exp > GENE_FILTER_THRESHOLD
    cnt_matrix <- cnt_matrix[gene_filter, ]
    mean_exp <- mean_exp[gene_filter]
    
    # Check if we have enough genes and cells
    if (nrow(cnt_matrix) < 100) {
      print(paste("Skipping", dataset_name, "- too few genes after filtering:", nrow(cnt_matrix)))
      next
    }
    
    if (ncol(cnt_matrix) < 10) {
      print(paste("Skipping", dataset_name, "- too few cells:", ncol(cnt_matrix)))
      next
    }
    
    print(paste("Dataset", dataset_name, "after filtering - Cells:", ncol(cnt_matrix), "Genes:", nrow(cnt_matrix)))
    
    # Run N_RANDOMIZATIONS for both biological and pseudo-replicates experiments
    for (randomization_id in 1:N_RANDOMIZATIONS) {
      print(paste("--- Running randomization", randomization_id, "of", N_RANDOMIZATIONS, "for", dataset_name, "---"))
      
      # Run biological replicates experiment
      print("Running biological replicates experiment...")
      bio_data <- run_bio_replicates_experiment(cnt_matrix, col_data, dataset_name, randomization_id)
      
      if (!is.null(bio_data)) {
        bio_results <- run_analysis_methods(bio_data$cnt_matrix, bio_data$col_data)
        bio_processed <- process_results(bio_results, dataset_name, "bio_replicates", randomization_id, mean_exp)
        all_results[[result_counter]] <- bio_processed
        result_counter <- result_counter + 1
      }
      
      # Run pseudo-replicates experiment
      print("Running pseudo-replicates experiment...")
      pseudo_data <- run_pseudo_replicates_experiment(cnt_matrix, col_data, dataset_name, randomization_id)
      pseudo_results <- run_analysis_methods(pseudo_data$cnt_matrix, pseudo_data$col_data)
      pseudo_processed <- process_results(pseudo_results, dataset_name, "pseudo_replicates", randomization_id, mean_exp)
      all_results[[result_counter]] <- pseudo_processed
      result_counter <- result_counter + 1
      
      print(paste("Completed randomization", randomization_id, "for", dataset_name))
    }
    
    print(paste("Completed all randomizations for", dataset_name))
    
  }, error = function(e) {
    print(paste("Error processing dataset", dataset_name, ":", e$message))
  })
}

end_time <- Sys.time()

# Combine all results
if (length(all_results) > 0) {
  final_results <- do.call("rbind", all_results)
  saveRDS(final_results, "multi_dataset_results.rds")
  
  # Print summary
  print("=== ANALYSIS SUMMARY ===")
  print(paste("Successfully analyzed", length(unique(final_results$dataset)), "datasets"))
  print(paste("Randomizations per dataset:", N_RANDOMIZATIONS))
  print(paste("Experiment types: biological replicates, pseudo-replicates"))
  print("Datasets processed:")
  for(ds in unique(final_results$dataset)) {
    bio_count <- sum(final_results$dataset == ds & final_results$experiment_type == "bio_replicates")
    pseudo_count <- sum(final_results$dataset == ds & final_results$experiment_type == "pseudo_replicates")
    print(paste("  -", ds, "- Bio:", bio_count / (length(method_names) * 3), "Pseudo:", pseudo_count / (length(method_names) * 3)))
  }
  
  print(paste("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
  print("Results saved to: multi_dataset_results.rds")
} else {
  print("No datasets were successfully processed!")
}

stop()

final_results = readRDS("celltype_results.rds")

# =============================================================================
# VISUALIZATION
# =============================================================================

# Calculate method ordering based on mean DEGs across all cell types
method_order <- final_results %>%
  group_by(method) %>%
  summarise(mean_deg = mean(n_deg, na.rm = TRUE)) %>%
  arrange(mean_deg) %>%
  pull(method)

# Create comprehensive plot
p1 <- final_results %>% 
  mutate(
    method = factor(method, levels = method_order),
    expression_group = factor(expression_group, levels = c("low", "medium", "high")),
    seq = ifelse(grepl("Pb", method), "Pseudo-bulk", "Single-cell")
  ) %>% 
  ggplot(aes(x = method, y = n_deg, fill = expression_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(cell_type ~ experiment_type, scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Differential Expression Analysis Results by Cell Type",
    subtitle = "Comparison of Bio-replicates vs Pseudo-replicates",
    x = "Method",
    y = "Number of DEGs",
    fill = "Expression Group"
  ) +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p1)

# Alternative plot - focus on differences between experiment types
p2 <- final_results %>% 
  select(cell_type, experiment_type, method, expression_group, n_deg) %>%
  pivot_wider(names_from = experiment_type, values_from = n_deg) %>%
  mutate(
    deg_difference = pseudo_replicates - bio_replicates,
    method = factor(method, levels = method_order)
  ) %>%
  ggplot(aes(x = method, y = deg_difference, fill = expression_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cell_type, scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Difference in DEG Detection: Pseudo-replicates - Bio-replicates",
    subtitle = "Positive values indicate more DEGs detected with pseudo-replicates",
    x = "Method",
    y = "Difference in Number of DEGs",
    fill = "Expression Group"
  ) +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p2)

# Summary statistics
print("=== SUMMARY STATISTICS ===")
summary_stats <- final_results %>%
  group_by(cell_type, experiment_type, method, expression_group) %>%
  summarise(
    n_deg = first(n_deg),
    prop_deg = first(prop_deg),
    .groups = "drop"
  )

print(summary_stats)

# Summary by cell type
print("=== SUMMARY BY CELL TYPE ===")
celltype_summary <- final_results %>%
  group_by(cell_type, experiment_type) %>%
  summarise(
    total_deg = sum(n_deg),
    mean_prop_deg = mean(prop_deg),
    .groups = "drop"
  )

print(celltype_summary)

# Save results (optional)
# write_csv(final_results, "de_analysis_celltype_results.csv")
# write_csv(summary_stats, "de_analysis_celltype_summary.csv")
