# Define GO clusters
# Define GO clusters
GO_CLUSTERS <- list(
  `Defense and Immune Response` = c(
    "defense response",
    "immune system process",
    "response to chemical",
    "myeloid leukocyte migration",
    "leukocyte migration",
    "cell surface toll-like receptor signaling pathway",
    "biological process involved in interspecies interaction between organisms"
  ),
  `Gene Expression and Regulation` = c(
    "gene expression",
    "nucleic acid metabolic process",
    "protein metabolic process",
    "positive regulation of metabolic process",
    "regulation of biological process",
    "biological regulation"
  ),
  `Cellular Adhesion and Communication` = c(
    "homophilic cell adhesion via plasma membrane adhesion molecules",
    "cell-cell adhesion via plasma-membrane adhesion molecules",
    "regulation of substrate adhesion-dependent cell spreading",
    "cell-cell signaling"
  ),
  `Developmental and Morphogenetic Processes` = c(
    "multicellular organism development",
    "system development",
    "muscle cell development",
    "muscle system process",
    "positive regulation of cell differentiation",
    "cellular component assembly involved in morphogenesis",
    "non-membrane-bounded organelle assembly"
  ),
  `Muscle and Movement Processes` = c(
    "muscle contraction",
    "actin filament-based movement"
  ),
  `Cell Death and Cell Cycle` = c(
    "cell death",
    "programmed cell death",
    "positive regulation of cell cycle process",
    "regulation of growth"
  ),
  `Localization and Transport` = c(
    "establishment of localization",
    "transport",
    "cellular localization",
    "localization",
    "regulation of protein secretion"
  ),
  `Signaling and Regulation` = c(
    "intracellular signaling cassette",
    "regulation of signaling",
    "transforming growth factor beta receptor superfamily signaling pathway",
    "cognition"
  ),
  `Broad Categories` = c(
    "biological_process",
    "cellular process"
  )
)

remove_redundant_terms <- function(data, enrichment_col = "enrichmentScore", core_col = "core_enrichment", desc_col = "Description", threshold = 0.5) {
  filtered_data <- data %>%
    mutate(genes = strsplit(!!sym(core_col), "/"))
  n_terms <- nrow(filtered_data)
  overlap_matrix <- matrix(0, nrow = n_terms, ncol = n_terms,
                           dimnames = list(filtered_data[[desc_col]], filtered_data[[desc_col]]))
  for (i in 1:n_terms) {
    for (j in i:n_terms) {
      shared_genes <- length(intersect(filtered_data$genes[[i]], filtered_data$genes[[j]]))
      total_genes <- length(union(filtered_data$genes[[i]], filtered_data$genes[[j]]))
      jaccard_index <- shared_genes / total_genes

      # Fill overlap matrix with Jaccard index
      overlap_matrix[i, j] <- jaccard_index
      overlap_matrix[j, i] <- jaccard_index
    }
  }
  redundant_terms <- c()
  for (i in 1:(n_terms - 1)) {
    for (j in (i + 1):n_terms) {
      if (overlap_matrix[i, j] > threshold) {
        redundant_terms <- c(redundant_terms, filtered_data[[desc_col]][j])
      }
    }
  }
  non_redundant_data <- filtered_data %>%
    filter(!(!!sym(desc_col) %in% redundant_terms))

  redundant_data <- filtered_data %>%
    filter(!!sym(desc_col) %in% redundant_terms)

  return(list(
    non_redundant_data = non_redundant_data,
    redundant_data = redundant_data
  ))
}


get_simplified_GOterms = function(by=.05) {
  sdevil = lapply(seq(0, 1, by =by), function(c) {
    gseGO_devil_s = clusterProfiler::simplify(gseGO_devil, cutoff=c)
    dplyr::tibble(model = "devil", n_simplified = nrow(gseGO_devil_s@result), n_total = nrow(gseGO_devil@result)) %>%
      dplyr::mutate(f = n_simplified / n_total, c=c)
  }) %>% do.call("bind_rows", .)

  sglm = lapply(seq(0, 1, by =by), function(c) {
    gseGO_glm_s = clusterProfiler::simplify(gseGO_glm, cutoff=c)
    dplyr::tibble(model = "glmGamPoi", n_simplified = nrow(gseGO_glm_s@result), n_total = nrow(gseGO_glm@result)) %>%
      dplyr::mutate(f = n_simplified / n_total, c=c)
  }) %>% do.call("bind_rows", .)

  snebula = lapply(seq(0, 1, by = .05), function(c) {
    gseGO_nebula_s = clusterProfiler::simplify(gseGO_nebula, cutoff=c)
    dplyr::tibble(model = "NEBULA", n_simplified = nrow(gseGO_nebula_s@result), n_total = nrow(gseGO_nebula@result)) %>%
      dplyr::mutate(f = n_simplified / n_total, c=c)
  }) %>% do.call("bind_rows", .)

  dplyr::bind_rows(sdevil, sglm, snebula)
}

plot_dotplot_GO = function(devil_res, glm_res, nebula_res) {

  devil_res <- devil_res %>% dplyr::mutate(method = "devil", DE_type = ifelse(enrichmentScore > 0, "Up-regulated", "Down-regulated"))
  glm_res <- glm_res %>% dplyr::mutate(method = "glmGamPoi", DE_type = ifelse(enrichmentScore > 0, "Up-regulated", "Down-regulated"))
  nebula_res <- nebula_res %>% dplyr::mutate(method = "nebula", DE_type = ifelse(enrichmentScore > 0, "Up-regulated", "Down-regulated"))

  combined_data <- bind_rows(devil_res, glm_res, nebula_res) %>%
    dplyr::filter(!(Description %in% c("biological_process", "cellular process")))

  combined_data$GO_cluster = lapply(combined_data$Description, function(go_term) {
    for (cluster_name in names(GO_CLUSTERS)) {
      if (go_term %in% GO_CLUSTERS[[cluster_name]]) {
        return(cluster_name)
      }
    }
    print(go_term)
    return("GO term not found in any cluster")
  }) %>% unlist()

  combined_data$GO_cluster <- str_wrap(combined_data$GO_cluster, width = 20)

  # Enrichment plot

  plot_GO = combined_data %>%
    dplyr::mutate(
      Description = factor(
        Description,
        levels = combined_data %>%
          arrange(factor(method, levels = c("devil", "glmGamPoi", "nebula")), enrichmentScore) %>%
          pull(Description) %>%
          unique()
      ),
      DE_type = factor(DE_type, levels = c("Up-regulated", "Down-regulated")) # Ensure DE_type is ordered
    ) %>%
    ggplot(aes(x = method, y = Description, size = setSize, color = p.adjust)) +
    geom_point() +
    facet_grid(GO_cluster~DE_type, space = "free", scales = "free") +
    scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p-value") +
    theme_bw() +
    labs(title = "", x = "", y = "Biological Process GO term", size = "Gene Count") +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0.5)  # Rotate labels horizontally
    )
  plot_GO
}




enrichmentGO <- function(rna_deg_data) {
  rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)
  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID

  gseGO <- clusterProfiler::gseGO(
    genes,
    ont = "BP",  # use BP for dotplot, 'All' to perfom gene classification
    OrgDb = org.Hs.eg.db,
    minGSSize = 10,
    maxGSSize = 500,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    eps = 1e-10,
    by ="fgsea",
    seed = 1234
  )
  return(gseGO)
  #return(gseGO@result %>% as.data.frame())
}


# Function to compute the area under the curve (AUC) using the trapezoidal rule
auc <- function(x, y) {
  # Check if the lengths of x and y are the same
  if (length(x) != length(y)) {
    stop("Vectors x and y must have the same length")
  }

  # Sort the points based on x values (in case they are not ordered)
  ord <- order(x)
  x_ordered <- x[ord]
  y_ordered <- y[ord]

  # Compute the area using the trapezoidal rule
  n <- length(x_ordered)
  area <- sum((y_ordered[-1] + y_ordered[-n]) * diff(x_ordered)) / 2

  return(area)
}

