###----------------------------------------------------------###
### Devil testing ###
###----------------------------------------------------------###

#devtools::install_github("caravagnalab/devil")

library(devil)
library(tidyverse)

#Create design matrix
design_matrix <- model.matrix(~ cluster, data = metadata)

### Parameters inference ###

devil_fit_retina <- devil::fit_devil(input_matrix = rna_counts,
                             design_matrix = design_matrix,
                             overdispersion = TRUE,
                             offset=0,
                             size_factors=TRUE,
                             verbose=TRUE,
                             max_iter=500,
                             eps=1e-4,
                             parallel = TRUE)

### Statistical test ###
contrast_vector <- c(0,1)

stat_test_res <- devil::test_de(devil.fit = devil_fit_retina,
                                contrast = contrast_vector,
                                pval_adjust_method = "BH",
                                max_lfc = 10,
                                clusters = NULL)


###-----------------------------------------------------------------------------

whole_res <- lapply(unique(seurat_scRetina$seurat_clusters), function(c) {
  print(c)
  metadata$cluster <- seurat_obj$seurat_clusters
  metadata$group <- metadata$cluster == c
  
  design_matrix <- model.matrix(~group, metadata)
  
  counts <- as.matrix(counts)
  gg <- ((counts[,metadata$group] %>% rowSums()) == 0)
  bad_genes <- gg[gg == T] %>% names()
  
  gg <- ((counts[,!metadata$group] %>% rowSums()) == 0)
  bad_genes <- c(bad_genes, gg[gg == T] %>% names())
  
  fit <- devil::fit_devil(as.matrix(counts[!(rownames(counts) %in% bad_genes),]), 
                          esign_matrix, size_factors = T, verbose = T, parallel = T)
  
  res <- devil::test_de(fit, contrast = c(0,1))
  top_res <- res %>%
    dplyr::filter(lfc > .2, adj_pval <= 1e-50) %>%
    dplyr::arrange(-lfc) %>%
    dplyr::slice(1:40)
  top_res$cluster <- c
  top_res
}) %>% do.call("bind_rows", .)


###-----------------------------------------------------------------------------
### Cell clusters anotation
###-----------------------------------------------------------------------------
#install.packages("scCATCH")
#library(scCATCH)

scObj <- scCATCH::createscCATCH(as.matrix(counts), as.character(seurat_obj$seurat_clusters))

## find marker gene for each cluster
markers <- scCATCH::cellmatch %>%
  dplyr::filter(species == "Human", gene %in% whole_res$name)

obj <- scCATCH::findmarkergene(scObj, species = "Human", marker = markers, if_use_custom_marker = TRUE, tissue = "Retina")

# find cell type for each cluster
obj <- scCATCH::findcelltype(object = obj)
celltype_data <- obj@celltype

metadata$seurat_cluster <- seurat_obj$seurat_clusters

s <- 1
metadata$inferred_cell_type <- lapply(metadata$seurat_cluster, function(s) {
  if (any(celltype$cluster == s)) {
    return(celltype_data$cell_type[celltype$cluster == s])
  } else {
    return("null")
  }
  
}) %>% unlist()
mm <- metadata %>% as.tibble()
