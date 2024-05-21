
DATASET_NAMES <- c("BaronPancreasData", 'ZhaoImmuneLiverDataBlood', 'DarmanisBrainData', 'ZhaoImmuneLiverDataLiver',
                   "BigBloodData", "BigLiverData")

read_data <- function(dataset_name, data_path=NULL) {
  if (dataset_name == "BaronPancreasData") {
    data <- scRNAseq::BaronPancreasData(which = "human")
    metadata <- data@colData
    metadata$cell_type <- metadata$label
    tissue <- "pancreas"
    counts <- data@assays@data[[1]]
  } else if (dataset_name == "DarmanisBrainData") {
    data <- scRNAseq::DarmanisBrainData()
    metadata <- data@colData
    metadata$cell_type <- metadata$cell.type
    metadata$donor <- metadata$experiment_sample_name
    tissue <- "brain"
    counts <- data@assays@data[[1]]
  } else if (dataset_name == 'ZhaoImmuneLiverDataBlood') {
    data <- scRNAseq::ZhaoImmuneLiverData()

    filter_cell <- grepl('blood', data$sample)
    metadata <- dplyr::tibble(
      cell_type = data$fine[filter_cell],
      donor = data$sample[filter_cell]
    )
    counts <- data@assays@data[[1]]
    counts <- counts[,filter_cell]
    tissue <- "blood"
  } else if (dataset_name == 'ZhaoImmuneLiverDataLiver') {
    data <- scRNAseq::ZhaoImmuneLiverData()

    filter_cell <- grepl('liver', data$sample)
    metadata <- dplyr::tibble(
      cell_type = data$fine[filter_cell],
      donor = data$sample[filter_cell]
    )
    counts <- data@assays@data[[1]]
    counts <- counts[,filter_cell]
    tissue <- "liver"
  } else if (dataset_name == 'BigBloodData') {
    seurat_data <- readRDS(data_path)
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    metadata <- dplyr::rename(metadata,c("cell_type" = "celltype.l2", "cell_type2" = "cell_type", "donor" = "donor_id"))
    tissue = "blood"
  } else if (dataset_name == 'BigLiverData') {
    seurat_data <- readRDS(data_path)
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    meta_features <- seurat_data@assays[["RNA"]]@meta.features
    rownames(counts) <- meta_features$feature_name
    metadata <- dplyr::rename(metadata,c("donor" = "donor_id"))
    tissue = "liver"
  }  else {
    stop("Dataset name not recognized")
  }

  return(list(counts=counts, metadata=metadata, tissue=tissue))
}

prep_seurat_object <- function(input_data, NPC=50, cluster_res=1) {
  counts <- input_data$counts
  metadata <- input_data$metadata
  metadata$cell_type %>% unique()
  metadata$donor %>% unique()

  # Filter cells
  total_counts <- colSums(counts)
  total_features <- colSums(counts > 0)

  mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
  feat100_filter <- total_features < 100
  feat_mad_filter <- total_features > 5 * mad(total_features)

  mitocondrial_genes <- grepl("^MT-", rownames(counts))
  mitocondiral_prop <- colSums(counts[mitocondrial_genes, ]) / colSums(counts)
  mit_prop_filter <- mitocondiral_prop > .1

  # ribosomal_genes <- grepl("^RPS", rownames(counts)) | grepl("^RPL", rownames(counts))
  # ribosomal_prop <- colSums(counts[ribosomal_genes, ]) / colSums(counts)
  # rib_prop_filter <- ribosomal_prop < .1

  #cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter | rib_prop_filter | mit_prop_filter
  cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter |  mit_prop_filter

  counts <- counts[, !cell_outliers_filter]
  metadata <- metadata[!cell_outliers_filter, ]
  rownames(metadata) <- colnames(counts)

  # Filter genes
  non_expressed_genes <- rowMeans(counts) <= 0.01
  counts <- counts[!non_expressed_genes, ]

  seurat_obj = SeuratObject::CreateSeuratObject(counts = counts, meta.data = as.tibble(metadata))
  seurat_obj <- NormalizeData(seurat_obj) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = NPC) %>%
    FindNeighbors(dims = 1:NPC) %>%
    FindClusters(resolution = cluster_res) %>%
    RunUMAP(dims=1:NPC)

  new_clusters <- as.numeric(seurat_obj$seurat_clusters)
  new_clusters <- factor(new_clusters, levels = sort(unique(new_clusters)))
  seurat_obj$seurat_clusters <- new_clusters

  seurat_obj
}

perform_analysis <- function(seurat_obj, method = "devil") {
  if (!(method %in% c('devil', "glmGamPoi", 'nebula'))) {stop('method not recognized')}

  if (method == 'devil') {
    c <- 1
    whole_res <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
      print(c)

      design_matrix <- model.matrix(~group, dplyr::tibble(group = seurat_obj$seurat_clusters == c))
      counts <- as.matrix(seurat_obj@assays$RNA$counts)

      gg <- ((counts[,seurat_obj$seurat_clusters == c] %>% rowSums()) == 0)
      bad_genes <- gg[gg == T] %>% names()
      counts <- counts[!(rownames(counts) %in% bad_genes),]

      fit <- devil::fit_devil(counts, design_matrix, verbose = T, size_factors = T, parallel.cores = 1)

      clusters <- as.numeric(as.factor(seurat_obj$donor))
      res <- devil::test_de(fit, contrast = c(0,1), clusters = clusters, max_lfc = Inf)
      res %>% dplyr::mutate(cluster = c)

    }) %>% do.call("bind_rows", .)
  
  } else if (method == "glmGamPoi") {
    c <- 1
    whole_res <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
      print(c)

      design_matrix <- model.matrix(~group, dplyr::tibble(group = seurat_obj$seurat_clusters == c))
      counts <- as.matrix(seurat_obj@assays$RNA$counts)

      gg <- ((counts[,seurat_obj$seurat_clusters == c] %>% rowSums()) == 0)
      bad_genes <- gg[gg == T] %>% names()
      counts <- counts[!(rownames(counts) %in% bad_genes),]

      fit <- glmGamPoi::glm_gp(counts, design_matrix, size_factors = T, verbose = T)
      res <- glmGamPoi::test_de(fit, contrast = c(0,1))
      res <- res %>% select(name, pval, adj_pval, lfc)

      res %>% dplyr::mutate(cluster = c)

    }) %>% do.call("bind_rows", .)
  } else if (method == 'nebula') {
    c <- 1
    whole_res <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
      print(c)
      
      design_matrix <- model.matrix(~group, dplyr::tibble(group = seurat_obj$seurat_clusters == c))
      counts <- as.matrix(seurat_obj@assays$RNA$counts)
    
      
      gg <- ((counts[,seurat_obj$seurat_clusters == c] %>% rowSums()) == 0)
      bad_genes <- gg[gg == T] %>% names()
      counts <- counts[!(rownames(counts) %in% bad_genes),]
      
      clusters <- as.numeric(as.factor(metadata$donor))
      data_g = group_cell(count=counts,id=clusters,pred=design_matrix)
      
      fit <- nebula::nebula(data_g$counts, id = data_g$clusters, pred = data_g$design_matrix, ncore = 1)
      
      res <- dplyr::tibble(
        name = fit$summary$gene,
        pval = fit$summary$p_groupTRUE,
        adj_pval = p.adjust(fit$summary$p_groupTRUE, "BH"),
        lfc=fit$summary$logFC_groupTRUE
      )
      
      res %>% dplyr::mutate(cluster = c)
    }) %>% do.call("bind_rows", .)
  }
  whole_res
}

prepScMayoInput <- function(de_res, count_matrix, seurat_clusters, n_markers=50, lfc_cut=1, pval_cut=1e-10, distinct_marker=TRUE) {
  if (!(all(colnames(de_res) == c("name", "pval", "adj_pval", "lfc", "cluster")))) {stop('input de_res have wrong column names')}

  cluster_values <- de_res$cluster %>% unique()
  remove_genes <- grepl("^ENS", de_res$name)
  de_res <- de_res[!remove_genes, ]

  scMayoInput <- dplyr::tibble()
  for (c in cluster_values) {
    top_res <- de_res %>%
      dplyr::filter(cluster == c, lfc > lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::arrange(-lfc) %>%
      dplyr::slice(1:n_markers)

    names <- rownames(count_matrix)
    if (all(grepl("ENSG", names))) {
      ensembl_names <- names
      symbol_names <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=ensembl_names,column='SYMBOL', keytype="ENSEMBL", multiVals="first")
    } else {
      symbol_names <- names
      ensembl_names <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=symbol_names,column="ENSEMBL", keytype='SYMBOL', multiVals="first")
    }

    gene_percentage <- dplyr::tibble()
    for (gene_name in top_res$name) {
      print(gene_name)
      if (length(which(gene_name == symbol_names)) > 0) {
        f <- gene_name == symbol_names
        f[is.na(f)] <- F
      } else {
        f <- gene_name == ensembl_names
        f[is.na(f)] <- F
      }

      g_counts <- count_matrix[f, seurat_clusters == c]
      non_g_counts <- count_matrix[f, !(seurat_clusters == c)]

      pct.1 <- sum(g_counts > 0) / length(g_counts)
      pct.2 <- sum(non_g_counts > 0) / length(non_g_counts)

      gene_percentage <- dplyr::bind_rows(gene_percentage, dplyr::tibble(name = gene_name, pct.1=pct.1, pct.2=pct.2))
    }

    top_res <- top_res %>% dplyr::left_join(gene_percentage, by="name")
    scMayoInput <- dplyr::bind_rows(scMayoInput, top_res)
  }

  scMayoInput <- scMayoInput %>%
    dplyr::rename(gene=name, p_val_adj=adj_pval, p_val=pval, pct.1=pct.1, pct.2=pct.2, cluster=cluster, avg_log2FC=lfc)

  if (distinct_marker) {
    input_scMayo$n <- lapply(input_scMayo$gene, function(x) { sum(input_scMayo$gene == x) }) %>% unlist()
    input_scMayo <- input_scMayo %>% dplyr::filter(n == 1) %>% dplyr::select(!n)
  }
  #rownames(scMayoInput) <- scMayoInput$gene
  scMayoInput
}

plot_scMayoOutput <- function(scMayoObj) {
  pred_res <- lapply(1:nrow(scMayoObj$res), function(i) {
    r <- scMayoObj$res[i,]
    s <- strsplit(r$ES.norm, split = ";") %>% unlist() %>% as.numeric()
    c <- strsplit(r$celltype, split = ";") %>% unlist()
    c <- lapply(c, function(x) {stringr::str_replace_all(x, "cell", "")}) %>% unlist()
    c <- lapply(c, function(x) {stringr::str_replace_all(tolower(x), " ", "")}) %>% unlist()
    s
    dplyr::tibble(cell_type_pred = c, score = s, cluster = r$cluster)
  }) %>% do.call('bind_rows', .)

  ground_truth <- computeGroundTruth(seurat_obj)

  ground_truth %>%
    dplyr::left_join(pred_res, by="cluster") %>%
    dplyr::mutate(ground_truth = as.factor(paste0(true_cell_type, ' (', cluster, ')'))) %>%
    ggplot(mapping = aes(x=cell_type_pred, y=ground_truth, fill=score)) +
    geom_tile() +
    scale_fill_continuous(type = "viridis") +
    theme_bw() +
    labs(y = "Annotated cell type", x = "Assigned cell type", fill="Score")
}

computeGroundTruth <- function(seurat_obj) {
  lapply(unique(seurat_obj$seurat_clusters), function(c) {
    table_c <- seurat_obj@meta.data %>%
      dplyr::filter(seurat_clusters == c) %>%
      dplyr::pull(cell_type) %>%
      table()
    table_c <- table_c / sum(table_c)
    max_idx <- which.max(table_c)
    dplyr::tibble(cluster = c, true_cell_type = names(table_c)[max_idx], pct = table_c[max_idx])
  }) %>% do.call('bind_rows', .)
}

