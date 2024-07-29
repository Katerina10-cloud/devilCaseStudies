rm(list=ls())
pkgs <- c("patchwork", "ggpubr","ggplot2", "dplyr","tidyr","tibble","reshape2", "stringr", "AnnotationDbi", "scMayoMap", "Seurat", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
source("utils.R")

set.seed(SEED)

args = commandArgs(trailingOnly=TRUE)

## Input data
dataset_name <- "pbmc"
tissue <- "blood"
seurat_obj <- readRDS(paste0('results/', dataset_name, '/seurat.RDS'))


# Upset Plot ####
#mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")
devil_res <- readRDS(paste0('results/', dataset_name, '/devil.RDS'))
nebula_res <- readRDS(paste0('results/', dataset_name, '/nebula.RDS'))
glm_res <- readRDS(paste0('results/', dataset_name, '/glmGamPoi.RDS'))

lfc_cut <- 1
n_markers <- 1000000
pval_cut <- .05
all_markers <- c()
c <- 2
ks_overlap <- dplyr::tibble()
c <- 3
for (c in unique(seurat_obj$seurat_clusters)) {
  for (m in c("devil", "nebula", "glmGamPoi")) {
    de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
    # if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    #   suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
    # }

    top_genes <- de_res %>%
      dplyr::group_by(cluster) %>%
      dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::arrange(-lfc) %>%
      dplyr::filter(cluster == c) %>%
      dplyr::slice(1:n_markers) %>% pull(name) %>% unname()

    all_markers <- c(all_markers, top_genes) %>% unique() %>% na.omit()
  }

  mm <- lapply(c("devil", "nebula", "glmGamPoi"), function(m) {
    de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
    # if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    #   suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
    # }

    top_genes <- de_res %>%
      dplyr::group_by(cluster) %>%
      dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::arrange(-lfc) %>%
      dplyr::filter(cluster == c) %>%
      dplyr::slice(1:n_markers) %>% pull(name) %>% unname()

    c((all_markers %in% top_genes) %>% as.numeric())
  }) %>% do.call('rbind', .)


  mm_tibble <- lapply(1:ncol(mm), function(i) {
    print(i)
    g <- all_markers[i]
    v <- mm[,i] %>% unlist()
    dplyr::tibble(gene = g, devil=v[1], nebula=v[2], glmGamPoi=v[3])
  }) %>% do.call("bind_rows", .) %>% as.data.frame()

  method_order <- mm_tibble[,2:4] %>% colSums() %>% sort(decreasing = T)
  method_colors_arranged <- c()

  for (m in names(method_order)) {
    method_colors_arranged <- c(method_colors_arranged, method_colors[which(names(method_colors) == m)])
  }

  #upset_plot <- UpSetR::upset(data.frame(mm_tibble), sets.bar.color = method_colors_arranged,
  #                            order.by = "freq", text.scale = 1.8)
  #upset_plot

  idx_cluster <- which(seurat_obj$seurat_clusters == c)
  idx_others <- which(!(seurat_obj$seurat_clusters == c))

  if (length(idx_others) > length(idx_cluster)) {
    set.seed(SEED)
    idx_others <- sample(idx_others, as.integer(length(idx_others) * .05), replace = F)
    #idx_others <- sample(idx_others, length(idx_cluster), replace = F)
  }

  idxs <- c(idx_cluster, idx_others)

  counts <- as.matrix(seurat_obj@assays$RNA$counts)
  design_matrix <- model.matrix(~group, dplyr::tibble(group = seurat_obj$seurat_clusters == c))
  dm <- design_matrix[idxs,]
  cc <- counts[,idxs]
  #cc <- cc[rowMeans(cc) > .005,]
  #rownames(cc) <- mapIds(org.Hs.eg.db, keys=rownames(cc),column="SYMBOL", keytype="ENSEMBL", multiVals="first")

  n <- mm_tibble$gene[1]
  ks_res <- lapply(mm_tibble$gene, function(n) {
    gene_idx <- which(rownames(cc) == n)[1]

    if (is.na(gene_idx)) {
      print(n)
      return(NULL)
    }
    if (is.na(n)) {return(NULL)}
    d <- dplyr::tibble(x = cc[gene_idx,], class = as.factor(rowSums(dm)))

    # get cdf of two samples
    empirical_cdf <- function(x) {
      # Sort the observations
      x_sorted <- sort(x)
      # Get the number of observations
      n <- length(x_sorted)

      # Define the empirical CDF function
      ecdf_function <- function(value) {
        # Compute the proportion of observations less than or equal to the given value
        sum(x_sorted <= value) / n
      }

      # Return the empirical CDF function
      return(ecdf_function)
    }

    ecdf_1 <- empirical_cdf(d %>% dplyr::filter(class == 1) %>% pull(x))
    ecdf_2 <- empirical_cdf(d %>% dplyr::filter(class == 2) %>% pull(x))

    max_c <- max(d$x)

    d_ks <- dplyr::bind_rows(
      dplyr::tibble(x = 0:max_c, class='1', y = unlist(lapply(0:max_c, ecdf_1))),
      dplyr::tibble(x = 0:max_c, class='2', y = unlist(lapply(0:max_c, ecdf_2)))
    )

    # d_ks %>%
    #   ggplot(mapping = aes(x=x, y=y, col=class)) +
    #   geom_point() +
    #   geom_line() +
    #   ggtitle(paste0("Gene ", nebula_genes$gene[gene_idx]))

    ks <- ks.test(
      x = d_ks %>% dplyr::filter(class == '1') %>% dplyr::pull(y),
      y = d_ks %>% dplyr::filter(class == '2') %>% dplyr::pull(y)
    )

    dplyr::tibble(name = rownames(cc)[gene_idx], pvalue = ks$p.value)
  }) %>% do.call('bind_rows', .)

  mm_tibble$KS <- as.numeric(ks_res$pvalue <= .05)
  colnames(mm_tibble)
  rrr <- lapply(c('devil', 'nebula', 'glmGamPoi'), function(m) {
    idx <- mm_tibble[,m] == 1
    perc <- ((mm_tibble[idx,m] == mm_tibble$KS[idx]) %>% sum()) / length(idx)
    dplyr::tibble(m = m, perc=perc * 100)
  }) %>% do.call('bind_rows', .) %>% dplyr::mutate(cluster = c)
  ks_overlap <- dplyr::bind_rows(ks_overlap, rrr)
}

ks_overlap %>%
  ggplot(mapping = aes(x=as.numeric(cluster), y=perc, col=m)) +
  geom_point() +
  geom_line()
ggsave("plot_figure/ks_overlap.pdf", dpi = 300, width = 10, height = 8, plot = last_plot())


