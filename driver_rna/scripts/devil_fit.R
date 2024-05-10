###----------------------------------------------------------###
### Devil testing ###
###----------------------------------------------------------###

#devtools::install_github("caravagnalab/devil")

library(devil)
library(tidyverse)

#Create design matrix
design_matrix <- model.matrix(~cell_type, data = metadata)

### Parameters inference ###

devil_fit_blood <- devil::fit_devil(input_matrix = rna_counts,
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

res <- devil::test_de(devil.fit = devil_fit_retina,
                      contrast = c(0,1),
                      pval_adjust_method = "BH",
                      max_lfc = 10,
                      clusters = NULL)


###----------------------------------------------------------------------------###
### Devil test on cluster ###
###----------------------------------------------------------------------------###

c <- seurat_obj$seurat_clusters[1]
whole_res <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
  print(c)
  metadata$cluster <- seurat_obj$seurat_clusters
  metadata$group <- metadata$cluster == c
  
  design_matrix <- model.matrix(~group, metadata)
  
  counts <- as.matrix(counts)
  gg <- ((counts[,metadata$group] %>% rowSums()) == 0)
  bad_genes <- gg[gg == T] %>% names()
  length(bad_genes)
  
  gg <- ((counts[,!metadata$group] %>% rowSums()) == 0)
  bad_genes <- c(bad_genes, gg[gg == T] %>% names())
  length(bad_genes)
  
  fit <- devil::fit_devil(as.matrix(counts[!(rownames(counts) %in% bad_genes),]), design_matrix, size_factors = T, verbose = T, parallel = T)
  #fit <- devil::fit_devil(as.matrix(counts), design_matrix, size_factors = T, verbose = T, parallel = T)
  
  clusters <- as.numeric(as.factor(metadata$donor))
  res <- devil::test_de(fit, contrast = c(0,1), clusters = clusters, max_lfc = Inf)
  
  remove_genes <- grepl("^ENS", res$name)
  res <- res[!remove_genes, ]
  min_pval <- res %>% dplyr::filter(lfc > .5) %>% dplyr::pull(adj_pval) %>% min()
  if (min_pval > 1e-50) {
    top_res <- res %>%
      dplyr::filter(lfc > .2, adj_pval <= min_pval * 100) %>%
      dplyr::arrange(-lfc) %>%
      dplyr::slice(1:40)
  } else {
    top_res <- res %>%
      dplyr::filter(lfc > .2, adj_pval <= 1e-50) %>%
      dplyr::arrange(-lfc) %>%
      dplyr::slice(1:40)
  }
  top_res$cluster <- c
  
  gene_percentage <- lapply(top_res$name, function(gene_name) {
    g_counts <- counts[gene_name, metadata$group]
    non_g_counts <- counts[gene_name, !metadata$group]
    
    pct.1 <- sum(g_counts > 0) / length(g_counts)
    pct.2 <- sum(non_g_counts > 0) / length(non_g_counts)
    
    dplyr::tibble(name = gene_name, pct.1=pct.1, pct.2=pct.2)
  }) %>% do.call("bind_rows", .)
  
  top_res <- top_res %>% dplyr::left_join(gene_percentage, by="name")
  top_res
}) %>% do.call("bind_rows", .)



###-----------------------------------------------------------------------------###
### Cell clusters anotation
###-----------------------------------------------------------------------------###
#install.packages("scCATCH")
library(scCATCH)

setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")
sc_retina <- readRDS("seurat_scRetina.rds")
load("results/blood/top_res/res_blood_top.Rdata")

# Create scCATCH object from Seurat object
data_blood <- GetAssayData(object = seurat_blood, layer = "data")
scObj <- scCATCH::createscCATCH(data = as.matrix(data_blood), as.character(seurat_blood@meta.data$seurat_clusters))


## find marker gene for each cluster
markers <- scCATCH::cellmatch %>%
  dplyr::filter(species == "Human", gene %in% res_retina$name)

#obj <- scCATCH::findmarkergene(scObj, species = "Human", marker = markers, if_use_custom_marker = TRUE, tissue = c("Eye", "Retina", "Retinal pigment epithelium"))

obj <- scCATCH::findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = c("Eye", "Retina", "Retinal pigment epithelium"))
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

celltype_data$cluster <- as.numeric(as.factor(celltype_data$cluster))


###-----------------------------------------------------------------------###
### Annotation using EasyCellType ###
###-----------------------------------------------------------------------###

#BiocManager::install("EasyCellType")
library(EasyCellType)

#Get the expressed markers for each cluster, then convert the gene symbols to Entrez IDs.
library(org.Hs.eg.db)
library(AnnotationDbi)

res_blood_top$entrezid <- mapIds(org.Hs.eg.db,
                           keys=res_blood_top$name, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
res_blood_top <- na.omit(res_blood_top)

markers_sort <- data.frame(gene=res_blood_top$entrezid, cluster=res_blood_top$cluster, 
                           score=res_blood_top$lfc) %>%
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 

input.d <- markers_sort %>% as.data.frame() %>% dplyr::select(gene,cluster,score)


annot.GSEA <- easyct(input.d, 
                     db="cellmarker", 
                     species="Human", 
                     tissue=c("Blood"), 
                     p_cut=0.7,
                     scoretype="pos",
                     test="GSEA")

data("cellmarker_tissue") 
cellmarker_tissue[["Human"]] %>% unique()
data("clustermole_tissue")
data("panglao_tissue")

plot_dot(test="GSEA", annot.GSEA)


###------------------------------------------------------------------###
### Annotation using scMayoMap ###
###------------------------------------------------------------------###

setwd("~/Documents/PhD_AI/sc_devil/results/blood/")
#devtools::install_github("chloelulu/scMayoMap")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
sapply(pkgs, require, character.only = TRUE)
library(scMayoMap)

dd <- scMayoMap::scMayoMapDatabase
dd$tissue %>% unique()

input_scMayo <- whole_res_blood %>% 
  dplyr::rename(gene=name, p_val_adj=adj_pval, p_val=pval, pct.1=pct.1, pct.2=pct.2, cluster=cluster, avg_log2FC=lfc) %>% 
  as.data.frame()
rownames(input_scMayo) <- input_scMayo$gene

#input_scMayo %>% pull(cluster) %>% table()

obj <- scMayoMap(data = input_scMayo, tissue = 'blood')
obj$res
plt <- scMayoMap.plot(scMayoMap.object = obj, width = 12, height = 8)




