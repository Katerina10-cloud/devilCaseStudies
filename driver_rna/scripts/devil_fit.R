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

setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")
sc_retina <- readRDS("seurat_scRetina.rds")
load("results/human_retina/res_retina.Rdata")

# Create scCATCH object from Seurat object
data_retina <- GetAssayData(object = sc_retina, layer = "data")
scObj <- scCATCH::createscCATCH(data = as.matrix(data_retina), as.character(sc_retina@meta.data$seurat_clusters))


## find marker gene for each cluster
#markers <- scCATCH::cellmatch %>%
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

###-----------------------------------------------------------------------###
### Annotation using EasyCellType ###
###-----------------------------------------------------------------------###

BiocManager::install("EasyCellType")
library(EasyCellType)

#Get the expressed markers for each cluster, then convert the gene symbols to Entrez IDs.
library(org.Hs.eg.db)
library(AnnotationDbi)

res_retina_all$entrezid <- mapIds(org.Hs.eg.db,
                           keys=res_retina_all$name, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
res_retina_all <- na.omit(res_retina_all)

markers_sort <- data.frame(gene=res_retina_all$entrezid, cluster=res_retina_all$cluster, 
                           score=res_retina_all$lfc) %>%
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 

input.d <- markers_sort %>% as.data.frame() %>% dplyr::select(gene,cluster,score)


annot.GSEA <- easyct(input.d, 
                     db="panglao", 
                     species="Human", 
                     tissue=c("Eye", "Brain"), 
                     p_cut=0.5,
                     test="GSEA")

data("cellmarker_tissue")
data("clustermole_tissue")
data("panglao_tissue")

plot_dot(test="GSEA", annot.GSEA)

###-------------------------------------------------------------------------------------------------------
res18 <- res %>%
  dplyr::filter(lfc >= 1.0, adj_pval <= 1e-50) %>%
  dplyr::arrange(-lfc)

res18$cluster <- "18"

#res <- res[-c(150,151),]
#res <- res %>%  slice(1:150)

save(res_retina_all, file = "res_retina_all.Rdata")

res18$cluster <- "17"
res_retina_list <- list(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13,res14,res15,res16,res17,res18)
res_retina_all <- do.call("rbind", res_retina_list)

### Calculate percentage of cells expressed particular gene ###
library(scCustomize)
percent_stats <- Percent_Expressing(
  seurat_object,
  features,
  threshold = 0,
  group_by = NULL,
  split_by = NULL,
  entire_object = FALSE,
  slot = deprecated(),
  layer = "data",
  assay = NULL
)

