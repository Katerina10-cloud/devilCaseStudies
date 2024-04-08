#!/usr/bin/env Rscript

### scRNA analysis ###

setwd("D:/PhD_AI/sc_devil/")

#setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")

library(tidyverse)
library(ggplot2)
library(Seurat)
library(scater)
library(scran)
#library(BPCells) #high performance single cell analysis
#library(SeuratObject)
#library(sva)
#library(BiocParallel)
#library(MatrixExtra)

PATH <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/"
sc_retina <- readRDS(file = paste(PATH, "/sc_retina.rds", sep = ""))
#sc_retina <- readRDS("/u/cdslab/kdavydzenka/fast/sc_multiome/sc_retina/scRNA_retina.rds")

#Extract data from Seurat object
rna_counts <- GetAssayData(object = sc_retina, layer = "counts")
meta_features <- sc_retina@assays[["RNA"]]@meta.features
metadata <- sc_retina@meta.data

#Data filtering for donor ID and cell type
metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_1"),
         cell_type %in% c("retinal progenitor cell", "retinal ganglion cell"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata <- metadata_filtered %>% select("cell_type") %>%
  mutate(cell_clusters = case_when(
    cell_type == "retinal progenitor cell"  ~ '0',
    #cell_type == "retinal rod cell"  ~ '1',
    cell_type == "retinal ganglion cell"  ~ '1')) %>%
  select("cell_clusters")

cell_barcodes <- rownames(metadata_filtered)
#gene_ids <- meta_features$feature_name
rownames(rna.mtx) <- meta_features$feature_name


#subsetting sparse matrix
rna_counts <- rna.mtx[i=1:36503, j=cell_barcodes, drop = FALSE]

#rna_counts <- rna_counts %>%
  #select(colnames(rna_counts) %in% rownames(metadata_filtered))
#atac_counts_ret <- subset(atac_counts_retina, subset = colnames(atac_counts_retina) %in% rownames(metadata_filtered_atac))

#Convert Large dgCMatrix to normal matrix
rna_counts <- as.matrix(rna.mtx)

### Non-Linear dimensionality reduction ###
seurat <- RunUMAP(seurat, dims = 1:10, n.neighbors = 30, min.dist = 0.3)

#saving plot cluster
options(bitmapType='cairo')
png(file="plots/umap.png", width = 480, height = 480)

DimPlot(seurat, 
        dims = c(1, 2),
        reduction = "umap",
        group.by = 'cell_type',
        split.by = 'cell_type')

dev.off()

