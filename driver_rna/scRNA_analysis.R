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
  filter(donor_id %in% c("Donor_1", "Donor_2", "Donor_3", "Donor_4", "Donor_5", "Donor_6", "Donor_7", "Donor_8", 
                         "Donor_9", "Donor_10", "Donor_11", "Donor_12"),
         cell_type %in% c("retinal progenitor cell", "retinal ganglion cell", "midget ganglion cell of retina", "Mueller cell"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata <- metadata_filtered %>% select("cell_type") %>%
  mutate(cell_clusters = case_when(
    cell_type == "retinal progenitor cell"  ~ '0',
    cell_type == "midget ganglion cell of retina"  ~ '1',  
    cell_type == "retinal ganglion cell"  ~ '1',
    cell_type == "Mueller cell"  ~ '1')) %>%
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


