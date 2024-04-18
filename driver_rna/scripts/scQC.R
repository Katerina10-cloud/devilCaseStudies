#!/usr/bin/env Rscript

#### Quality control of data for filtering cells using Seurat and Scater packages####

require(tidyverse)
require(Seurat)
require(scater)
require(swissknife)
require(Matrix)
require(ggplot2)
require(gridExtra)

setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")

sc_retina <- readRDS("/u/cdslab/kdavydzenka/fast/sc_multiome/sc_retina/scRNA_retina.rds")

## Rename features in Seurat object##
meta_features <- sc_retina@assays[["RNA"]]@meta.features
#counts <- GetAssayData(sc_retina,assay = "RNA",layer = "counts")

rownames(sc_retina@assays$RNA@counts) <- meta_features$feature_name
rownames(sc_retina@assays$RNA@data) <- meta_features$feature_name

## Calculate the proportion of mitochondrial reads and add to the metadata table ##
rna_counts <- GetAssayData(object = sc_retina, layer = "counts")
mt.genes <- rownames(sc_retina)[grep("^MT-",rownames(sc_retina))]
percent.mito <- colSums(rna_counts[mt.genes,])/Matrix::colSums(rna_counts)*100
sc_retina <- AddMetaData(sc_retina, percent.mito, col.name = "percent.mito")

## Calculate ribosomal proportion ##
rp.genes <- rownames(sc_retina)[grep("^RP",rownames(sc_retina))]
percent.ribo <- colSums(rna_counts[rb.genes,])/Matrix::colSums(rna_counts)*100
sc_retina <- AddMetaData(sc_retina, percent.ribo, col.name = "percent.ribo")

## Plot QC ##
#mt.genes <- grep(pattern = "MT-", x = rownames(rna.mtx), value = TRUE)
#mt_percentage.per.cell <- 100*colSums(rna.mtx[mt.genes,])/colSums(rna.mtx)
#rp_percentage.per.cell <- 100*colSums(rna.mtx[rp.genes,])/colSums(rna.mtx)

p1 <- hist(nReads.per.cell, breaks = 100) %>% abline(v=1000, col="red", lwd=2)
p2 <- hist(nGenes.per.cell, breaks=100) %>% abline(v=1000, col="red", lwd=2) %>% abline(v=6000, col="red", lwd=2)
p3 <- hist(mt_percentage.per.cell, breaks=100) %>% abline(v=5, col="red", lwd=2)
p4 <- hist(percent.ribo, breaks=100) %>% abline(v=5, col="red", lwd=2)


### Filtering ###

## Cell-level filtering ##
#filter the cell with low (low quality libraries) and high (putative doublets) gene detection
nGenes.per.cell <- Matrix::colSums(rna.mtx>0)

cells.to.keep <- mt_percentage.per.cell<5 & nGenes.per.cell>1000 & nGenes.per.cell<8000 & nReads.per.cell>100

metadata <- subset(metadata_filtered, 
                   metadata_filtered$nFeatures_RNA >= 3000 & metadata_filtered$nFeatures_RNA <= 8000)
rna_counts <- rna_counts[,colnames(rna_counts) %in% rownames(metadata)]

## Gene-level filtering ##
#Remove absent genes from dataset:
absent_genes=which(rowSums(sc_retina)==0)
sc_retina=sc_retina[-absent_genes,]

#Remove other genes
rna_counts <- rna_counts[!rownames(rna_counts) %in% mt.genes, ]
rna_counts <- rna_counts[!rownames(rna_counts) %in% rp.genes, ]

# Remove genes that do not present cell-to-cell fluctuations above what is expected due to technical variation
# use the mean-variance trend fit and keep only genes falling above the fitted line

#Using Single cell experiment function
# create SCE
sce <- SingleCellExperiment(list(counts = rna_counts))

# calculate sizeFactors
libsizes <- colSums(rna_counts)
sizeFactors(sce) <- libsizes / mean(libsizes)

# select variable genes
varGenes <- swissknife::selVarGenes(sce, assay.type="counts")
plot <- swissknife::plotSelVarGenes(varGenes)


