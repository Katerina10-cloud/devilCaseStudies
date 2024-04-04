#!/usr/bin/env Rscript

#### Quality control of data for filtering cells using Seurat and Scater packages####

require(tidyverse)
require(Seurat)
require(scater)
require(Matrix)
library(ggplot2)
library(gridExtra)

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
library(gridExtra)

p1 <- hist(nFeature_RNA, breaks=100) %>% abline(v=1000, col="red", lwd=2)
p2 <- hist(mt_percentage.per.cell, breaks=100) %>% abline(v=5, col="red", lwd=2)
p3 <- hist(rp_percentage.per.cell, breaks=100) %>% abline(v=5, col="red", lwd=2)

plot1 <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave(plot1,
       filename = "plot1.png",
       device = "pdf",
       height = 6, width = 5, units = "in")

### Filtering ###
#Cell-level filtering#
#filter the cell with low (low quality libraries) and high (putative doublets) gene detection
nGenes.per.cell <- Matrix::colSums(rna.mtx>0)
p2 <- hist(nGenes.per.cell, breaks=100) %>% abline(v=1000, col="red", lwd=2)


#Gene-level filtering#
#Remove absent genes from dataset:
absent_genes=which(rowSums(sc_retina)==0)
sc_retina=sc_retina[-absent_genes,]

#Remove other genes
rna_counts <- rna_counts[!rownames(rna_counts) %in% mt.genes, ]
rna_counts <- rna_counts[!rownames(rna_counts) %in% rp.genes, ]
