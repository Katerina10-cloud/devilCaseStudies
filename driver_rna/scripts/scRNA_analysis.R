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

### Data filtering for donor ID and cell type ###

#Differentiated vs non diff. neurons
metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_1", "Donor_2", "Donor_3", "Donor_4", "Donor_5", "Donor_6", "Donor_7", "Donor_8", 
                         "Donor_9", "Donor_10", "Donor_11", "Donor_12"),
         cell_type %in% c("retinal progenitor cell", "Mueller cell", "amacrine cell", "retinal rod cell", "diffuse bipolar 1 cell", 
                          "GABAergic amacrine cell", "retinal bipolar neuron", "retinal cone cell"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata2 <- subset(metadata_filtered,
                    metadata_filtered$nFeatures_RNA >= 1000 & metadata_filtered$nFeatures_RNA <= 8000)

metadata2 <- metadata2 %>% 
  mutate(cell_clusters = case_when(
    cell_type == "retinal progenitor cell"  ~ '0',
    cell_type == "Mueller cell"  ~ '0',  
    cell_type == "amacrine cell"  ~ '1',
    cell_type == "retinal rod cell"  ~ '1',
    cell_type == "diffuse bipolar 1 cell" ~ '1',
    cell_type == "GABAergic amacrine cell" ~ '1',
    cell_type == "retinal bipolar neuron" ~ '1',
    cell_type == "retinal cone cell" ~ '1')) %>% 
  dplyr::select(c("cell_clusters", "cell_type"))


# Tissue specific filtering
metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_1", "Donor_2", "Donor_3", "Donor_4", "Donor_5", "Donor_6", "Donor_7", "Donor_8"),
         tissue %in% c("macula lutea", "peripheral region of retina"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata <- metadata_filtered %>% select("tissue") %>%
  mutate(tissue_clusters = case_when(
    tissue == "macula lutea"  ~ '1',
    tissue == "peripheral region of retina"  ~ '0')) %>%
  select("tissue_clusters")


# Cell specific filtering
metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_1", "Donor_2", "Donor_3", "Donor_4", "Donor_5", "Donor_6", "Donor_7", "Donor_8", 
                         "Donor_9", "Donor_10", "Donor_11", "Donor_12"),
         cell_type %in% c("retinal progenitor cell", "Mueller cell", "retinal ganglion cell", "midget ganglion cell of retina"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata2 <- metadata2 %>% 
  mutate(cell_clusters = case_when(
    cell_type == "retinal progenitor cell"  ~ '0',
    cell_type == "Mueller cell"  ~ '0',  
    cell_type == "retinal ganglion cell"  ~ '1',
    cell_type == "midget ganglion cell of retina" ~ '1')) %>% 
  dplyr::select(c("cell_clusters", "cell_type"))


# Filtering based on development stage
metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_1", "Donor_2", "Donor_3", "Donor_4", "Donor_5", "Donor_6", "Donor_7", "Donor_8", 
                         "Donor_9", "Donor_10", "Donor_11", "Donor_12"),
         cell_type %in% c("retinal progenitor cell", "Mueller cell"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata2 <- metadata_filtered %>% 
  mutate(dev_stage = case_when(
    development_stage == "9th week post-fertilization human stage" & cell_type == "retinal progenitor cell"  ~ '0',
    development_stage == "11th week post-fertilization human stage" & cell_type == "retinal progenitor cell"~ '0',
    development_stage == "12th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '0',
    development_stage == "13th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '0',
    development_stage == "14th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '0',
    development_stage == "15th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '1',
    development_stage == "17th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '1',
    development_stage == "20th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '1',
    development_stage == "21st week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '1',
    development_stage == "24th week post-fertilization human stage" & cell_type == "retinal progenitor cell" ~ '1',
    cell_type == "Mueller cell" ~ '1')) %>% 
  dplyr::select(c("dev_stage", "cell_type"))



cell_barcodes <- rownames(metadata_filtered)
#gene_ids <- meta_features$feature_name
rownames(rna.mtx) <- meta_features$feature_name


#subsetting sparse matrix
rna_counts <- rna.mtx[i=1:36503, j=cell_barcodes, drop = FALSE]

#Convert Large dgCMatrix to normal matrix
rna_counts <- as.matrix(rna.mtx)

### Non-Linear dimensionality reduction ###
seurat <- RunUMAP(seurat, dims = 1:10, n.neighbors = 30, min.dist = 0.3)



### Gene set Enrichement analysis ###

# Filter DEG #
resSig_up <- subset(res, res$adj_pval < 0.05 & lfc > 0.5)
resSig_down <- subset(res, res$adj_pval < 0.05 & lfc < -0.5)
deg <- rbind(resSig_down, resSig_up)

# Convert Gene symbols to EntrezID #
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my_symbols <- deg$name
gene_list <- select(hs,
       keys = my_symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
gene_list <- na.omit(gene_list)

# GSEA using fsea #
#BiocManager::install("fgsea")
#BiocManager::install("reactome.db")
library(reactome.db)
library(fgsea)
library(data.table)
library(ggplot2)
library(tidyverse)

# Preparing input #
gene_list <- gene_list[!duplicated(gene_list$SYMBOL),]
gene_list <- gene_list[gene_list$SYMBOL %in% deg$name,]
gene_list_rank <- as.vector(deg$lfc)
names(gene_list_rank) <- gene_list$ENTREZID
gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

# Running fgsea #
pathways <- reactomePathways(names(gene_list_rank))
fseaRes <- fgsea::fgsea(pathways = pathways,
                         stats = gene_list_rank,
                         minSize = 15,
                         maxSize = 500)

# Make table plot for a bunch of selected pathways #
topPathwaysUp <- fseaRes[ES > 0][head(order(pval), n=5), pathway]
topPathwaysDown <- fseaRes[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plot_fgsea <- plotGseaTable(pathways[topPathways], gene_list_rank, fseaRes, gseaParam = 0.5)
plot_fgsea

# Remove redundant terms #
collapsePathways <- collapsePathways(fseaRes[order(pval)][padj < 0.01],
                                     pathways, gene_list_rank)
mainPathways <- fseaRes[pathway %in% collapsePathways$mainPathways][order(-NES), pathway]
plot_fgsea <- plotGseaTable(pathways[mainPathways], gene_list_rank, fseaRes, gseaParam = 0.5)
plot_fgsea

