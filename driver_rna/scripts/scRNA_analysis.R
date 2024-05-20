###---------------------------------------------------------###
### scRNA analysis ###
###---------------------------------------------------------###

#setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")
setwd("~/Documents/PhD_AI/sc_devil/data/muscle/")
#install.packages("Seurat")

#/u/cdslab/kdavydzenka/fast/sc_multiome/blood
#/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/results/blood
#/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/results/human_retina/top_res
# seurat_blood.rds
#/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/devilCaseStudies/cell_types_analysis/utils.R


library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(swissknife)

#PATH <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/"
#sc_retina <- readRDS(file = paste(PATH, "/sc_retina.rds", sep = ""))
sc_retina <- readRDS("/u/cdslab/kdavydzenka/fast/sc_multiome/sc_retina/scRNA_retina.rds")
seurat_blood <- readRDS("/u/cdslab/kdavydzenka/fast/sc_multiome/blood/blood_seurat.rds")
seurat_muscl <- readRDS("seurat_rna_muscl.RDS")

### Rename features in Seurat object ###
meta_features <- seurat_blood@assays[["RNA"]]@meta.features
rownames(seurat_blood@assays$RNA@counts) <- meta_features$feature_name
rownames(seurat_blood@assays$RNA@data) <- meta_features$feature_name

### Extract data from Seurat object ###
counts <- GetAssayData(object = seurat_blood, layer = "counts")
meta_features <- sc_retina@assays[["RNA"]]@meta.features
metadata <- seurat_blood@meta.data

rownames(seurat_rna@meta.data) <- seurat_rna@meta.data$cellID 

dim_red_rna <- seurat_rna@reductions

### Data filtering (donor ID, cell type) ###

# Human retina #

metadata_filtered <- metadata %>%
  filter(sequencing_platform %in% c("unknown"))

metadata_filtered <- metadata_filtered[ !(metadata_filtered$donor_id %in% c("17-010", "17-011", "Hu11", "SC",
                                                          "donor1-hafler", "donor1-scheetz", "donor2-scheetz", "donor2-hafler", "donor3-hafler",
                                                          "donor3-scheetz", "sanes_Pt2", "H1", "H9", "R-00646_03")), ]
metadata_filtered <- metadata_filtered %>% 
  mutate(cluster = case_when(
    cell_type == "GABAergic amacrine cell"  ~ '0',
    cell_type == "H1 horizontal cell"  ~ '1',  
    cell_type == "H2 horizontal cell"  ~ '2',
    cell_type == "Mueller cell"  ~ '3',
    cell_type == "S cone cell" ~ '4',
    cell_type == "amacrine cell" ~ '5',
    cell_type == "astrocyte" ~ '6',
    cell_type == "glycinergic amacrine cell" ~ '7',
    cell_type == "microglial cell" ~ '8',
    cell_type == "midget ganglion cell of retina" ~ '9',
    cell_type == "parasol ganglion cell of retina" ~ '10',
    cell_type == "retinal bipolar neuron" ~ '11',
    cell_type == "retinal cone cell" ~ '12',
    cell_type == "retinal ganglion cell" ~ '13',
    cell_type == "retinal pigment epithelial cell" ~ '14',
    cell_type == "retinal rod cell" ~ '15',
    cell_type == "retinal bipolar cell" ~ '16',
    cell_type == "starburst amacrine cell" ~ '17',
    cell_type == "rod bipolar cell" ~ '18'
    ))

metadata_filtered$cluster <- as.factor(metadata_filtered$cluster)

### Calculate the proportion of mitochondrial reads and add to the metadata table ###
mt.genes <- rownames(rna_counts)[grep("^MT-",rownames(rna_counts))]
percent.mito <- colSums(rna_counts[mt.genes,])/Matrix::colSums(rna_counts)*100
sc_retina <- AddMetaData(sc_retina, percent.mito, col.name = "percent.mito")

### Calculate ribosomal proportion ###
rp.genes <- rownames(rna_counts)[grep("^RP",rownames(rna_counts))]
percent.ribo <- colSums(rna_counts[rb.genes,])/Matrix::colSums(rna_counts)*100
sc_retina <- AddMetaData(sc_retina, percent.ribo, col.name = "percent.ribo")


###---------------------------------------------------------------------------###
### Data Filtering ###
###---------------------------------------------------------------------------###

rna_counts <- rna_counts[,colnames(rna_counts) %in% rownames(metadata)]

## Cell-level filtering ##

#filter the cell with low (low quality libraries) and high (putative doublets) gene detection
total_counts <- colSums(rna_counts)
total_features <- colSums(rna_counts > 0)

mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
feat100_filter <- total_features < 100
feat_mad_filter <- total_features > 5 * mad(total_features)

ribosomal_genes <- grepl("^RPS", rownames(rna_counts)) | grepl("^RPL", rownames(rna_counts))
ribosomal_prop <- colSums(rna_counts[ribosomal_genes, ]) / colSums(rna_counts)
rib_prop_filter <- ribosomal_prop > 0.1

cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter | rib_prop_filter

rna_counts <- rna_counts[, !cell_outliers_filter]
metadata <- metadata[!cell_outliers_filter, ]


## Gene-level filtering ##
#Remove absent genes from dataset:
not_expressed_genes <- which(rowSums(rna_counts) < 10)
rna_counts <- rna_counts[-not_expressed_genes,]

#Remove other genes
rna_counts <- rna_counts[!rownames(rna_counts) %in% mt.genes, ]
rna_counts <- rna_counts[!rownames(rna_counts) %in% rp.genes, ]

# Remove genes that do not present cell-to-cell fluctuations above what is expected due to technical variation
# use the mean-variance trend fit and keep only genes falling above the fitted line

#Using Single cell experiment function
# create SCE
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = rna_counts))

# calculate sizeFactors
libsizes <- colSums(rna_counts)
sizeFactors(sce) <- libsizes / mean(libsizes)

# select variable genes
varGenes <- swissknife::selVarGenes(sce, assay.type="counts")
var_genes <- varGenes[["geneInfo"]]
rna_counts <- rna_counts[rownames(rna_counts) %in% rownames(var_genes), ]


### Create Seurat object ###

seurat_retina <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
seurat_retina <- Seurat::NormalizeData(seurat_retina)
seurat_retina <- Seurat::FindVariableFeatures(seurat_retina)
seurat_retina <- Seurat::ScaleData(seurat_retina)

# Perform clustering
seurat_retina <- Seurat::RunPCA(
  seurat_retina,
  features = Seurat::VariableFeatures(object = seurat_retina))

seurat_retina <- Seurat::FindNeighbors(seurat_retina, dims = 1:20)
seurat_retina <- Seurat::FindClusters(seurat_retina)

# Optional: Run UMAP

seurat_retina <- Seurat::RunUMAP(seurat_retina, dims=1:20)

#n.neighbors = 30, min.dist = 0.3



###-----------------------------------------------------------------###
### Gene set Enrichement analysis ###
###-----------------------------------------------------------------###

# Filter DEG #
resSig_up <- subset(res, res$adj_pval < 0.05 & lfc > 0.5)
resSig_down <- subset(res, res$adj_pval < 0.05 & lfc < -0.5)
deg <- rbind(resSig_down, resSig_up)

# Convert Gene symbols to EntrezID #
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my_symbols <- deg$name
gene_list <- AnnotationDbi::select(hs,
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

### Running fgsea ###
pathways <- reactomePathways(names(gene_list_rank))
fseaRes <- fgsea::fgsea(pathways = pathways,
                         stats = gene_list_rank,
                         minSize = 15,
                         maxSize = 350)

# Make table plot for a bunch of selected pathways #
topPathwaysUp <- fseaRes[ES > 0][head(order(pval), n=7), pathway]
topPathwaysDown <- fseaRes[ES < 0][head(order(pval), n=7), pathway]
topPathways <- c(topPathwaysUp, topPathwaysDown)
plot_fgsea <- plotGseaTable(pathways[topPathways], gene_list_rank, fseaRes, gseaParam = 0.5)
plot_fgsea

# Remove redundant terms #
collapsePathways <- collapsePathways(fseaRes[order(pval)][padj < 0.01],
                                     pathways, gene_list_rank)
mainPathways <- fseaRes[pathway %in% collapsePathways$mainPathways][order(-NES), pathway]
plot_fgsea <- plotGseaTable(pathways[mainPathways], gene_list_rank, fseaRes, gseaParam = 0.5)
plot_fgsea


### GO Enrichment clusterProfiler###
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

res_gseGO <- gseGO(geneList = gene_list_rank, ont = "BP", OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)

# Select enriched pathways #
res_gse <- res_gseGO@result
res_gse <- res_gse %>%
  filter(Description %in% c("neurotransmitter transport", "acetylcholine receptor signaling pathway", "photoreceptor cell development",
                            "synapse assembly", "regulation of protein secretion", "DNA metabolic process", "mitotic cell cycle", 
                            "chromosome organization", "Wnt signaling pathway", "DNA replication" )) %>% 
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))


  

