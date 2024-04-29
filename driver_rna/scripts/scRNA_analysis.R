###---------------------------------------------------------###
### scRNA analysis ###
###---------------------------------------------------------###

#setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")

library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(swissknife)

#PATH <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/"
#sc_retina <- readRDS(file = paste(PATH, "/sc_retina.rds", sep = ""))
sc_retina <- readRDS("/u/cdslab/kdavydzenka/fast/sc_multiome/sc_retina/scRNA_retina.rds")

### Rename features in Seurat object ###
meta_features <- sc_retina@assays[["RNA"]]@meta.features
rownames(sc_retina@assays$RNA@counts) <- meta_features$feature_name
rownames(sc_retina@assays$RNA@data) <- meta_features$feature_name

### Extract data from Seurat object ###
rna_counts <- GetAssayData(object = sc_retina, layer = "counts")
meta_features <- sc_retina@assays[["RNA"]]@meta.features
metadata <- sc_retina@meta.data

### Data filtering (donor ID, cell type) ###

metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_2", "Donor_3", "Donor_4", "Donor_5", "Donor_6", "Donor_8", "Donor_10", "Donor_11", "Donor_12"),
         #cell_type %in% c("retinal progenitor cell", "Mueller cell", "amacrine cell", "retinal rod cell", 
                          #"diffuse bipolar 1 cell", "GABAergic amacrine cell", "retinal bipolar neuron", "retinal cone cell"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))


metadata_filtered <- metadata_filtered[ !(metadata_filtered$cell_type %in% c("diffuse bipolar 4 cell", "diffuse bipolar 6 cell",
                                                                             "S cone cell", "OFFx cell", "OFF parasol ganglion cell")), ]

metadata_filtered <- metadata_filtered %>% 
  mutate(cluster = case_when(
    cell_type == "amacrine cell"  ~ '1',
    cell_type == "retinal cone cell"  ~ '2',  
    cell_type == "retinal rod cell"  ~ '3',
    cell_type == "Mueller cell"  ~ '4',
    cell_type == "retinal ganglion cell" ~ '5',
    cell_type == "retinal horizontal cell" ~ '6',
    cell_type == "retinal bipolar neuron" ~ '7',
    cell_type == "ON-bipolar cell" ~ '8',
    cell_type == "rod bipolar cell" ~ '9',
    cell_type == "retinal progenitor cell" ~ '10',
    cell_type == "H1 horizontal cell" ~ '11',
    cell_type == "H2 horizontal cell" ~ '12',
    cell_type == "starburst amacrine cell" ~ '13',
    cell_type == "midget ganglion cell of retina" ~ '14',
    cell_type == "GABAergic amacrine cell" ~ '15',
    cell_type == "glycinergic amacrine cell" ~ '16',
    cell_type == "diffuse bipolar 1 cell" ~ '17',
    cell_type == "diffuse bipolar 2 cell" ~ '18',
    cell_type == "diffuse bipolar 3a cell" ~ '19',
    cell_type == "diffuse bipolar 3b cell" ~ '20',
    cell_type == "flat midget bipolarvcell" ~ '21',
    cell_type == "invaginating midget bipolar cell" ~ '22',
    cell_type == "ON parasol ganglion cell" ~ '23'))


metadata2 <- metadata2 %>% 
  mutate(cell_clusters = case_when(
    cell_type == "retinal progenitor cell"  ~ '0',
    cell_type == "Mueller cell"  ~ '0',  
    cell_type == "retinal ganglion cell"  ~ '1',
    cell_type == "midget ganglion cell of retina" ~ '1')) %>% 
  dplyr::select(c("cell_clusters", "cell_type"))


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

rna_counts <- rna_counts[,colnames(rna_counts) %in% rownames(metadata_filtered)]

## Cell-level filtering ##

#filter the cell with low (low quality libraries) and high (putative doublets) gene detection
total_counts <- colSums(rna_counts)
total_features <- colSums(rna_counts > 0)

mad5_filter <- total_counts > median(total_counts) + 5 * mad(total_counts)
feat100_filter <- total_features < 100
feat_mad_filter <- total_features > 5 * mad(total_features)

ribosomal_genes <- grepl("^RPS", rownames(rna_counts)) | grepl("^RPL", rownames(rna_counts))
ribosomal_prop <- colSums(rna_counts[ribosomal_genes, ]) / colSums(rna_counts)
rib_prop_filter <- ribosomal_prop > .1

cell_outliers_filter <- mad5_filter | feat100_filter | feat_mad_filter | rib_prop_filter

rna_counts <- rna_counts[, !cell_outliers_filter]
metadata_filtered <- metadata_filtered[!cell_outliers_filter, ]


## Gene-level filtering ##
#Remove absent genes from dataset:
not_expressed_genes <- which(rowSums(rna_counts)<=0.05)
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


### Create Seurat object

seurat_scRetina = SeuratObject::CreateSeuratObject(counts = rna_counts, meta.data = metadata)
seurat_scRetina <- Seurat::NormalizeData(seurat_scRetina)
seurat_scRetina <- Seurat::FindVariableFeatures(seurat_scRetina)
seurat_scRetina <- Seurat::ScaleData(seurat_scRetina)

# Perform clustering
seurat_scRetina <- Seurat::RunPCA(
  seurat_scRetina,
  features = Seurat::VariableFeatures(object = seurat_scRetina)
)

seurat_scRetina <- Seurat::FindNeighbors(seurat_scRetina)
seurat_scRetina <- Seurat::FindClusters(seurat_scRetina)

# Optional: Run UMAP

seurat_scRetina <- Seurat::RunUMAP(seurat_scRetina, dims=1:20)

#seurat <- RunUMAP(seurat, dims = 1:10, n.neighbors = 30, min.dist = 0.3)




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

res_gse <- res_gse %>%
  filter(Description %in% c("DNA metabolic process", "mitotic cell cycle", "DNA replication",
                            "canonical Wnt signaling pathway", "chromosome organization",
                            "neurotransmitter secretion", "dopamine receptor signaling pathway",
                            "synapse organization", "eye photoreceptor cell development","myotube cell development" )) %>% 
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))


res_gse <- res_gse %>% mutate(gene_clusters = case_when(
  NES > 0  ~ 'up-regulated',
  NES < 0  ~ 'down-regulated'))
  


