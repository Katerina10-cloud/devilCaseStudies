###-------------------------------------------------------------###
### snATACseq + snRNA data analysis ###
###-------------------------------------------------------------###

### Extract scRNA data ###

# File conversion from AnnData/H5AD to h5Seurat #
#remotes::install_github("mojaveazure/seurat-disk")
# setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")
# atac_muscl.RDS  scrna_muscl.h5ad
#metadata_rna_muscl.Rdata  
#seurat_rna_muscl.RDS

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/sc_devil/")

library(SeuratDisk)
library(tidyverse)

# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
SeuratDisk::Convert("rna_muscl.h5ad", dest = "rna_muscl_seurat.h5seurat")
SeuratObject <- SeuratDisk::LoadH5Seurat("rna_muscl_seurat.h5seurat", misc = FALSE, meta.data = FALSE)

#BiocManager::install("rhdf5")
library(rhdf5)

# to see the structure of the file
rhdf5::h5ls("my_obj.h5ad")
#extract metadata
metadata_rna_muscl <- rhdf5::h5read("rna_muscl.h5ad", "/obs/")
class(metadata)

# add metadata to Seurat object
seurat_rna <- AddMetaData(seurat_rna, metadata = metadata_rna$age_pop, col.name = "age_pop")


### Extract ATAC data ###
peak_counts <- readRDS("peak_counts.RDS")

peak_counts <- atac@assays@data@listData[["PeakMatrix"]]
granges <- atac@rowRanges

genRang_annot <- annot %>% mutate(ranges = paste(annot$seqname,annot$start,annot$end, sep = ":"))
rownames(peak_counts) <- genomic_ranges$ranges

annot <- annot %>% mutate(ranges = paste(annot$seqname,annot$start,annot$end, sep = ":"))
peak_counts <- peak_counts[ rownames(peak_counts) %in% annot$ranges,]
peak_counts <- peak_counts[ ,colnames(peak_counts) %in% rownames(metadata_filtered) ]

metadata_atac <- metadata_atac[order(match(rownames(metadata_atac), colnames(peak_counts))),]

save(genomic_ranges, file = "genomic_ranges.Rdata")
saveRDS(peak_counts, file = "peak_counts.RDS")
saveRDS(granges, file = "granges.RDS")

### Create Seurat object ###

atac_seurat <- CreateSeuratObject(
  counts = peak_counts,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata_atac)

# Create Seurat object with selected cell types #
metadata <- metadata_rna[ (metadata_rna$tech_id %in% c("1") & metadata_rna$cell_cluster %in% c("13", "14")) , ]
metadata <- metadata_rna[ (metadata_rna$tech_id %in% c("1")), ]

metadata <- metadata %>% 
  mutate(age_cluster = case_when(
    age_group == "0"  ~ '1',
    age_group == "1" ~ '0'
    ))  
metadata$age_cluster <- as.factor(metadata$age_cluster)

metadata_atac <- metadata_atac[ (metadata_atac$group %in% c("young", "old") & metadata_atac$cell_type %in% c("Type I", "Type II")) , ]

metadata_atac <- metadata_atac %>% 
  mutate(age_cluster = case_when(
    group == "old"  ~ '1',
    group == "young" ~ '0'
  ))  

peak_counts <- as.matrix(peak_counts)
peak_counts <- peak_counts[ ,colnames(peak_counts) %in% rownames(metadata_atac) ]

metadata_atac <- metadata_atac %>% 
  mutate(cell_type = case_when(
    cell_type == "Adypocyte"  ~ 'Adypocyte',
    cell_type == "EC" ~ 'EC',
    cell_type == "Erytrocyte" ~ 'Erytrocyte',
    cell_type == "FAP" ~ 'FAP',
    cell_type == "Lymphocyte" ~ 'Lymphocyte',
    cell_type == "Mast cell" ~ 'Mast cell',
    cell_type == "MuSC" ~ 'MuSC',
    cell_type == "Myeloid cell" ~ 'Myeloid cell',
    cell_type == "Pericyte" ~ 'Pericyte',
    cell_type == "SMC" ~ 'SMC',
    cell_type == "Schwann cell" ~ 'Schwann cell',
    cell_type == "Specialized MF" ~ 'Specialized MF',
    cell_type == "Fibroblast" ~ 'Fibroblast',
    cell_type == "Type I" ~ 'Myonuclei TI',
    cell_type == "Type II" ~ 'Myonuclei TII'
  ))  


### Peaks Anotation ###

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

grange_object <- readRDS("granges.rds")

annotation <- ChIPseeker::annotatePeak(
  peak = grange_object,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  level = "gene",
  assignGenomicAnnotation = TRUE,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                "Downstream", "Intergenic"),
  annoDb = "org.Hs.eg.db",
  addFlankGeneInfo = FALSE,
  flankDistance = 5000,
  sameStrand = FALSE,
  ignoreOverlap = FALSE,
  ignoreUpstream = FALSE,
  ignoreDownstream = FALSE,
  overlap = "TSS",
  verbose = TRUE,
  columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
)

plotAnnoBar(annotation)

annot <- as.data.frame(annotation)
annot <- annot[ (annot$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")), ]
annot <- annot %>% mutate(ranges = paste(annot$seqnames, annot$start, annot$end, sep = ":"))


# Select DE genes #
sum(res_atac$adj_pval < 0.05 & res_atac$lfc > 0.2, na.rm=TRUE) #up_regulated
sum(res_atac$adj_pval < 0.05 & res_atac$lfc < -0.2, na.rm=TRUE) #down-reg

atac_up <- subset(res_atac, res_atac$adj_pval < 0.05 & res_atac$lfc > 0.2)
atac_down <- subset(res_atac, res_atac$adj_pval < 0.05 & res_atac$lfc < -0.2)
atac_deg <- rbind(atac_up, atac_down)

rna_up <- subset(res, res$adj_pval < 0.05 & res$lfc > 0.2)
rna_down <- subset(res, res$adj_pval < 0.05 & res$lfc < -0.2)
rna_deg <- rbind(rna_up, rna_down)

granAnn <- annot %>% select(SYMBOL, ranges)
colnames(granAnn)[2] <- "name"
res_atac <- merge(res_atac, granAnn, by = "name")
colnames(res_atac)[5] <- "geneID"
colnames(rna_deg)[1] <- "geneID"

atac_up <- subset(res_atac_dup, res_atac_dup$lfc > 0)
atac_up <- atac_up %>% group_by(geneID) %>% arrange(lfc)
atac_up <- atac_up[order(atac_up$lfc, decreasing = TRUE), ] 
atac_up_nodup <- atac_up[!duplicated(atac_up$geneID), ]
atac_up_nodup <- atac_up_nodup[!atac_up_nodup$geneID %in% c("PPARA", "PER2"),]

atac_down <- subset(res_atac_dup, res_atac_dup$lfc < 0)
atac_down <- atac_down %>% group_by(geneID) %>% arrange(lfc)
atac_down <- atac_down[order(atac_down$lfc), ] 
atac_down_nodup <- atac_down[!duplicated(atac_down$geneID), ]

res_atac_nodup <- rbind(atac_up_nodup, atac_down_nodup)
res_atac_nodup <- res_atac_nodup[!duplicated(res_atac_nodup$geneID), ]

deg_atac_nodup <- deg_atac_nodup[ deg_atac_no_dup$geneID %in% rna_deg$geneID,]
rna_deg <- rna_deg[ rna_deg$geneID %in% atac_deg_nodup$geneID,]

atac_up <- subset(res_atac_nodup, res_atac_nodup$adj_pval < 0.05 & res_atac_nodup$lfc >= 1)
atac_down <- subset(res_atac_nodup, res_atac_nodup$adj_pval < 0.05 & res_atac_nodup$lfc <= -1.0)
atac_deg <- rbind(atac_up, atac_down)

rna_up <- subset(rna_deg, rna_deg$adj_pval < 0.05 & rna_deg$lfc >= 1)
rna_down <- subset(rna_deg, rna_deg$adj_pval < 0.05 & rna_deg$lfc <= -1.0)
rna_deg <- rbind(rna_up, rna_down)

atac_deg <- atac_deg[ atac_deg$geneID %in% rna_deg$geneID,]
rna_deg <- rna_deg[ rna_deg$geneID %in% atac_deg$geneID,]

atac_deg <- atac_deg %>% select(geneID, adj_pval, lfc)
colnames(atac_deg) <- c("geneID", "adj_pval_snATAC", "lfc_snATAC")

rna_deg <- rna_deg %>% select(geneID, adj_pval, lfc)
colnames(rna_deg) <- c("geneID", "adj_pval_snRNA", "lfc_snRNA")

res_atac_rna <- merge(atac_deg, rna_deg, by = "geneID")



###-----------------------------------------------------------------###
### Gene set Enrichement analysis ###
###-----------------------------------------------------------------###
rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

# Convert Gene symbols to EntrezID #
hs <- org.Hs.eg.db
my_symbols <- res_atac_rna$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
gene_list <- na.omit(gene_list)

# GSEA using fsea #

# Preparing input #
gene_list <- gene_list[!duplicated(gene_list$SYMBOL),]
gene_list <- gene_list[gene_list$SYMBOL %in% res_atac_rna$geneID,]
gene_list_rank <- as.vector(res_atac_rna$lfc_snRNA)
names(gene_list_rank) <- gene_list$ENTREZID
gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

### GO Enrichment clusterProfiler###

res_gseGO <- gseGO(geneList = gene_list_rank, ont = "BP", OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)

# Select enriched pathways #
res_gse <- res_gseGO@result
res_gse <- res_gse %>%
  filter(Description %in% c("regulation of supramolecular fiber organization", "regulation of actin cytoskeleton organization",
                            "actin filament-based movement", "muscle contraction", "muscle cell proliferation",
                            "skeletal system development")) %>% 
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))

### Visualize fgsea enrichment results ###

res_gse$log_padjust <- -log10(res_gse$p.adjust)

plot1 <- ggdotchart(res_gse, x = "Description", y = "log_padjust",
                    color = "gene_clusters",
                    palette = c("blue", "#FC4E07"),
                    sorting = "descending",
                    rotate = TRUE,
                    group = "gene_clusters",
                    dot.size = 6,
                    add = "segments",
                    title = "Enrichment of key identified genes",
                    xlab = "Pathways",
                    ylab = "-log10(padjust)",
                    ggtheme = theme_pubr()
)
plot1

