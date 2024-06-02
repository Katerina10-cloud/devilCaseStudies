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

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis/results/devil/")

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
metadata$patient_id <- as.factor(metadata$patient_id)

metadata_atac <- metadata_atac[ (metadata_atac$group %in% c("young", "old") & metadata_atac$cell_type %in% c("Type I", "Type II")) , ]

metadata_atac <- metadata_atac %>% 
  mutate(age_cluster = case_when(
    group == "old"  ~ '1',
    group == "young" ~ '0'
  ))  

peak_counts <- as.matrix(peak_counts)
peak_counts <- peak_counts[ ,colnames(peak_counts) %in% rownames(metadata_atac) ]

metadata_atac$cell_type <- as.character(metadata_atac$cell_type)
conditions <- c("Type I", "Type II")
replacement_values <- c("Myonuclei TI", "Myonuclei TII")
metadata_atac$cell_type <- replace(metadata_atac$cell_type, metadata_atac$cell_type %in% conditions, replacement_values)
metadata_atac$cell_type <- as.factor(metadata_atac$cell_type)

metadata_rna <- metadata_rna %>% 
  mutate(cell_type = case_when(
    cell_cluster == "0"  ~ 'Adypocyte',
    cell_cluster == "1" ~ 'EC',
    cell_cluster == "2" ~ 'Erytrocyte',
    cell_cluster == "3" ~ 'FAP',
    cell_cluster == "4" ~ 'Lymphocyte',
    cell_cluster == "5" ~ 'Mast cell',
    cell_cluster == "6" ~ 'MuSC',
    cell_cluster == "7" ~ 'Myeloid cell',
    cell_cluster == "8" ~ 'Pericyte',
    cell_cluster == "9" ~ 'SMC',
    cell_cluster == "10" ~ 'Schwann cell',
    cell_cluster == "11" ~ 'Specialized MF',
    cell_cluster == "12" ~ 'Fibroblast',
    cell_cluster == "13" ~ 'Myonuclei TI',
    cell_cluster == "14" ~ 'Myonuclei TII'
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


### Select DE genes ###

#res_rna_neb[is.na(res_rna_neb)] <- 1

sum(res_neb_atac$adj_pval < 0.05 & res_neb_atac$lfc > 0.2, na.rm=TRUE) #up_regulated
sum(res_neb_atac$adj_pval < 0.05 & res_neb_atac$lfc < -0.2, na.rm=TRUE) #down-reg

atac_up <- subset(res_atac, res_atac$adj_pval < 0.05 & res_atac$lfc > 0.2)
atac_down <- subset(res_atac, res_atac$adj_pval < 0.05 & res_atac$lfc < -0.2)
atac_deg <- rbind(atac_up, atac_down)

rna_up <- subset(res_rna, res_rna$adj_pval < 0.05 & res_rna$lfc > 0.2)
rna_down <- subset(res_rna, res_rna$adj_pval < 0.05 & res_rna$lfc < -0.2)
rna_deg <- rbind(rna_up, rna_down)

granAnn <- annot %>% select(SYMBOL, ranges)
colnames(granAnn)[2] <- "name"
res_atac <- merge(atac_deg, granAnn, by = "name")
colnames(res_atac)[5] <- "geneID"
colnames(rna_deg)[1] <- "geneID"

res_atac_dup <- res_atac
atac_up <- subset(res_atac_dup, res_atac_dup$lfc > 0)
atac_up <- atac_up %>% group_by(geneID) %>% arrange(lfc)
atac_up <- atac_up[order(atac_up$lfc, decreasing = TRUE), ] 
atac_up_nodup <- atac_up[!duplicated(atac_up$geneID), ]
#atac_up_nodup <- atac_up_nodup[!atac_up_nodup$geneID %in% c("PPARA", "PER2"),]

atac_down <- subset(res_atac_dup, res_atac_dup$lfc < 0)
atac_down <- atac_down %>% group_by(geneID) %>% arrange(lfc)
atac_down <- atac_down[order(atac_down$lfc), ] 
atac_down_nodup <- atac_down[!duplicated(atac_down$geneID), ]

res_atac_nodup <- rbind(atac_up_nodup, atac_down_nodup)
res_atac_nodup <- res_atac_nodup[!duplicated(res_atac_nodup$geneID), ]

atac_deg_nodup <- atac_deg_nodup[ atac_deg_nodup$geneID %in% rna_deg$geneID,]
rna_deg <- rna_deg[ rna_deg$geneID %in% atac_deg_nodup$geneID,]

atac_up <- subset(res_atac_nodup, res_atac_nodup$adj_pval < 0.05 & res_atac_nodup$lfc >= 1)
atac_down <- subset(res_atac_nodup, res_atac_nodup$adj_pval < 0.05 & res_atac_nodup$lfc <= -1.0)
atac_deg_nodup <- rbind(atac_up, atac_down)
atac_deg <- atac_deg_nodup

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


save(res_atac_rna, file = "overlap_atac_rna_nebula.Rdata")



###-----------------------------------------------------------------###
### Gene set Enrichement analysis ###
###-----------------------------------------------------------------###
rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

# Convert Gene symbols to EntrezID #
hs <- org.Hs.eg.db
my_symbols <- overlap_devil$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
gene_list <- na.omit(gene_list)


# GSEA using fsea #

# Preparing input #
gene_list <- gene_list[!duplicated(gene_list$SYMBOL),]
gene_list <- gene_list[gene_list$SYMBOL %in% overlap_devil$geneID,]
gene_list_rank <- as.vector(overlap_devil$lfc_snRNA)
names(gene_list_rank) <- gene_list$ENTREZID
gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

-----------------------------------------------------------
### GO Enrichment clusterProfiler###
-----------------------------------------------------------  

res_gseGO <- clusterProfiler::gseGO(geneList = gene_list_rank, ont = "BP", OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)

# Select enriched pathways #
res_gse <- res_gseGO@result
res_gse <- res_gse %>%
  filter(Description %in% c("cellular response to salt", "negative regulation of metabolic process", "apoptotic process",
                            "actin filament-based process", "muscle contraction", "muscle system process", "regulation of supramolecular fiber organization", 
                            "muscle cell proliferation")) %>% 
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))

### Visualize fgsea enrichment results ###

#res_gse$log_padjust <- -log10(res_gse$p.adjust)

plot1 <- ggdotchart(res_gse, x = "Description", y = "log_padjust",
                    color = "gene_clusters",
                    palette = c("blue", "#FC4E07"),
                    sorting = "descending",
                    rotate = TRUE,
                    group = "gene_clusters",
                    dot.size = "setSize",
                    add = "segments",
                    title = "Enriched BP pathways in myonuclei cells (old cohort)",
                    xlab = "Biologic Pathways",
                    ylab = "-log10(padjust)",
                    ggtheme = theme_pubr()
)
plot1 + theme(legend.position = "right")+
  theme_bw()



