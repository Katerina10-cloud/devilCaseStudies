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

#install.packages("Signac")
library(Signac)

atac_seurat <- CreateSeuratObject(
  counts = peak_counts,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata_atac)


## Create a Gene Activity matrix ##
gene.activities <- GeneActivity(atac_seurat)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
atac_seurat[['RNA']] <- CreateAssayObject(counts = gene.activities)

atac_seurat <- NormalizeData(
  object = atac_seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac_seurat$nCount_RNA)
)

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

#snrna_muscl.RDS
#metadata_snrna_muscl.Rdata

### Peaks Anotation ###

#BiocManager::install("ChIPseeker")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("EnsDb.Hsapiens.v75")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

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
annot <- annot %>% mutate(ranges = paste(annot$start, annot$end, sep = "-"))
annot <- annot %>% mutate(ranges = paste(annot$seqnames, annot$ranges))

save(annot, file = "grange_annot.Rdata")

sum(res1$lfc > 0.2 & res1$padj < 0.05, na.rm=TRUE)


# ATAC data #
metadata_rownames <- as.data.frame(atac@colData@rownames)
colnames(metadata_rownames) <- "cellID"

sample <- as.data.frame(atac@colData@listData[["Sample"]])
colnames(sample) <- "sample"

patient <- as.data.frame(atac@colData@listData[["sample"]])
colnames(patient) <- "patient"

group <- as.data.frame(atac@colData@listData[["group"]])
colnames(group) <- "group"

umap1 <- as.data.frame(atac@colData@listData[["UMAP_1"]])
colnames(umap1) <- "UMAP_1"

umap2 <- as.data.frame(atac@colData@listData[["UMAP_2"]])
colnames(umap2) <- "UMAP_2"

age <- as.data.frame(atac@colData@listData[["age"]])
colnames(age) <- "age"

annot <- as.data.frame(atac@colData@listData[["Annotation"]])
colnames(annot) <- "cell_type"

metadata_atac <- cbind(metadata_rownames, sample, patient, group, umap1, umap2, age, annot)
rownames(metadata_atac) <- metadata_atac$cellID

save(metadata_atac, file = "metadata_atac.Rdata")

# RNA data #

annot <- do.call("cbind", metadata_rna_muscl[["Annotation"]]) %>% as.data.frame()
colnames(annot) = c("cell_type", "cell_cluster")

cellID <- as.data.frame(metadata_rna_muscl[["_index"]])
names(cellID) = "cellID"

age <- as.data.frame(metadata_rna_muscl[["age"]])
names(age) = "age"

age_pop <- do.call("cbind", metadata_rna_muscl[["age_pop"]]) %>% as.data.frame()
colnames(age_pop) = c("age_pop", "age_group")

nCount_RNA <- as.data.frame(metadata_rna_muscl[["nCount_RNA"]])
names(nCount_RNA) = "nCount_RNA"

nFeature_RNA <- as.data.frame(metadata_rna_muscl[["nFeature_RNA"]])
names(nFeature_RNA) = "nFeature_RNA"

ident <- do.call("cbind", metadata_rna_muscl[["orig.ident"]]) %>% as.data.frame()
ident <- ident %>% select("categories")
names(ident) = "ident"

percent.mt <- as.data.frame(metadata_rna_muscl[["percent.mt"]])
names(percent.mt) <- "percent.mt"

sample <- do.call("cbind", metadata_rna_muscl[["sample"]]) %>% as.data.frame()
sample <- sample %>% select("categories")
names(sample) <- "sample"

tech <- do.call("cbind", metadata_rna_muscl[["tech"]]) %>% as.data.frame()
colnames(tech) <- c("tech", "tech_id")

metadata_rna <- cbind(cellID, age, age_pop, annot,ident,sample,tech,nCount_RNA,nFeature_RNA,percent.mt)
rownames(metadata_rna) <- metadata_rna$cellID

save(metadata_rna, file = "metadata_rna_muscl.Rdata")

# Plot Umap #
#install.packages("viridis")
library(viridis)

umap_rna <- umap_rna[ rownames(umap_rna) %in% rownames(metadata),]

p_umap_rna <- umap_rna %>% ggplot(aes(x = umap_rna$umap_1, y = umap_rna$umap_2, color = metadata$cell_type))+
  geom_point()+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snRNA (~ 212 000 cells)")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_minimal() +
  #theme(legend.position = "bottom")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)
p_umap_rna + scale_color_brewer(palette = "Spectral")
p_umap_rna + scale_color_viridis()
#Seurat::LabelClusters(plot =p_umap_rna, id=metadata$cell_cluster, color="black")
p_umap_rna

p_umap_atac <- metadata_atac %>% ggplot(aes(x = metadata_atac$UMAP_1, y = metadata_atac$UMAP_2, color = metadata_atac$cell_type))+
  geom_point()+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snATAC (~ 90 000 cells)")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_minimal() +
  #theme(legend.position = "bottom")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)

p_umap_rna + p_umap_atac
corr_plot

#Vienn diagram

# load Venn diagram package 
library("VennDiagram") 

# create pairwise Venn diagram 
overlap <- draw.pairwise.venn(area1=1507, area2=1034,cross.area=340, 
                   category=c("snATAC","snRNA"),fill=c("Red","Blue"))
overlap

# Select DE genes #
sum(res_atac$adj_pval < 0.05 & res_atac$lfc > 0.2, na.rm=TRUE) #up_regulated
sum(res_atac$adj_pval < 0.05 & res_atac$lfc < -0.2, na.rm=TRUE) #down-reg

atac_up <- subset(res_atac, res_atac$adj_pval < 0.05 & res_atac$lfc > 0.2)
atac_down <- subset(res_atac, res_atac$adj_pval < 0.05 & res_atac$lfc < -0.2)
atac_deg <- rbind(atac_up, atac_down)

rna_up <- subset(res, res$adj_pval < 0.05 & res$lfc > 0.2)
rna_down <- subset(res, res$adj_pval < 0.05 & res$lfc < -0.2)
rna_deg <- rbind(rna_up, rna_down)

granAnn <- genRang_annot %>% select(SYMBOL, ranges)
colnames(granAnn)[2] <- "name"
atac_deg <- merge(atac_deg, granAnn, by = "name")
colnames(atac_deg)[5] <- "geneID"
colnames(rna_deg)[1] <- "geneID"

atac_deg_nodup <- atac_deg[!duplicated(atac_deg$geneID), ]
atac_deg_nodup <- atac_deg_nodup[ atac_deg_nodup$geneID %in% rna_deg$geneID,]
rna_deg <- rna_deg[ rna_deg$geneID %in% atac_deg_nodup$geneID,]

atac_up <- subset(atac_deg_nodup, atac_deg_nodup$adj_pval < 0.05 & atac_deg_nodup$lfc >= 1)
atac_down <- subset(atac_deg_nodup, atac_deg_nodup$adj_pval < 0.05 & atac_deg_nodup$lfc <= -1.0)
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

### Correlation analysis ###
#install.packages('smplot2')
library(smplot2)

corr_plot <- ggplot(mapping = aes(x = res_atac_rna$lfc_snRNA, y = res_atac_rna$lfc_snATAC)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  labs(title = "Correlation of significantly DAG (snATAC-seq) and GE (snRNA-seq)")+
  xlab("LogFC2 snRNA") +
  ylab ("LogFC2 snATAC") +
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  sm_statCorr()


library(ggplot2)
cor_coefs <- cor.test(res_atac_rna$lfc_snRNA, res_atac_rna$lfc_snATAC)
corr_plot <- ggplot(data = res_atac_rna, aes(x = lfc_snRNA, y = lfc_snATAC)) + 
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  annotate("text", x = 10, y = 4, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 10, y = 3.5, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))
corr_plot

dim_red_rna <- readRDS("dim_red_rna.RDS")
umap_rna <- as.data.frame(dim_red_rna[["umap"]]@cell.embeddings)



