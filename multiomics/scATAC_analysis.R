###-------------------------------------------------------------###
### snATACseq + snRNA data analysis ###
###-------------------------------------------------------------###

### Extract scRNA data ###

# File conversion from AnnData/H5AD to h5Seurat #
#remotes::install_github("mojaveazure/seurat-disk")
# setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/multiomics")
# atac_muscl.RDS  scrna_muscl.h5ad

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/sc_devil/data/muscle")

library(SeuratDisk)
library(tidyverse)

# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
SeuratDisk::Convert("scrna_muscl.h5ad", dest = "scrna_muscl2.h5seurat", overwrite = TRUE)
SeuratObject <- SeuratDisk::LoadH5Seurat("scrna_muscl2.h5seurat", misc = FALSE, meta.data = FALSE)


#BiocManager::install("rhdf5")
library(rhdf5)

# to see the structure of the file
rhdf5::h5ls("my_obj.h5ad")
#extract metadata
metadata_rna_muscl <- rhdf5::h5read("scrna_muscl.h5ad", "/obs/")
class(metadata)

# add metadata to Seurat object
rna_seurat_muscl <- AddMetaData(rna_seurat_muscl, metadata = metadata_snrna_muscl$index, col.name = "index")


### Extract ATAC data ###
atac <- readRDS("atac_muscl.RDS")

metadata <- atac@colData@listData
metadata_atac <- cbind(sample, annotation, age, sex)
rownames(metadata_atac) <- colnames(peak_counts)

peak_counts <- atac@assays@data@listData[["PeakMatrix"]]
granges <- atac@rowRanges

#genomic_ranges <- genomic_ranges %>% mutate(ranges = paste(genomic_ranges$seqnames, genomic_ranges$ranges))
rownames(peak_counts) <- genomic_ranges$ranges

save(genomic_ranges, file = "genomic_ranges.Rdata")
save(peak_counts, file = "peak_counts.Rdata")
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
metadata_filtered <- metadata_atac[ (metadata_atac$cell_type %in% c("Type I", "Type II")), ]

metadata_filtered <- metadata_filtered %>% 
  mutate(cluster = case_when(
    cell_type == "Type I"  ~ '0',
    cell_type == "Type II"  ~ '1'))  

metadata_filtered$cluster <- as.factor(metadata_filtered$cluster)

peak_counts <- as.matrix(peak_counts)
peak_counts <- peak_counts[,colnames(peak_counts) %in% rownames(metadata_filtered)]

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




