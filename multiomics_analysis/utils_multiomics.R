DATASET_NAMES <- c("MuscleRNA", "MuscleATAC")

read_data <- function(dataset_name) {
  if (dataset_name == "MuscleRNA") {
    seurat_data <- readRDS(data_path)
    counts <- Seurat::GetAssayData(object = seurat_data, layer = "counts")
    metadata <- seurat_data@meta.data
    tissue <- "muscle"
  } else if (dataset_name == "MuscleATAC") {
    atac <- readRDS(data_path)
    counts <- atac@assays@data@listData[["PeakMatrix"]]
    granges <- atac@rowRanges
    granges %>% mutate(ranges = paste(granges$seqname,granges$start,granges$end, sep = ":"))
    rownames(peak_counts) <- granges$ranges
    metadata <- atac@colData %>% as.data.frame()
    tissue <- "muscle"
  }  else {
    stop("Dataset name not recognized")
  }
  
  return(list(counts=counts, metadata=metadata, grange=granges, tissue=tissue))
}


grange_annot <- function(input_data) {
  counts <- input_data$counts
  matadata <- input_data$metadata
  grange <- input_data$grange
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  grange_annot <- ChIPseeker::annotatePeak(
    peak = grange,
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
  tissue = "muscle"
  return(list(counts=counts, metadata=metadata, grange=grange_annot, tissue=tissue))
}

