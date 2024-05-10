###-------------------------------------------------------------###
### scATACseq data analysis ###
###-------------------------------------------------------------###

# File conversion from AnnData/H5AD to h5Seurat #

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)

# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

Convert("scrna_muscl.h5ad", dest = "scrna_muscl_2.h5seurat", overwrite = TRUE)
SeuratObject <- LoadH5Seurat("scrna_muscl_2.h5seurat", misc = FALSE, meta.data = FALSE)

#setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/multiomics")
# atac_muscl.RDS  scrna_muscl.h5ad

BiocManager::install("rhdf5")
library(rhdf5)

# to see the structure of the file
h5ls("my_obj.h5ad")
#extract metadata
metadata <- h5read("scrna_muscl.h5ad", "/obs/")
class(metadata)
metadata <- do.call("cbind", metadata) %>% as.data.frame()
