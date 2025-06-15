
rm(list = ls())

UTILS_DIR = "../de_analysis/nullpower/utils/"
DATA_DIR = "../de_analysis/nullpower/datasets/"

source(file.path(UTILS_DIR, "utils.R"))
source(file.path(UTILS_DIR,"edgeR.R"))
source(file.path(UTILS_DIR,"limma.R"))
source(file.path(UTILS_DIR,"glmGamPoi.R"))
source(file.path(UTILS_DIR,"nebula.R"))
source(file.path(UTILS_DIR,"devil.R"))
library(Seurat)
library(ggplot2)
library(tidyverse)
library(muscData)

list.func <- list(
  glmgp.mult,
  edger.mult,
  limma.mult,
  glmgp.cell.mult,
  nebula.mult,
  devil.base,
  devil.mixed,
  glmgp.cell.fixed
)


data = muscData::Crowell19_4vs4(metadata = TRUE)
seurat.obj <- data[["EH3297"]]

# Keep only controls
seurat.obj[group_id == "Vehicle",]
#seurat.obj = readRDS("../cell_types_analysis/datasets/pbmc.rds"))
seurat.obj@metadata

cell_types = seurat.obj@meta.data %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(-n) %>% 
  slice_head(n = 5) %>% 
  dplyr::pull(cell_type)

ct.idx <- 1
ct <- cell_types[ct.idx]


# Subset Seurat object to selected cell type
# Filter by cell type: keep only cells where cell_type == ct
cell_ids <- rownames(seurat.obj@meta.data)[seurat.obj@meta.data$cell_type == ct]
sub.seurat.obj <- seurat.obj[, cell_ids]

# Get unique donors
unique_donors <- unique(sub.seurat.obj@meta.data$donor_id)

# Randomly assign donors to groups
set.seed(123)  # for reproducibility
donor_assignment <- data.frame(
  donor_id = unique_donors,
  group = sample(rep(c(0, 1), length.out = length(unique_donors)))
)

# Assign group info to metadata
sub.seurat.obj@meta.data$tx_cell <- donor_assignment$group[
  match(sub.seurat.obj@meta.data$donor_id, donor_assignment$donor_id)
]

# Optional: check how many cells per group
table(sub.seurat.obj@meta.data$tx_cell)

cnt.input = as.matrix(sub.seurat.obj@assays$RNA@counts)
sub.seurat.obj@meta.data$id = sub.seurat.obj@meta.data$donor_id
col.data.select = sub.seurat.obj@meta.data

# Filter genes
cnt.input = cnt.input[rowMeans(cnt.input) > .01,]

list.result.method <- list()
timings <- c()

list.func[[7]]

s <- Sys.time()

for (int.test in 1:length(list.func)){
  print(int.test)
  list.result.method[[int.test]] <- list.func[[int.test]](
    cnt.input, 
    col.data.select
  )[,4:5]
  
  timings <- c(timings, unique(list.result.method[[int.test]][,2]))
  list.result.method[[int.test]] <- list.result.method[[int.test]][,1]
}
e <- Sys.time()



df.result <- do.call(rbind, list.result.method)
rownames(df.result) <- c('glmGamPoi (Pb)', 'edgeR (Pb)', 'limma (Pb)', 'glmGamPoi (cell)', 'Nebula', 'Devil (base)', 'Devil (mixed)', "glmGamPoi (fixed)")

lapply(1:nrow(df.result), function)

df.result[1,] %>% hist()
df.result[2,] %>% hist()
df.result[3,] %>% hist()
df.result[4,] %>% hist()
df.result[5,] %>% hist()
df.result[6,] %>% hist()
df.result[7,] %>% hist()
df.result[8,] %>% hist()




design_matrix <- model.matrix(~1+tx_cell, data = df)
clusters = as.factor(df$id)

s <- Sys.time()
fit <- devil::fit_devil(count, design_matrix, size_factors=FALSE, verbose=F, parallel.cores=1, init_overdispersion = 100, offset = 1e-6, max_iter = 200, tolerance = 1e-3)
#fit <- devil::fit_devil(count, design_matrix, size_factors=TRUE, verbose=F, parallel.cores=1, init_overdispersion = 100, offset = 1e-6, max_iter = 200, tolerance = 1e-3)
# fit <- devil::fit_devil(count, design_matrix, size_factors=T, verbose=T, min_cells=-1, avg_counts=-1, parallel.cores=4)
e <- Sys.time()
fit$input_parameters$parallel = FALSE
delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters)

beta <- fit$beta[2,]
pval <- test$pval
tval <- qnorm(1-pval/2) * sign(beta)
se <- beta/tval

result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
result$delta_time <- delta_time
colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
return(result %>% as.matrix())


if (dataset_name == "Crowell19") {
  data <- muscData::Crowell19_4vs4(metadata = TRUE)
  sce <- data[["EH3297"]]
  sce$donor_id <- sce$sample_id
  sce.sub <- sce[, sce$group_id == "Vehicle"]
  cnt_mat = sce.sub@assays@data$counts
  meta = sce.sub@colData
} else if (dataset_name == "BacherTCellData") {
  sce = scRNAseq::BacherTCellData()
  sce$donor_id <- sce$sample_id
  sce.sub <- sce[, sce$diagnosis == "Healthy"]
  cnt_mat = sce.sub@assays@data$counts
  meta = sce.sub@colData
  all(colnames(cnt_mat) == rownames(meta))
} else if (dataset_name == "HuCortexData") {
  sce = scRNAseq::HuCortexData(mode = c("ctx"))
  sce$donor_id <- sce$Sample
  sce.sub <- sce[, sce$Mode == "ctx"]
  cnt_mat = sce.sub@assays@data$counts
  meta = sce.sub@colData
  all(colnames(cnt_mat) == rownames(meta))
} else if (dataset_name == "LedergorMyelomaData") {
  # LedergorMyelomaData
  sce = scRNAseq::LedergorMyelomaData()
  sce$donor_id <- sce$Subject_ID
  sce.sub <- sce[, sce$Condition == "Control"]
  cnt_mat = sce.sub@assays@data$counts
  meta = sce.sub@colData
  all(colnames(cnt_mat) == rownames(meta))
} else if (dataset_name == "GrunHSCData") {
  # GrunHSCData
  sce = scRNAseq::GrunHSCData()
  sce$donor_id <- sce$sample
  sce.sub <- sce[, sce$protocol == "sorted hematopoietic stem cells"]
  cnt_mat = sce.sub@assays@data$counts
  meta = sce.sub@colData
  all(colnames(cnt_mat) == rownames(meta))
} else {
  stop(paste("Unknown dataset:", dataset_name))
}






# LedergorMyelomaData
sce = scRNAseq::ZilionisLungData(c("mouse"))
sce$donor_id <- sce$Subject_ID
sce.sub <- sce[, sce$Condition == "Control"]
cnt_mat = sce.sub@assays@data$counts
meta = sce.sub@colData
all(colnames(cnt_mat) == rownames(meta))