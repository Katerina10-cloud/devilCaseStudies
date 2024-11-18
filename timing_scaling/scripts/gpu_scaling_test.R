rm(list = ls())
require(devil)
require(magrittr)
source('utils.R')
set.seed(123456)

MIN_ITER = 1
N_CELL_TYPES = 2

data = prep_data(N_CELL_TYPES = N_CELL_TYPES)
cnt <- data$cnt
design_matrix <- data$design_matrix

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

for (p_gene in c(.1, .5, 1)) {
  print(paste0("p genes = ", p_genes))
  n.sub.genes <- as.integer(p_genes * n.genes)
  sub.genes.idx <- sample(1:n.genes, n.sub.genes)

  for (p_cells in c(.1, .25, .5, 1)) {
    print(paste0("p cells = ", p_cells))
    n.sub.cell <- as.integer(p_cells * n.cells)
    sub.cell.idx <- sample(1:n.cells, n.sub.cell)

    c <- cnt[sub.genes.idx,sub.cell.idx]
    d <- design_matrix[sub.cell.idx, ]

    print("inference starting ...")

    b.devil.gpu <- bench::mark(devil_gpu_res <<- my_fit_devil(
      c,
      d,
      overdispersion = FALSE,
      init_overdispersion = FALSE,
      offset=0,
      size_factors=TRUE,
      verbose=FALSE,
      max_iter=100,
      tolerance=1e-3,
      eps=0,
      CUDA = TRUE,
      batch_size = 512,
      parallel.cores=1
    ), min_iterations = MIN_ITER, memory = T)

    res_name = paste0("results/macaque_brain/gpu_devil_", p_gene, "_pgene_", p_cells, "_pcells.rds")
    saveRDS(b.devil.gpu, file = res_name)
  }
}
