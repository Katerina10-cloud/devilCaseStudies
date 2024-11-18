rm(list = ls())
require(devil)
require(magrittr)
source('utils.R')
set.seed(123456)

MIN_ITER = 3
N_CELL_TYPES = 2

data = prep_data_small()
cnt <- data$cnt
design_matrix <- data$design_matrix

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

for (p_genes in c(.1, .5, 1)) {
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

    b.devil.gpu <- bench::mark(my_fit_devil(
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
    b.devil.gpu$result <- NULL
    
    res_name = paste0("results/baronPancreas/gpu_devil_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.devil.gpu, file = res_name)
  }
}

# Save beta and overdispersion for a specific object
p_genes <- 1
p_cells <- 1

print(paste0("p genes = ", p_genes))
n.sub.genes <- as.integer(p_genes * n.genes)
sub.genes.idx <- sample(1:n.genes, n.sub.genes)

print(paste0("p cells = ", p_cells))
n.sub.cell <- as.integer(p_cells * n.cells)
sub.cell.idx <- sample(1:n.cells, n.sub.cell)

c <- cnt[sub.genes.idx,sub.cell.idx]
d <- design_matrix[sub.cell.idx, ]

print("Devil fitting...")

s <- Sys.time()
fit.devil <- my_fit_devil(
  c,
  d,
  overdispersion = TRUE,
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
)
e <- Sys.time()
print(e-s)

devil.res <- devil::test_de(fit.devil, contrast = c(0,1))
devil.final.res <- devil.res %>%
  cbind(fit.devil$beta) %>%
  cbind(dplyr::tibble(theta = fit.devil$overdispersion))
rm(fit.devil, devil.res)

saveRDS(devil.final.res, paste0("results/baronPancreas/fits/gpu_devil_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds"))
