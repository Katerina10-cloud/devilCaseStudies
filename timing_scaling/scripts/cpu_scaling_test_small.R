rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source("scripts/utils.R")
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

    print("inference starting devil ...")

    b.devil <- bench::mark(devil::fit_devil(
      as.matrix(c),
      d,
      overdispersion = T,
      size_factors = T,
      verbose = T,
      parallel.cores = 1,
      offset = 1e-6
    ), min_iterations = MIN_ITER, memory = T)
    b.devil$result <- NULL

    res_name = paste0("results/baronPancreas/cpu_devil_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.devil, file = res_name)

    print("inference starting glm ...")

    b.glmGamPoi <- bench::mark(
      glmGamPoi::glm_gp(c, d, size_factors = 'normed_sum', overdispersion = T), min_iterations = MIN_ITER, memory = T
    )
    b.glmGamPoi$result <- NULL

    res_name = paste0("results/baronPancreas/cpu_glmGamPoi_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.glmGamPoi, file = res_name)
  }
}

# Save beta and overdispersion for a specific object
p_genes <- 1
p_cells <- 1

# print(paste0("p genes = ", p_genes))
# n.sub.genes <- as.integer(p_genes * n.genes)
# sub.genes.idx <- sample(1:n.genes, n.sub.genes)
# 
# print(paste0("p cells = ", p_cells))
# n.sub.cell <- as.integer(p_cells * n.cells)
# sub.cell.idx <- sample(1:n.cells, n.sub.cell)

c <- cnt
d <- design_matrix

print("Devil fitting...")

s <- Sys.time()
fit.devil <- devil::fit_devil(
  as.matrix(c),
  d,
  overdispersion = T,
  size_factors = T,
  verbose = T,
  parallel.cores = 1,
  offset = 1e-6
)
e <- Sys.time()
print(e-s)

devil.res <- devil::test_de(fit.devil, contrast = c(0,1))
devil.final.res <- devil.res %>%
  cbind(fit.devil$beta) %>%
  cbind(dplyr::tibble(theta = fit.devil$overdispersion))
rm(fit.devil, devil.res)

print("glmGamPoi fitting...")

s <- Sys.time()
fit.glm <- glmGamPoi::glm_gp(
  as.matrix(c),
  d,
  overdispersion = T,
  size_factors = "normed_sum",
  offset = 1e-6
)
e <- Sys.time()
print(e-s)

glm.res <- glmGamPoi::test_de(fit.glm, contrast = c(0,1))
glm.final.res <- glm.res %>%
  cbind(fit.glm$Beta) %>%
  cbind(dplyr::tibble(theta = fit.glm$overdispersions))
rm(glm.res, fit.glm)

print("Saving results...")

saveRDS(devil.final.res, paste0("results/baronPancreas/fits/cpu_devil_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds"))
saveRDS(glm.final.res, paste0("results/baronPancreas/fits/cpu_glmGamPoi_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds"))
