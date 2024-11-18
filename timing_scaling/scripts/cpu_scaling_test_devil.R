rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source('scripts/utils.R')
set.seed(123456)

MIN_ITER = 1
N_CELL_TYPES = 2

data = prep_data(N_CELL_TYPES = N_CELL_TYPES)

cnt <- data$cnt
design_matrix <- data$design_matrix

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

print("Macaque brain dataset")
print(paste0(".  n genes = ", n.genes))
print(paste0(".  n cells = ", n.cells))

for (p_genes in c(.1, .5)) {
  print(paste0("p genes = ", p_genes))

  for (p_cells in c(.1, .25, .5, .75)) {
    print(paste0("p cells = ", p_cells))

    input <- filter_input(cnt, design_matrix, p_genes, p_cells)
    c = as.matrix(input$c)
    d = input$d

    print("inference starting devil ...")

    b.devil <- bench::mark(devil::fit_devil(
  		c,
  		d,
  		overdispersion = T,
  		size_factors = T,
  		verbose = T,
  		parallel.cores = 1,
  		offset = 1e-6
	  ), min_iterations = MIN_ITER, memory = T)
    b.devil$result <- NULL

    res_name = paste0("results/macaque_brain/cpu_devil_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.devil, file = res_name)
  }
}

# Save beta and overdispersion for a specific object
p_genes <- .25
p_cells <- .5

input <- filter_input(cnt, design_matrix, p_genes, p_cells)

print(paste0("p genes = ", p_genes))
print(paste0("p cells = ", p_cells))

c = as.matrix(input$c)
d = input$d

print("Devil fitting...")

s <- Sys.time()
fit.devil <- devil::fit_devil(
  c,
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

print("Saving results...")

saveRDS(glm.final.res, paste0("results/macaque_brain/fits/cpu_glmGamPoi_", p_genes, "_pgene_", p_cells, "_pcells_", N_CELL_TYPES, "_celltypes.rds"))
