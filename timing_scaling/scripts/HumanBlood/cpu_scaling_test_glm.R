rm(list = ls())
require(devil)
require(Seurat)
require(magrittr)
source('scripts/utils.R')
set.seed(123456)

MIN_ITER = 1
N_CELL_TYPES = 3
CONTRAST = c(0,0,0,1)

data = prep_HumanBlood_data(N_CELL_TYPES = N_CELL_TYPES)

cnt <- data$cnt
design_matrix <- data$design_matrix

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

print("Human blood dataset")
print(paste0(".  n genes = ", n.genes))
print(paste0(".  n cells = ", n.cells))

# Save beta and overdispersion for a specific object
print("Single fit test...")
n_genes <- 1000
n_cells <- 1e6

input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.cells = n_cells, n.sub.genes = n_genes)

print(paste0("N genes = ", n_genes))
print(paste0("N cells = ", n_cells))

c = as.matrix(input$c)
d = input$d

print("glmGamPoi fitting...")

s <- Sys.time()
fit.glm <- glmGamPoi::glm_gp(
  c,
  d,
  overdispersion = T,
  size_factors = "normed_sum",
  offset = 1e-6
)
e <- Sys.time()
print(e-s)

glm.res <- glmGamPoi::test_de(fit.glm, contrast = CONTRAST)
glm.final.res <- glm.res %>%
  cbind(fit.glm$Beta) %>%
  cbind(dplyr::tibble(theta = fit.glm$overdispersions))
rm(glm.res, fit.glm)

print("Saving results...")

saveRDS(glm.final.res, paste0("results/HumanBlood/fits/cpu_glmGamPoi_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds"))

print("Subsampling test...")

for (n_genes in c(100, 1000, 10000)) {
  for (n_cells in c(1000, 100000, 1e6)) {
    print(paste0("p genes = ", n_genes))
    print(paste0("p cells = ", n_cells))

    input <- filter_input(cnt, design_matrix, NULL, NULL, n.sub.genes = n_genes, n.sub.cells = n_cells)
    c = as.matrix(input$c)
    d = input$d

    print("inference starting glm ...")

    b.glmGamPoi <- bench::mark(
      glmGamPoi::glm_gp(c, d, size_factors = 'normed_sum', overdispersion = T), min_iterations = MIN_ITER, memory = T
    )
    b.glmGamPoi$result <- NULL

    res_name = paste0("results/HumanBlood/cpu_glmGamPoi_", n_genes, "_ngene_", n_cells, "_ncells_", N_CELL_TYPES, "_celltypes.rds")
    saveRDS(b.glmGamPoi, file = res_name)
  }
}


