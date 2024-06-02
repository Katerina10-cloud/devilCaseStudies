rm(list = ls())
require(devil)
require(glmGamPoi)
require(nebula)
require(TENxPBMCData)
require(foreach)
require(future)
require(rngtools)
require(magrittr)

set.seed(123456)

cnt <- Seurat::ReadMtx("dataset/filtered_feature_bc_matrix/matrix.mtx.gz", 
		       "dataset/filtered_feature_bc_matrix/barcodes.tsv.gz", 
		       "dataset/filtered_feature_bc_matrix/features.tsv.gz")

cnt <- cnt %>% as.matrix()
cnt <- cnt[rowMeans(cnt) > .001,]

n.genes <- dim(cnt)[1]
n.cells <- dim(cnt)[2]

print(dim(cnt))

d.timings <- dplyr::tibble()
for (p in c(.1,.25,.5,.75,1)) {
  print(p)
  n.sub.cell <- as.integer(p * n.cells)
  sub.cell.idx <- sample(1:n.cells, n.sub.cell)
  
  c <- cnt[,sub.cell.idx] %>% as.matrix()
  d <- model.matrix(~1, dplyr::tibble(1:n.sub.cell))
  ids <- rep(1, n.sub.cell)
  
  sf <- devil:::calculate_sf(c)
  
  b.devil <- bench::mark(devil::fit_devil(c, d, overdispersion = T, size_factors = T, avg_counts = -1, min_cells = -1, parallel.cores = 1), min_iterations = 3, memory = F) %>% 
    dplyr::select(min, median, n_itr, total_time, time) %>% dplyr::mutate(algo = 'devil', p=p, n.cell=n.sub.cell)
  b.glmGamPoi <- bench::mark(glmGamPoi::glm_gp(c, d, size_factors = 'normed_sum', overdispersion = T), min_iterations = 3, memory = F) %>% 
    dplyr::select(min, median, n_itr, total_time, time) %>% dplyr::mutate(algo = 'glmGamPoi', p=p, n.cell=n.sub.cell)
  b.nebula <- bench::mark(nebula::nebula(c, ids, d, cpc = -1, mincp = -1, ncore = 1, offset = sf), min_iterations = 3, memory = F) %>% 
    dplyr::select(min, median, n_itr, total_time, time) %>% dplyr::mutate(algo = 'nebula', p=p, n.cell=n.sub.cell)
  
  b.total <- dplyr::bind_rows(b.devil, b.glmGamPoi, b.nebula) %>% 
    dplyr::mutate(min=as.numeric(min, units='secs'), median=as.numeric(median, units='secs'), total_time=as.numeric(total_time, units='secs')) 
  
  # d <- model.matrix(~out-1, dplyr::tibble(out = rnorm(n.sub.cell, mean = 1, sd = .001)))
  # 
  # fit.devil <- devil::fit_devil(c, d, overdispersion = T, size_factors = T, avg_counts = -1, min_cells = -1, parallel.cores = 1)
  # fit.glmGamPoi <- glmGamPoi::glm_gp(c, d, size_factors = "normed_sum", overdispersion = T)
  # fit.nebula <- nebula::nebula(c, ids, d, cpc = -1, mincp = -1, ncore = 1, offset = sf)
  
  d.timings <- dplyr::bind_rows(d.timings, b.total)
  saveRDS(d.timings, "timings.rds") 
}

saveRDS(d.timings, "timings.rds")
