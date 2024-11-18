rm(list = ls())
require(devil)
require(glmGamPoi)
require(nebula)
require(TENxPBMCData)
require(foreach)
require(future)
require(rngtools)
require(magrittr)
source('utils.R')
set.seed(123456)

MIN_ITER = 1

for (n_cond in c(1,2,4)) {
  dds <- readRDS(paste0("datasets/",n_cond,"_cond.RDS"))

  cnt = as.matrix(dds$countData)
  design_matrix = dds$design

  n.genes <- dim(cnt)[1]
  n.cells <- dim(cnt)[2]

  d.timings <- NULL

  for (p_genes in c(.25, 1)) {
    print(paste0("p genes = ", p_genes))
    n.sub.genes <- as.integer(p_genes * n.genes)
    sub.genes.idx <- sample(1:n.genes, n.sub.genes)

    for (p_cells in c(.25,1)) {
      print(paste0("p cells = ", p_cells))
      n.sub.cell <- as.integer(p_cells * n.cells)
      sub.cell.idx <- sample(1:n.cells, n.sub.cell)

      print("subsetting data ...")

      c <- cnt[sub.genes.idx,sub.cell.idx]
      if (n_cond == 1) {
        d <- matrix(1, nrow = n.sub.cell)
      } else {
        d <- design_matrix[sub.cell.idx, ]
      }

      ids <- rep(1, n.sub.cell)

      print("inference starting ...")

      b.devil <- bench::mark(devil_res <<- my_fit_devil(
        c,
        d,
        overdispersion = FALSE,
        init_overdispersion = TRUE,
        offset=0,
        size_factors=TRUE,
        verbose=FALSE,
        max_iter=100,
        tolerance=1e-3,
        eps=0,
        CUDA = FALSE,
        batch_size = 2**13,
        parallel.cores=1
      ), min_iterations = MIN_ITER, memory = T) %>%
        dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'devil', p.genes=p_genes, p.cell=p_cells, n.cell=n.sub.cell, n.genes = n.sub.genes)

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
      ), min_iterations = MIN_ITER, memory = T) %>%
        dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'devil_gpu', p.genes=p_genes, p.cell=p_cells, n.cell=n.sub.cell, n.genes = n.sub.genes)

      if (n.sub.cell <= 150000) {
        b.glmGamPoi <- bench::mark(res_glm <<- glmGamPoi::glm_gp(c, d, size_factors = 'normed_sum', overdispersion = F), min_iterations = MIN_ITER, memory = T) %>%
          dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'glmGamPoi', p.genes=p_genes, p.cell=p_cells, n.cell=n.sub.cell, n.genes = n.sub.genes)

        b.total <- dplyr::bind_rows(b.devil, b.devil.gpu, b.glmGamPoi) %>% #, b.nebula) %>%
          dplyr::mutate(min=as.numeric(min, units='secs'), median=as.numeric(median, units='secs'), total_time=as.numeric(total_time, units='secs'))

        good_idx <- ((res_glm$Beta >= -100 & res_glm$Beta <= 100) %>% rowSums()) == ncol(d)

        print("Correlation between CPU devil and GLM:")
        print(cor.test(devil_res$beta[good_idx,], res_glm$Beta[good_idx,])$estimate)

        print("Correlation between GPU devil and GLM:")
        print(cor.test(devil_gpu_res$beta[good_idx,], res_glm$Beta[good_idx,])$estimate)
      } else {
        b.total <- dplyr::bind_rows(b.devil, b.devil.gpu) %>% #, b.nebula) %>%
          dplyr::mutate(min=as.numeric(min, units='secs'), median=as.numeric(median, units='secs'), total_time=as.numeric(total_time, units='secs'))
      }

      if (is.null(d.timings)) {
        d.timings <- b.total
      } else {
        d.timings <- dplyr::bind_rows(d.timings, b.total)
      }
      saveRDS(d.timings, paste0("results/", n_cond, "_timing_beta_only.rds"))
    }
  }
  saveRDS(d.timings, paste0("results/", n_cond, "_timing_beta_only.rds"))
}
