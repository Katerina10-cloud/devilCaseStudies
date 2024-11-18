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

for (n_cond in c(4,2,1)) {
  dds <- readRDS(paste0("datasets/",n_cond,"_cond.RDS"))

  cnt = as.matrix(dds$countData)
  design_matrix = dds$design

  n.genes <- dim(cnt)[1]
  n.cells <- dim(cnt)[2]

  d.timings <- NULL

  for (p_genes in c(.1,.25,.5,.75,1)) {
    print(paste0("p genes = ", p_genes))
    n.sub.genes <- as.integer(p_genes * n.genes)
    sub.genes.idx <- sample(1:n.genes, n.sub.genes)

    for (p_cells in c(.1,.25,.5,.75,1)) {
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

      #sf <- devil:::calculate_sf(c)

      print("inference starting ...")

      # print(class(c))
      # print(c[1,1])
      # print(class(c[1,1]))
      # print(class(d))
      # print(class(d[1,1]))
      # print(d[1,1])

      # s <- Sys.time()
      # devil_res <- my_fit_devil(
      #   c,
      #   d,
      #   overdispersion = TRUE,
      #   offset=0,
      #   size_factors=TRUE,
      #   verbose=TRUE,
      #   max_iter=100,
      #   tolerance=1e-3,
      #   eps=0,
      #   CUDA = FALSE,
      #   batch_size = 2**13,
      #   parallel.cores=1
      # )
      # e <- Sys.time()
      # print(paste0("devil cpu time = ", as.numeric(e - s, units = 'secs')))
      #
      # s <- Sys.time()
      # gpu_devil_res <- my_fit_devil(
      #   c,
      #   d,
      #   overdispersion = TRUE,
      #   offset=0,
      #   size_factors=TRUE,
      #   verbose=TRUE,
      #   max_iter=100,
      #   #tolerance=1e-3,
      #   tolerance=1e-3,
      #   eps=0,
      #   CUDA = TRUE,
      #   #batch_size = dim(c)[1],
      #   batch_size = 1024,
      #   parallel.cores=1
      # )
      # e <- Sys.time()
      # print(paste0("devil GPU time = ", as.numeric(e - s, units = 'secs')))

      b.devil <- bench::mark(my_fit_devil(
        c,
        d,
        overdispersion = TRUE,
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

      b.devil.gpu <- bench::mark(my_fit_devil(
        c,
        d,
        overdispersion = TRUE,
        offset=0,
        size_factors=TRUE,
        verbose=FALSE,
        max_iter=100,
        tolerance=1e-3,
        eps=0,
        CUDA = TRUE,
        batch_size = 1024,
        parallel.cores=1
      ), min_iterations = MIN_ITER, memory = T) %>%
        dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'devil_gpu', p.genes=p_genes, p.cell=p_cells, n.cell=n.sub.cell, n.genes = n.sub.genes)


      b.glmGamPoi <- bench::mark(glmGamPoi::glm_gp(c, d, size_factors = 'normed_sum', overdispersion = T), min_iterations = MIN_ITER, memory = T) %>%
        dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'glmGamPoi', p.genes=p_genes, p.cell=p_cells, n.cell=n.sub.cell, n.genes = n.sub.genes)

      #b.nebula <- bench::mark(nebula::nebula(c, ids, d, cpc = -1, mincp = -1, ncore = 1, offset = sf), min_iterations = 3, memory = T) %>%
      #  dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'nebula', p=p, n.cell=n.sub.cell)

      b.total <- dplyr::bind_rows(b.devil, b.devil.gpu, b.glmGamPoi) %>% #, b.nebula) %>%
        dplyr::mutate(min=as.numeric(min, units='secs'), median=as.numeric(median, units='secs'), total_time=as.numeric(total_time, units='secs'))

      if (is.null(d.timings)) {
        d.timings <- b.total
      } else {
        d.timings <- dplyr::bind_rows(d.timings, b.total)
      }
      saveRDS(d.timings, paste0("results/", n_cond, "_timing.rds"))
    }
  }
  saveRDS(d.timings, paste0("results/", n_cond, "_timing.rds"))
}
#
# for (d_name in c("pbmc4k", "pbmc68k")) {
#   data <- TENxPBMCData(dataset = d_name)
#   cnt <- data@assays@data$counts
#
#   cnt <- cnt %>% as.matrix()
#   cnt <- cnt[rowMeans(cnt) > .01,]
#
#   n.genes <- dim(cnt)[1]
#   n.cells <- dim(cnt)[2]
#
#   d.timings <- NULL
#   for (p in c(.1,.25,.5,.75,1)) {
#     print(p)
#     n.sub.cell <- as.integer(p * n.cells)
#     sub.cell.idx <- sample(1:n.cells, n.sub.cell)
#
#     c <- cnt[,sub.cell.idx] %>% as.matrix()
#     d <- model.matrix(~1, dplyr::tibble(1:n.sub.cell))
#     d <- cbind(d, rnorm(nrow(d), 0, 1))
#     ids <- rep(1, n.sub.cell)
#
#     sf <- devil:::calculate_sf(c)
#
#     s <- Sys.time()
#     devil_res <- my_fit_devil(
#       c,
#       d,
#       overdispersion = TRUE,
#       offset=0,
#       size_factors=TRUE,
#       verbose=FALSE,
#       max_iter=200,
#       tolerance=1e-3,
#       eps=0,
#       CUDA = FALSE,
#       batch_size = 2**13,
#       parallel.cores=1
#     )
#     e <- Sys.time()
#     print(paste0("devil cpu time = ", e - s))
#
#     s <- Sys.time()
#     gpu_devil_res <- my_fit_devil(
#       c,
#       d,
#       overdispersion = TRUE,
#       offset=0,
#       size_factors=TRUE,
#       verbose=FALSE,
#       max_iter=200,
#       tolerance=1e-3,
#       eps=0,
#       CUDA = TRUE,
#       batch_size = dim(c)[1],
#       parallel.cores=1
#     )
#     e <- Sys.time()
#     print(paste0("devil cpu time = ", e - s))
#
#     b.devil <- bench::mark(my_fit_devil(
#       c,
#       d,
#       overdispersion = TRUE,
#       offset=0,
#       size_factors=TRUE,
#       verbose=FALSE,
#       max_iter=200,
#       tolerance=1e-3,
#       eps=0,
#       CUDA = FALSE,
#       batch_size = 2**13,
#       parallel.cores=1
#     ), min_iterations = 3, memory = T) %>%
#       dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'devil', p=p, n.cell=n.sub.cell)
#
#     b.devil.gpu <- bench::mark(my_fit_devil(
#       c,
#       d,
#       overdispersion = TRUE,
#       offset=0,
#       size_factors=TRUE,
#       verbose=FALSE,
#       max_iter=200,
#       tolerance=1e-3,
#       eps=0,
#       CUDA = TRUE,
#       batch_size = 2**13,
#       parallel.cores=1
#     ), min_iterations = 3, memory = T) %>%
#       dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'devil_gpu', p=p, n.cell=n.sub.cell)
#
#
#     b.glmGamPoi <- bench::mark(glmGamPoi::glm_gp(c, d, size_factors = 'normed_sum', overdispersion = T), min_iterations = 3, memory = T) %>%
#       dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'glmGamPoi', p=p, n.cell=n.sub.cell)
#
#     #b.nebula <- bench::mark(nebula::nebula(c, ids, d, cpc = -1, mincp = -1, ncore = 1, offset = sf), min_iterations = 3, memory = T) %>%
#     #  dplyr::select(min, median, n_itr, total_time, time, mem_alloc) %>% dplyr::mutate(algo = 'nebula', p=p, n.cell=n.sub.cell)
#
#     b.total <- dplyr::bind_rows(b.devil, b.devil.gpu, b.glmGamPoi) %>% #, b.nebula) %>%
#       dplyr::mutate(min=as.numeric(min, units='secs'), median=as.numeric(median, units='secs'), total_time=as.numeric(total_time, units='secs'))
#
#     if (is.null(d.timings)) {
#     	d.timings <- b.total
#     } else {
#     	d.timings <- dplyr::bind_rows(d.timings, b.total)
#     }
#     saveRDS(d.timings, paste0("results/", d_name, "_timing.rds"))
#   }
#
#   saveRDS(d.timings, paste0("results/", d_name, "_timing.rds"))
# }
#
#
#
#
