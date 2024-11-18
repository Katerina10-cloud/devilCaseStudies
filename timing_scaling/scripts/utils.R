
prep_data <- function(N_CELL_TYPES, min_count_per_cell = 100, min_cell_per_gene = 100) {
  data = readRDS("datasets/macaque_brain.rds")

  metadata <- data@meta.data
  cell_types = names(sort(table(metadata$cell_type), decreasing = TRUE))[1:N_CELL_TYPES]

  print("Using the following cell types:")
  print(cell_types)

  cell_idx = which(metadata$cell_type %in% cell_types == TRUE)

  metadata <- metadata[cell_idx,]
  metadata$cell_type <- factor(metadata$cell_type, levels = unique(metadata$cell_type))
  #cnt <- as.matrix(data@assays$RNA@counts)[,cell_idx]
  cnt <- data@assays$RNA@counts[,cell_idx]

  cell_idx <- which((colSums(cnt) > min_count_per_cell) == TRUE)
  cnt <- cnt[,cell_idx]
  
  gene_idx = which((rowSums(cnt) > min_cell_per_gene) == TRUE)
  cnt <- cnt[gene_idx,]
  
  metadata <- metadata[cell_idx,]
  design_matrix <- model.matrix(~cell_type, metadata)

  return(list(cnt = cnt, design_matrix = design_matrix))
}

filter_input <- function(cnt, design_matrix, p_genes, p_cells) {
  n.genes <- dim(cnt)[1]
  n.cells <- dim(cnt)[2]
  
  n.sub.genes <- as.integer(p_genes * n.genes)
  n.sub.cells <- as.integer(p_cells * n.cells)
  
  # Find cells with highest expression
  cell.idxs <- order(colSums(cnt), decreasing = TRUE)[1:n.sub.cells]
  cnt <- cnt[,cell.idxs]
  design_matrix <- design_matrix[cell.idxs,]
  
  # Find genes with highest expression
  gene.idxs <- order(rowSums(cnt), decreasing = TRUE)[1:n.sub.genes]
  cnt <- cnt[gene.idxs,]
  
  print(paste0("Lowest RNA counts for a gene = ", min(rowSums(cnt))))
  list(c=cnt, d = design_matrix)
}

prep_data_small <- function(min_count_per_cell = 100, min_cell_per_gene = 100) {
  data = readRDS("datasets/baronPancreas.rds")
  cnt <- as.matrix(data$counts)
  metadata <- data$metadata
  cell_idx <- (metadata$label %in% c("alpha", "beta"))
  
  cnt <- cnt[, cell_idx]
  cnt <- cnt[rowSums(cnt > 0) > 100, ]
  metadata <- metadata[cell_idx,]
  
  design_matrix <- model.matrix(~label, data = metadata)
  
  return(list(cnt = cnt, design_matrix = design_matrix))
}

my_fit_devil <- function(
    input_matrix,
    design_matrix,
    overdispersion = TRUE,
    init_overdispersion = TRUE,
    offset=1e-6,
    size_factors=TRUE,
    verbose=FALSE,
    max_iter=200,
    tolerance=1e-3,
    #eps=1e-6,
    CUDA = TRUE,
    batch_size = 1024L,
    parallel.cores=NULL) {

  max.cores <- parallel::detectCores()
  if (is.null(parallel.cores)) {
    n.cores = max.cores
  } else {
    if (parallel.cores > max.cores) {
      message(paste0("Requested ", parallel.cores, " cores, but only ", max.cores, " available."))
    }
    n.cores = min(max.cores, parallel.cores)
  }

  # Check if CUDA is available
  CUDA_is_available <- FALSE
  if (CUDA) {
    message("Check CUDA availability function need to be implemented")
    CUDA_is_available <- TRUE
  }

  # Add epsilon to input_matrix to avoid non invertible matrices
  # start_time <- Sys.time()
  # 
  # input_matrix <- input_matrix + eps
  # input_mat <- devil:::handle_input_matrix(input_matrix, verbose=verbose)
  # 
  # end_time <- Sys.time()
  # message("INPUT MATRIX HANDLING:")
  # message(as.numeric(difftime(end_time, start_time, units = "secs")))

  gene_names <- rownames(input_matrix)
  # counts_per_cell <- rowMeans(input_mat)
  # cell_per_genes <- rowSums(input_mat > 0)
  # filter_genes <- (counts_per_cell <= avg_counts) | (cell_per_genes <= min_cells)
  # n_low_genes <- sum(filter_genes)
  # if (n_low_genes > 0) {
  #   message(paste0("Removing ", n_low_genes, " lowly expressed genes."))
  #   input_mat <- matrix(input_mat[!filter_genes, ], ncol = nrow(design_matrix), nrow = sum(!filter_genes))
  #   input_matrix <- input_matrix[!filter_genes, ]
  #   gene_names <- gene_names[!filter_genes]
  # }

  start_time <- Sys.time()
  if (size_factors) {
    if (verbose) { message("Compute size factors") }
    sf <- devil:::calculate_sf(input_matrix, verbose = verbose)
  } else {
    sf <- rep(1, nrow(design_matrix))
  }
  end_time <- Sys.time()
  message("Size factor computing:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  start_time <- Sys.time()
  offset_matrix = devil:::compute_offset_matrix(offset, input_matrix, sf)
  #offset_matrix = devil:::compute_offset_matrix(offset, input_matrix, sf)[1,]
  #offset_matrix = matrix(offset_matrix, nrow = 1)
  end_time <- Sys.time()
  message("Offset matrix computing:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))
  
  estimate_dispersion <- function (y, offset_matrix) {
    xim <- 1/mean(DelayedMatrixStats::colMeans2(exp(offset_matrix),
                                                useNames = T))
    bv <- DelayedMatrixStats::rowVars(y, useNames = T)
    bm <- DelayedMatrixStats::rowMeans2(y, useNames = T)
    disp <- (bv - xim * bm)/bm^2
    ifelse(is.na(disp) | disp < 0, 100, disp)
  }

  start_time <- Sys.time()
  if (init_overdispersion) {
    dispersion_init <- c(estimate_dispersion(input_matrix, offset_matrix))
  } else {
    dispersion_init <- rep(100, nrow(input_mat))
  }

  end_time <- Sys.time()
  message("Dispersion init:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  ngenes <- nrow(input_matrix)
  nfeatures <- ncol(design_matrix)
  
  

  start_time <- Sys.time()
  if (verbose) { message("Initialize beta estimate") }
  groups <- devil:::get_groups_for_model_matrix(design_matrix)
  end_time <- Sys.time()
  message("get_groups_for_model_matrix:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))
  
  init_beta <- function (y, design_matrix, offset_matrix) {
    qrx <- qr(design_matrix)
    Q <- qr.Q(qrx)[seq_len(nrow(design_matrix)), , drop = FALSE]
    R <- qr.R(qrx)
    norm_log_count_mat <- t(log1p((y/exp(offset_matrix[1, ]))))
    t(solve(R, as.matrix(t(Q) %*% norm_log_count_mat)))
  }
  
  init_beta_groups <- function(y, groups, offset_matrix) {
    #norm_Y <- y / exp(offset_matrix)
    norm_Y <- y / exp(offset_matrix[1,])
    do.call(cbind, lapply(unique(groups), function(gr){
      log(DelayedMatrixStats::rowMeans2(norm_Y, cols = groups == gr, useNames=TRUE))
    }))
  }

  start_time <- Sys.time()
  # if (!is.null(groups)) {
  #   beta_0 <- devil:::init_beta_groups(input_matrix, groups, offset_matrix)
  #   #beta_0 <- devil:::init_beta(input_mat, design_matrix, offset_matrix)
  # } else {
  #   beta_0 <- devil:::init_beta(input_matrix, design_matrix, offset_matrix)
  # }
  beta_0 <- init_beta(input_matrix, design_matrix, offset_matrix)
  #beta_0 <- init_beta_groups(input_matrix, groups, offset_matrix)
  end_time <- Sys.time()
  message("init_beta:")
  message(as.numeric(difftime(end_time, start_time, units = "secs")))

  #beta_0 <- devil:::init_beta(input_mat, design_matrix, offset_matrix)

  if (CUDA & CUDA_is_available) {
    message("Messing with CUDA! Implementation still needed")

    message("Fit beta CUDA")

    start_time <- Sys.time()
    res_beta_fit <- devil:::beta_fit_gpu(input_matrix, design_matrix, beta_0, offset_matrix, dispersion_init, max_iter = max_iter, eps = tolerance, batch_size = batch_size)
    end_time <- Sys.time()
    message("BETA GPU RUNTIME:")
    message(as.numeric(difftime(end_time, start_time, units = "secs")))

    beta = res_beta_fit$mu_beta
    beta <- beta[1:ngenes,]
    if (is.null(dim(beta))) {
      beta = matrix(beta, ncol = 1)
    }
    
    iterations=res_beta_fit$iter

  } else {
    start_time <- Sys.time()
    if (verbose) { message("Fit beta coefficients") }
    tmp <- parallel::mclapply(1:ngenes, function(i) {
      devil:::beta_fit(input_mat[i,], design_matrix, beta_0[i,], offset_matrix[1,], dispersion_init[i], max_iter = max_iter, eps = tolerance)
      #devil:::beta_fit(input_mat[i,], design_matrix, beta_0[i,], offset_matrix[i,], 1, max_iter = max_iter, eps = tolerance)
    }, mc.cores = n.cores)
    beta <- lapply(1:ngenes, function(i) { tmp[[i]]$mu_beta }) %>% do.call("rbind", .)
    rownames(beta) <- gene_names
    end_time <- Sys.time()
    message("BETA CPU RUNTIME:")
    message(as.numeric(difftime(end_time, start_time, units = "secs")))

    iterations <- lapply(1:ngenes, function(i) { tmp[[i]]$iter }) %>% unlist()
  }


  if (overdispersion) {
    if (verbose) { message("Fit overdispersion") }

    theta <- parallel::mclapply(1:ngenes, function(i) {
      devil:::fit_dispersion(beta[i,], design_matrix, input_matrix[i,], offset_matrix[i,], tolerance = tolerance, max_iter = max_iter)
    }, mc.cores = n.cores) %>% unlist()

  } else {
    theta = rep(0, ngenes)
  }

  return(
    list(
      beta=beta,
      overdispersion=theta,
      iterations=iterations,
      size_factors=sf,
      offset_matrix=offset_matrix,
      design_matrix=design_matrix,
      input_matrix=input_matrix,
      input_parameters=list(max_iter=max_iter, tolerance=tolerance, parallel.cores=n.cores)
    )
  )
}


make_example_dds <- function(n_genes, n_replicates, n_conditions){

  dispMeanRel <-  function(x) 4/x + 0.1

  beta <- matrix(rnorm(n_genes*n_conditions, mean = 4, sd = 2), ncol = n_conditions)
  dispersion <- dispMeanRel(2^(beta[, 1]))
  colData <- data.frame(condition = factor(rep(paste0("cond", 1:n_conditions), n_replicates)))
  if (n_conditions == 1) {
    x <- matrix(data = 1, nrow = n_replicates, ncol = 1)
  } else {
    x <- model.matrix.default(~colData$condition)
  }

  mu <- t(2^(x %*% t(beta)))

  countData <- matrix(rnbinom(mu, mu = mu, size = 1/dispersion), ncol = ncol(mu))
  mode(countData) <- "integer"

  #design <- as.formula("~ condition", env = .GlobalEnv)
  object <- list(countData = countData, design = x)
  object
}
