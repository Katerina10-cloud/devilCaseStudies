require(glmGamPoi)
require(SingleCellExperiment)

glmgp.mult <- function(count, df){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(id, tx_cell),
    verbose=FALSE
  )

  s <- Sys.time()
  fit <- glmGamPoi::glm_gp(sce.pb, design=~1+tx_cell)
  #fit <- glmGamPoi::glm_gp(sce.pb, design=~1+tx_cell, size_factors='normed_sum')
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- glmGamPoi::test_de(fit, reduced_design=~1)

  beta <- fit$Beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

glmgp.cell.mult <- function(count, df){
  design_matrix <- model.matrix(~1+tx_cell, data = df)

  s <- Sys.time()
  fit <- glmGamPoi::glm_gp(count, design_matrix, on_disk=FALSE, size_factors=FALSE)
  #fit <- glmGamPoi::glm_gp(count, design_matrix, on_disk=FALSE, size_factors='normed_sum')
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- glmGamPoi::test_de(fit, reduced_design=~1)

  beta <- fit$Beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
