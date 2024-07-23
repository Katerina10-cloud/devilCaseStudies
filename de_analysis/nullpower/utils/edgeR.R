
require(glmGamPoi)
require(SingleCellExperiment)
require(edgeR)

edger.mult <- function(count, df){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
  sce.pb <- glmGamPoi::pseudobulk_sce(
    sce.obj,
    group_by=vars(id, tx_cell),
    verbose=FALSE
  )

  design <- model.matrix(~1+tx_cell, data=colData(sce.pb))
  edger.obj <- edgeR::DGEList(counts(sce.pb))
  edger.obj <- edgeR::estimateDisp(edger.obj, design)
  s <- Sys.time()
  fit <- edgeR::glmQLFit(y=edger.obj, design=design)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- edgeR::glmTreat(fit, coef=2)

  beta <- test$coefficients[,2]
  pval <- test$table[,'PValue']
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
