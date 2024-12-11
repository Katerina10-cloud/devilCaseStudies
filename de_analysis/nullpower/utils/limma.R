
limma.mult <- function(count, df){
  sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
  sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(id, tx_cell),
    verbose=FALSE
  )

  design <- model.matrix(~1+tx_cell, data=colData(sce.pb))
  s <- Sys.time()
  edger.obj <- edgeR::DGEList(counts(sce.pb))
  v <- limma::voom(edger.obj, design)
  vfit <- limma::lmFit(v, design)
  efit <- limma::eBayes(vfit)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()

  beta <- efit$coefficients[,2] * log(2)
  pval <- efit$p.value[,2]
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
