
devil.base <- function(count, df){
  design_matrix <- model.matrix(~1+tx_cell, data = df)

  s <- Sys.time()
  #fit <- devil::fit_devil(count, design_matrix, size_factors=T, verbose=F, parallel.cores=1, min_cells=-1, avg_counts=-1)
  fit <- devil::fit_devil(count, design_matrix, size_factors=FALSE, verbose=F, parallel.cores=1, min_cells=-1, avg_counts=-1)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- devil::test_de(fit, contrast=as.array(c(0,1)))

  beta <- fit$beta[2,]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}

devil.mixed <- function(count, df) {
  design_matrix <- model.matrix(~1+tx_cell, data = df)
  clusters = as.factor(df$id)

  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, size_factors=F, verbose=T, min_cells=-1, avg_counts=-1, parallel.cores=4)
 # fit <- devil::fit_devil(count, design_matrix, size_factors=T, verbose=T, min_cells=-1, avg_counts=-1, parallel.cores=4)
  e <- Sys.time()
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters)

  beta <- fit$beta[2,]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
