library(devil)

devil.base <- function(count, df){
  design_matrix <- model.matrix(~1+tx_cell, data = df)

  s <- Sys.time()
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=TRUE, verbose=F, 
                          parallel.cores=1, init_overdispersion = 100, 
                          offset = 1e-6, max_iter = 200, tolerance = 1e-3)
  e <- Sys.time()
  fit$input_parameters$parallel = FALSE
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- devil::test_de(fit, contrast=as.array(c(0,1)))

  beta <- fit$beta[,2]
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
  fit <- devil::fit_devil(count, design_matrix, 
                          size_factors=TRUE, verbose=F, parallel.cores=1, 
                          init_overdispersion = 100, offset = 1e-6, 
                          max_iter = 200, tolerance = 1e-3)
  e <- Sys.time()
  fit$input_parameters$parallel = FALSE
  delta_time <- difftime(e, s, units = "secs") %>% as.numeric()
  test <- devil::test_de(fit, contrast=as.array(c(0,1)), clusters=clusters)

  beta <- fit$beta[,2]
  pval <- test$pval
  tval <- qnorm(1-pval/2) * sign(beta)
  se <- beta/tval

  result <- dplyr::as_tibble(cbind(beta, se, tval, pval))
  result$delta_time <- delta_time
  colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)', 'Time')
  return(result %>% as.matrix())
}
