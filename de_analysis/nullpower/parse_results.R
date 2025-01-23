
rm(list = ls())
require(tidyverse)
require(patchwork)

method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

# FPR test ####
beta <- 0.5
res <- dplyr::tibble()

for (is.pb in c(TRUE, FALSE)) {
  if (is.pb) {
    head_foler_null = "nullpower/null_subject"
    head_foler_pow = "nullpower/pow_subject"
  } else {
    head_foler_null = "nullpower/null_cell"
    head_foler_pow = "nullpower/pow_cell"
  }

  ll <- list.files(head_foler_null, full.names = T)
  #ll <- ll[grepl("bca", ll)]
  if (length(ll) > 0) {
    for (i in 1:length(ll)) {
      l <- ll[i]
      print(l)
      l <- unlist(str_split(l, "/"))[[3]]
      author <- unlist(str_split(l, ".n."))[[1]]
      n.patients <- as.numeric(unlist(str_split(l, ".n."))[[2]])
      n.genes <- as.numeric(unlist(strsplit(unlist(str_split(l, ".ngenes."))[[2]], ".ct."))[[1]])
      ct.index <- as.numeric(unlist(strsplit(unlist(strsplit(l, ".ct."))[2], ".prob"))[[1]])
      prob.de <- as.numeric(unlist(strsplit(unlist(strsplit(l, ".probde."))[2], ".iter"))[[1]])
      i.iter <- as.numeric(unlist(strsplit(unlist(strsplit(l, ".iter."))[2], ".csv"))[[1]])
      n_genes_de <- n.genes
      n_genes_non_de <- n.genes / prob.de - n.genes

      dnull <- read.delim(ll[i], sep=",") %>% dplyr::mutate(DE = FALSE, iter = i.iter)
      dpow <- read.delim(paste0(head_foler_pow,"/", l), sep=",") %>% dplyr::mutate(DE = TRUE, iter = i.iter)

      d <- dplyr::bind_rows(dpow, dnull)
      colnames(d) <- c("Gene", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)", "DE", "Iter")

      c <- colnames(d)[2]
      n_algo <- ncol(d) - 2
      
      r <- d %>% 
        tidyr::pivot_longer(!c(Gene, DE, Iter)) %>% 
        dplyr::select(DE, name, value) %>% 
        dplyr::group_by(name) %>% 
        dplyr::mutate(value = ifelse(is.na(value), 1, value)) %>% 
        dplyr::mutate(padj = p.adjust(value, "BH")) %>% 
        mutate(
          predicted = padj <= 0.05,
          TP = as.numeric(DE & predicted),    # True Positive
          TN = as.numeric(!DE & !predicted),  # True Negative
          FP = as.numeric(!DE & predicted),   # False Positive
          FN = as.numeric(DE & !predicted)    # False Negative
        ) %>%
        summarise(
          TP = sum(TP),
          TN = sum(TN),
          FP = sum(FP),
          FN = sum(FN),
          numerator = (TP * TN) - (FP * FN),
          denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
          MCC = ifelse(denominator == 0, 0, numerator / denominator)
        ) %>% 
        dplyr::select(name, MCC)

      # c <- colnames(d)[4]
      # r <- lapply(colnames(d)[2:n_algo], function(c) {
      #   rr <- dplyr::tibble()
      #   for (i in unique(d$Iter)) {
      #     pvals <- d[,c][d$Iter==i]
      #     pvals_null <- pvals[(n_genes_de+1):(n.genes/prob.de)]
      #     KS_test <- ks.test(pvals_null, punif)$statistic
      # 
      #     pvals[is.na(pvals)] <- runif(1,0,1)
      #     pvals_adj <- p.adjust(pvals, method = "BH")
      # 
      #     pred = pvals_adj <= .05
      #     #pred = pvals <= .05
      #     gt <- d$DE[d$Iter==i]
      # 
      #     PP = sum(pred == TRUE)
      #     PN = sum(pred == FALSE)
      #     TP = sum((pred == TRUE) * (gt == TRUE))
      #     FP = sum((pred == TRUE) * (gt == FALSE))
      #     TN = sum((pred == FALSE) * (gt == FALSE))
      #     FN = sum((pred == FALSE) * (gt == TRUE))
      #     P = sum(gt == TRUE)
      #     N = sum(gt == FALSE)
      # 
      #     FPR = FP / N
      #     TPR = TP / P
      #     FNR = FN / P
      #     TNR = TN / N
      #     PPV = TP / PP
      #     NPV = TN / PN
      #     FOR = FN / PN
      #     FDR = 1 - PPV
      #     F1 = (2 * TP) / (2 * TP + FP + FN)
      #     Fbeta = ((1 + beta^2) * TP) / ((1 + beta^2) * TP + beta^2 * FN + FP)
      # 
      #     if (is.nan(PPV)) {PPV <- 0}
      #     if (is.nan(FDR)) {FDR <- 0}
      # 
      #     MCC = sqrt((TPR * TNR * PPV * NPV)) - sqrt(FNR * FPR * FOR * FDR)
      #     ACC = (TP + TN) / (P + N)
      #     rr <- dplyr::bind_rows(rr, dplyr::tibble(
      #       name = c,
      #       FPR=FPR,
      #       TPR=TPR,
      #       ACC=ACC, FNR=FNR, TNR=TNR, PPV=PPV, NPV=NPV,
      #       FOR=FOR, FDR=FDR, MCC=MCC, F1=F1, Fbeta=Fbeta, KS_test=KS_test,
      #       iter=i))
      #   }
      #   rr
      # }) %>% do.call("bind_rows", .)

      r <- r %>% dplyr::mutate(ct.index = ct.index, is.pb=is.pb, author=author, patients=n.patients, ngenes=n.genes, i.iter=i.iter)

      res <- dplyr::bind_rows(res, r)

    }
  }
}

saveRDS(res, "nullpower/final_res/results.rds")
