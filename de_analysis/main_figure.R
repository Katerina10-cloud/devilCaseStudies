setwd("~/GitHub/de_analysis/")
rm(list = ls())
require(tidyverse)
require(patchwork)
source("plots.R")

author <- "hsc"
method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

p1 <- all_null_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
  ggtitle("Cell-wise") + labs(x = "Uniform quantile")
p2 <- all_null_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
  ggtitle("Patient-wise") + labs(x = "Uniform quantile")
p3 <- all_pow_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
  ggtitle("") + labs(x = "Uniform quantile")
p4 <- all_pow_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
  ggtitle("") + labs(x = "Uniform quantile")

pB <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme_bw() &
  theme(
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    #title = element_blank(),
    text = element_text(size = 10)
  )
pB
# null_plot_pb <- wrap_plots(all_null_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(size = 12))
# pow_plot_pb <- wrap_plots(all_pow_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(size = 12))
#
# null_plot_not_pb <- wrap_plots(all_null_plots(author, FALSE, algos = method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(size = 12))
# pow_plot_not_pb <- wrap_plots(all_pow_plots(author, FALSE, algos = method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(size = 12))
#
# p <- wrap_plots(null_plot_not_pb, null_plot_pb, pow_plot_not_pb, pow_plot_pb) &
#   theme_bw() &
#   theme(
#     legend.position = 'none', plot.title = element_blank(),
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     #title = element_blank(),
#     text = element_text(size = 10)
#   )
# p
ggsave("main_fig/example_pow_and_null.pdf", dpi=300, width = 105, height = 90, units = 'mm')

# FPR test ####
#cell_indexes <- c(1:6)
#authors <- c("hsc")
#ngenes <- c(500, 1000, 2000)
#ngenes <- c(100, 500, 1000)

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
  if (length(ll) > 0) {
    for (i in 1:length(ll)) {
      l <- ll[i]
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

      c <- colnames(d)[4]
      r <- lapply(colnames(d)[2:n_algo], function(c) {
        rr <- dplyr::tibble()
        for (i in unique(d$Iter)) {
          pvals <- d[,c][d$Iter==i]
          pvals_null <- pvals[(n_genes_de+1):(n.genes/prob.de)]
          KS_test <- ks.test(pvals_null, punif)$statistic

          pvals[is.na(pvals)] <- runif(1,0,1)
          pvals_adj <- p.adjust(pvals, method = "BH")

          pred = pvals_adj <= .05
          #pred = pvals <= .05
          gt <- d$DE[d$Iter==i]

          PP = sum(pred == TRUE)
          PN = sum(pred == FALSE)
          TP = sum((pred == TRUE) * (gt == TRUE))
          FP = sum((pred == TRUE) * (gt == FALSE))
          TN = sum((pred == FALSE) * (gt == FALSE))
          FN = sum((pred == FALSE) * (gt == TRUE))
          P = sum(gt == TRUE)
          N = sum(gt == FALSE)

          FPR = FP / N
          TPR = TP / P
          FNR = FN / P
          TNR = TN / N
          PPV = TP / PP
          NPV = TN / PN
          FOR = FN / PN
          FDR = 1 - PPV
          F1 = (2 * TP) / (2 * TP + FP + FN)
          Fbeta = ((1 + beta^2) * TP) / ((1 + beta^2) * TP + beta^2 * FN + FP)

          if (is.nan(PPV)) {PPV <- 0}
          if (is.nan(FDR)) {FDR <- 0}

          MCC = sqrt((TPR * TNR * PPV * NPV)) - sqrt(FNR * FPR * FOR * FDR)
          ACC = (TP + TN) / (P + N)
          rr <- dplyr::bind_rows(rr, dplyr::tibble(
            name = c,
            FPR=FPR,
            TPR=TPR,
            ACC=ACC, FNR=FNR, TNR=TNR, PPV=PPV, NPV=NPV,
            FOR=FOR, FDR=FDR, MCC=MCC, F1=F1, Fbeta=Fbeta, KS_test=KS_test,
            iter=i))
        }
        rr
      }) %>% do.call("bind_rows", .)

      r <- r %>% dplyr::mutate(ct.index = ct.index, is.pb=is.pb, author=author, patients=n.patients, ngenes=n.genes)

      res <- dplyr::bind_rows(res, r)

    }
  }
}

a <- "bca"

pfalse <- res %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
  #dplyr::filter(name %in% c("glmGamPoi (cell)")) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  ggplot(mapping = aes(x=name, y=MCC, col=name, fill=name)) +
  geom_violin(col='black') +
  geom_boxplot(col='black', width=.1, fill="white") +
  #geom_jitter() +
  ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  #scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  theme_bw() +
  labs(x = "", y = "MCC", fill="Algorithm") +
  theme(legend.position = 'none', text = element_text(size = 12), axis.text.x = element_blank()) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3))

ptrue <- res %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  ggplot(mapping = aes(x=name, y=MCC, fill=name)) +
  geom_violin(col="black") +
  geom_boxplot(col='black', width=.1, fill="white") +
  ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = method_colors) +
  labs(x = "", y = "MCC", fill="Algorithm") +
  theme(legend.position = 'none') +
  theme(legend.position = 'none', text = element_text(size = 12), axis.text.x = element_blank()) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=3))

pfalse
ptrue
#ggsave(filename = paste0("plot_figure/true_comparison.svg"), dpi=300, plot = ptrue, width = 145, height = 125, units='mm')



# Table
table_true <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(Ngenes, Npatients, name, Dataset) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  dplyr::summarise(medianMCC = stats::quantile(MCC, .5)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Ngenes, Npatients, Dataset) %>%
  dplyr::mutate(is_max = medianMCC == max(medianMCC)) %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=name, y=1, fill=medianMCC, label=round(medianMCC, 2))) +
  geom_tile() +
  geom_text(color='white') +
  scale_fill_gradient(low="#a6bddb", high = "#0570b0") +
  labs(x = "Algorithm", y = "Median MCC") +
  ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size = 12), axis.text.y = element_blank())

table_false <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(Ngenes, Npatients, name, Dataset) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  dplyr::summarise(medianMCC = stats::quantile(MCC, .5)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Ngenes, Npatients, Dataset) %>%
  dplyr::mutate(is_max = medianMCC == max(medianMCC)) %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=name, y=1, fill=medianMCC, label=round(medianMCC, 2))) +
  geom_tile() +
  geom_text(color='white') +
  scale_fill_gradient(low="lightblue", high = "steelblue4") +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  labs(x = "", y = "Median MCC") +
  ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size = 12), axis.text.y = element_blank())


res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
  # dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  # dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(ngenes, patients, name, author) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  dplyr::summarise(medianMCC = stats::quantile(MCC, .5)) %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=ngenes, y=medianMCC, linetype = as.factor(patients), col=name)) +
  geom_point(position=position_dodge(width=5)) +
  geom_line(position=position_dodge(width=5)) +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  labs(x = "N genes", y = "Median MCC", col="Algorithm", linetype = "N patients") +
  #ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12))

res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
  # dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  # dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(ngenes, patients, name, author) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=ngenes, y=medianMCC, ymax = highMCC, ymin=lowMCC,linetype = as.factor(patients), col=name, shape=as.factor(patients))) +
  geom_pointrange(position=position_dodge(width=5)) +
  geom_line(position=position_dodge(width=5)) +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  ggtitle(paste0(a, " dataset")) +
  labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
  #ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12))


ptrue <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
  # dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  # dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(ngenes, patients, name, author) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  #dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, linetype = as.factor(patients), col=name, shape=as.factor(patients))) +
  # geom_pointrange(position=position_dodge(width=5)) +
  # geom_line(position=position_dodge(width=5)) +
  geom_boxplot() +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  ggtitle(paste0("Patient-wise")) +
  labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
  #ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = 'bottom')

pfalse <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
  # dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  # dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(ngenes, patients, name, author) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  #dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, linetype = as.factor(patients), col=name, shape=as.factor(patients))) +
  # geom_pointrange(position=position_dodge(width=5)) +
  # geom_line(position=position_dodge(width=5)) +
  geom_boxplot() +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  ggtitle(paste0("Cell-wise")) +
  labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
  #ggh4x::facet_nested(Npatients+Ngenes~Dataset, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12))

pfalse <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(ngenes, patients, name, author, Npatients) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  #dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
  # geom_pointrange(position=position_dodge(width=5)) +
  # geom_line(position=position_dodge(width=5)) +
  geom_boxplot() +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  ggtitle(paste0("Cell-wise")) +
  labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
  #facet_wrap(~patients, ncol = 2) +
  ggh4x::facet_nested(~Npatients, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12))

ptrue <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(ngenes, patients, name, author, Npatients) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  #dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
  # geom_pointrange(position=position_dodge(width=5)) +
  # geom_line(position=position_dodge(width=5)) +
  geom_boxplot() +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  ggtitle(paste0("Patient-wise")) +
  labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
  #facet_wrap(~patients, ncol = 2) +
  ggh4x::facet_nested(~Npatients, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12))


ptrue
pfalse


ggsave(filename = "main_fig/p_legend.svg", dpi=300, width = 200, height = 100, plot = ptrue, units = 'mm')
ggsave(filename = "main_fig/p_legend.pdf", dpi=300, width = 200, height = 100, plot = ptrue, units = 'mm')
ptrue <- ptrue + theme(legend.position = 'none')
pfalse <- pfalse + theme(legend.position = 'none')
ggsave(filename = "main_fig/pfalse.svg", dpi=300, width = 110, height = 100, plot = pfalse, units = 'mm')
ggsave(filename = "main_fig/ptrue.svg", dpi=300, width = 110, height = 100, plot = ptrue, units = 'mm')
ggsave(filename = "main_fig/pfalse.pdf", dpi=300, width = 200, height = 100, plot = pfalse, units = 'mm')
ggsave(filename = "main_fig/ptrue.pdf", dpi=300, width = 200, height = 100, plot = ptrue, units = 'mm')
# ggsave(filename = "main_fig/table_true.pdf", dpi=300, width = 35, height = 110, plot = table_true, units = 'mm')
# ggsave(filename = "main_fig/table_false.pdf", dpi=300, width = 35, height = 110, plot = table_false, units = 'mm')
#(null_plot_not_pb | pow_plot_not_pb | pfalse | table_false) / (null_plot_pb | pow_plot_pb | ptrue | table_true) & theme(legend.position = "none", text = element_text(size = 12), title = element_blank())


timing_res <- readRDS(paste0("nullpower/timing_results/", a,".rds"))

timing_plot <- timing_res %>%
  dplyr::filter(algo %in% c("limma (Pb)", "Devil (base)", "glmGamPoi (cell)", "Nebula")) %>%
  dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
  dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>%
  ggplot(mapping = aes(x=cell_order, y=timings, col=algo)) +
  geom_boxplot() +
  scale_color_manual(values = method_colors) +
  labs(x = "N cells", y="Time (s)", col = "Algorithm") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ggtitle("")

ggsave(filename = "main_fig/timing.svg", dpi=300, width = 80, height = 100, plot = last_plot(), units = 'mm')
ggsave(filename = "main_fig/timing.pdf", dpi=300, width = 80, height = 100, plot = last_plot(), units = 'mm')

layout <- "
###AA
###AA
EEFFG
EEFFG
"

pB <- p1 + p2 + p3 + p4 + plot_layout(axis_titles = "collect") &
  theme_bw() &
  theme(
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    #title = element_blank(),
    text = element_text(size = 10)
  )

final_plot <- patchwork::wrap_plots(A = pB, E = ptrue, F = pfalse, G = timing_plot, design = layout) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none", text = element_text(size = 12))
ggsave(filename = "~/Dropbox/2023. DEVIL/figures/main2_v0.svg", plot = final_plot, dpi = 300, width = 13, height = 9, units = "in")
