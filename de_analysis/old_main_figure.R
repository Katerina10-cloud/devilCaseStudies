
rm(list = ls())
require(tidyverse)
require(patchwork)
source("plots.R")

AUTHOR <- author <- "yazar"
method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

d <- dplyr::bind_rows(
  all_null_plots(author, FALSE, algos =method_cellwise, only_tibble=TRUE) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(is_de = FALSE, method = "Cell-wise"),
  all_pow_plots(author, FALSE, algos =method_cellwise, only_tibble = TRUE) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(is_de = TRUE, method = "Cell-wise")
)

threshold <- 0.05
d %>% 
  dplyr::group_by(name, ct.index, n.genes, n.samples, method, i.iter) %>% 
  dplyr::mutate(padj = p.adjust(observed_p_value, "BH")) %>% 
  mutate(
    predicted = padj <= threshold,
    TP = as.numeric(is_de & predicted),    # True Positive
    TN = as.numeric(!is_de & !predicted),  # True Negative
    FP = as.numeric(!is_de & predicted),   # False Positive
    FN = as.numeric(is_de & !predicted)    # False Negative
  ) %>%
  summarise(
    TP = sum(TP),
    TN = sum(TN),
    FP = sum(FP),
    FN = sum(FN),
    numerator = (TP * TN) - (FP * FN),
    denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    mcc = ifelse(denominator == 0, 0, numerator / denominator)
  ) %>% 
  ggplot(mapping = aes(x=as.factor(n.genes), y=mcc, col=name, fill=name)) +
  geom_boxplot() +
  facet_grid(~n.samples)

dplyr::bind_rows(d1[[1]], d1[[2]]) %>% 
  dplyr::filter(observed_p_value <= .05) %>% 
  dplyr::group_by(name, ct.index, n.samples) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(mapping = aes(x=name, y=n)) +
  geom_col() +
  facet_grid(n.samples ~ ct.index)


dplyr::bind_rows(d3[[1]], d3[[2]]) %>% 
  dplyr::filter(observed_p_value <= .05) %>% 
  dplyr::group_by(name, ct.index, n.samples) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(mapping = aes(x=name, y=n)) +
  geom_col() +
  facet_grid(n.samples ~ ct.index)

# pB <- dplyr::bind_rows(d1, d2, d3, d4) %>%
#   dplyr::group_by(xtype, ytype) %>%
#   dplyr::sample_n(500) %>%
#   ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
#   geom_point(position = position_dodge(width = .02)) +
#   geom_line(linewidth = .8, position = position_dodge(width = .02)) +
#   #geom_line(linewidth = 1, position = position_dodge(width = 0.01)) +
#   scale_color_manual(values = method_colors) +
#   facet_grid(ytype~xtype, scales = "free") +
#   theme_bw() +
#   labs(x = "Uniform quantiles", y="", col="Algorithm") +
#   #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
#   scale_color_manual(values = method_colors) +
#   #facet_wrap(~paste0(n.genes, " genes")) +
#   theme(legend.position = "bottom") +
#   scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n=3))
# pB


method_levels <- c("limma", "glmGamPoi", "glmGamPoi (cell)", "Nebula", "NEBULA", "Devil (mixed)", "Devil (base)", "Devil", "devil")
#dplyr::bind_rows(d1, d2, d3, d4) %>%
pB <- dplyr::bind_rows(d1, d2) %>%
  dplyr::mutate(name = dplyr::if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = dplyr::if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = dplyr::if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::group_by(xtype, ytype) %>%
  dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
  ggridges::geom_density_ridges(alpha = .7, scale = 1) +
  scale_color_manual(values = sort(method_colors)) +
  scale_fill_manual(values = sort(method_colors)) +
  facet_grid(~xtype, scales = "free") +
  theme_bw() +
  labs(x = "p value", y="", col="Algorithm") +
  scale_color_manual(values = method_colors) +
  #facet_wrap(~paste0(n.genes, " genes")) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
  theme(legend.position = "none") +
  theme()
pB

pC <- dplyr::bind_rows(d3, d4) %>%
  dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::group_by(xtype, ytype) %>%
  dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
  ggridges::geom_density_ridges(alpha = .7, scale = 1) +
  scale_color_manual(values = sort(method_colors)) +
  scale_fill_manual(values = sort(method_colors)) +
  facet_grid(~xtype, scales = "free") +
  theme_bw() +
  labs(x = "-log10 p value", y="", col="Algorithm") +
  scale_color_manual(values = method_colors) +
  #facet_wrap(~paste0(n.genes, " genes")) +
  theme(legend.position = "bottom") +
  scale_x_continuous(transform = "log10") +
  #scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
  theme(legend.position = "none") +
  theme()
pB
pC
# p1 <- all_null_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
#   ggtitle("Cell-wise") + labs(x = "Uniform quantile")
# p2 <- all_null_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
#   ggtitle("Patient-wise") + labs(x = "Uniform quantile")
# p3 <- all_pow_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
#   ggtitle("") + labs(x = "Uniform quantile")
# p4 <- all_pow_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20))[[1]] +
#   ggtitle("") + labs(x = "Uniform quantile")
#
# pB <- (p1 | p2) / (p3 | p4) +
#   plot_layout(guides = "collect") &
#   theme_bw() &
#   theme(
#     legend.position = 'none',
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     #title = element_blank(),
#     text = element_text(size = 10)
#   )
# pB
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
#ggsave("main_fig/example_pow_and_null.pdf", dpi=300, width = 105, height = 90, units = 'mm')

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

a <- AUTHOR

pfalse <- res %>%
  na.omit() %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
  dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::filter(author == a) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(Npatients, ngenes) %>%
  #dplyr::mutate(mean_MCC_devil = median(MCC[name == 'Devil (base)'])) %>%
  # dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  #dplyr::group_by(ngenes, patients, name, author, Npatients) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  #dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
  #ggplot(mapping = aes(x=name, y=MCC, col=name, yintercept=mean_MCC_devil)) +
  # geom_pointrange(position=position_dodge(width=5)) +
  # geom_line(position=position_dodge(width=5)) +
  geom_boxplot() +
  #geom_hline(mapping = aes(yintercept = mean_MCC_devil), linetype='dashed', color="darkslategray") +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  #ggtitle(paste0("Cell-wise")) +
  labs(x = "N genes", y = "MCC", col="", linetype = "N patients", shape = "N patients") +
  #facet_wrap(~patients, ncol = 2) +
  #ggh4x::facet_nested(~Npatients, scales = "free_y") +
  ggh4x::facet_nested(~Npatients, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "top")
pfalse

ptrue <- res %>%
  na.omit() %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
  dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::filter(author == a) %>%
  dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = paste0(patients, " patients"), Ngenes = paste0(ngenes, " genes")) %>%
  dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
  dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
  dplyr::group_by(Npatients, ngenes) %>%
  dplyr::mutate(mean_MCC_devil = median(MCC[name == 'Devil (mixed)'])) %>%
  #dplyr::summarise(medianMCC = median(MCC)) %>%
  #dplyr::summarise(medianMCC = stats::quantile(MCC, .5), lowMCC = stats::quantile(MCC, .25), highMCC = stats::quantile(MCC, .75)) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(is_max = ifelse(!is_max, "gray", name)) %>%
  ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
  #ggplot(mapping = aes(x=name, y=MCC, col=name)) +
  # geom_pointrange(position=position_dodge(width=5)) +
  # geom_line(position=position_dodge(width=5)) +
  geom_boxplot() +
  #geom_hline(mapping = aes(yintercept = mean_MCC_devil), linetype='dashed', color="darkslategray") +
  #geom_pointrange() +
  #geom_line() +
  scale_color_manual(values = method_colors) +
  #scale_fill_manual(values = c("sienna3", "plum4")) +
  #ggtitle(paste0("Patient-wise")) +
  labs(x = "N genes", y = "MCC", col="", linetype = "N patients", shape = "N patients") +
  #facet_wrap(~patients, ncol = 2) +
  ggh4x::facet_nested(~Npatients, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "top")

ptrue
pfalse

# kolmogorov smirnof plots - FALSE ####
dx <- .15

r_ks <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == FALSE, name %in% method_cellwise)

r_ks_tot <- lapply(unique(r_ks$name), function(n) {
  MCCs <- r_ks %>%
    dplyr::filter(name==n) %>%
    dplyr::pull(MCC)
  ECDFs <- lapply(MCCs, function(mcc) {
    sum(MCCs <= mcc) / length(MCCs)
  }) %>% unlist()

  dplyr::tibble(MCCs =MCCs, ECDFs=ECDFs, name=n)
}) %>% do.call("rbind", .)

ks_pvals <- lapply(method_cellwise[method_cellwise!="Devil (base)"], function(m) {
  print(m)
  pval = ks.test(filter(r_ks_tot, name=="Devil (base)")$MCCs, filter(r_ks_tot, name==m)$MCCs)$p.value
  if (pval == 0) {
    pval <- "< 2e-16"
  } else {
    pval <- paste0("= ", round(pval, 2))
  }
  dplyr::tibble(m=m, pval=pval)
}) %>% do.call("bind_rows", .)

ks_false <- r_ks_tot %>%
  dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
  geom_line(linewidth = .8, position = position_dodge(width = .02)) +
  scale_color_manual(values = sort(method_colors)) +
  #ggtitle(paste0("Cell-wise")) +
  labs(x = "MCC", y = "Probability", col="") +
  theme_bw() +
  theme(text = element_text(size = 12))
ks_false

L = ggplot_build(ks_false)$layout$panel_params[[1]]
Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

ks_false <- ks_false +
  annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
  annotate(geom='label', x=Lx, y=Ly - dx, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
  annotate(geom='label', x=Lx, y=Ly - 2*dx, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
  theme(legend.position = "left")

ks_false
# kolmogorov smirnof plots - TRUE ####
r_ks <- res %>%
  na.omit() %>%
  dplyr::filter(author == a) %>%
  dplyr::filter(is.pb == TRUE, name %in% method_patientwise)

r_ks_tot <- lapply(unique(r_ks$name), function(n) {
  MCCs <- r_ks %>%
    dplyr::filter(name==n) %>%
    dplyr::pull(MCC)
  ECDFs <- lapply(MCCs, function(mcc) {
    sum(MCCs <= mcc) / length(MCCs)
  }) %>% unlist()

  dplyr::tibble(MCCs =MCCs, ECDFs=ECDFs, name=n)
}) %>% do.call("rbind", .)

ks_pvals <- lapply(method_patientwise[method_patientwise!="Devil (mixed)"], function(m) {
  print(m)
  pval = ks.test(filter(r_ks_tot, name=="Devil (mixed)")$MCCs, filter(r_ks_tot, name==m)$MCCs)$p.value
  true_pval <- pval
  if (pval <= 1e-6) {
    pval <- "< 1e-6"
  } else {
    pval <- paste0("= ", round(pval, 3))
  }
  dplyr::tibble(m=m, pval=pval, true_pval=true_pval)
}) %>% do.call("bind_rows", .)

ks_true <- r_ks_tot %>%
  dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::mutate(name = factor(name, levels = method_levels)) %>%
  ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
  geom_line(linewidth = .8, position = position_dodge(width = .02)) +
  scale_color_manual(values = method_colors) +
  #ggtitle(paste0("Patient-wise")) +
  labs(x = "MCC", y = "Probability", col="") +
  theme_bw() +
  theme(text = element_text(size = 12))

L = ggplot_build(ks_true)$layout$panel_params[[1]]
Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

ks_true <- ks_true +
  annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
  annotate(geom='label', x=Lx, y=Ly - dx, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
  annotate(geom='label', x=Lx, y=Ly - 2*dx, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
  theme(legend.position = "left")


ks_true
ks_false


#ggsave(filename = "main_fig/p_legend.svg", dpi=300, width = 200, height = 100, plot = ptrue, units = 'mm')
#ggsave(filename = "main_fig/p_legend.pdf", dpi=300, width = 200, height = 100, plot = ptrue, units = 'mm')
#ggsave(filename = "main_fig/pfalse.svg", dpi=300, width = 110, height = 100, plot = pfalse, units = 'mm')
#ggsave(filename = "main_fig/ptrue.svg", dpi=300, width = 110, height = 100, plot = ptrue, units = 'mm')
#ggsave(filename = "main_fig/pfalse.pdf", dpi=300, width = 200, height = 100, plot = pfalse, units = 'mm')
#ggsave(filename = "main_fig/ptrue.pdf", dpi=300, width = 200, height = 100, plot = ptrue, units = 'mm')
# ggsave(filename = "main_fig/table_true.pdf", dpi=300, width = 35, height = 110, plot = table_true, units = 'mm')
# ggsave(filename = "main_fig/table_false.pdf", dpi=300, width = 35, height = 110, plot = table_false, units = 'mm')
#(null_plot_not_pb | pow_plot_not_pb | pfalse | table_false) / (null_plot_pb | pow_plot_pb | ptrue | table_true) & theme(legend.position = "none", text = element_text(size = 12), title = element_blank())


timing_res <- readRDS(paste0("nullpower/timing_results/", a,".rds"))

timing_plot <- timing_res %>%
  dplyr::filter(algo %in% c("Devil (base)", "glmGamPoi (cell)", "Nebula")) %>%
  dplyr::rename(name = algo) %>%
  dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
  dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
  dplyr::rename(algo = name) %>%
  dplyr::group_by(author, is.pb, n.sample, n.gene, int.ct, iter) %>%
  dplyr::mutate(time_fold = timings/timings[algo == "devil"]) %>%
  dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
  dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>%
  dplyr::filter(algo %in% c("glmGamPoi", "NEBULA")) %>%
  #ggplot(mapping = aes(x=cell_order, y=timings, col=algo)) +
  ggplot(mapping = aes(x=cell_order, y=time_fold, col=algo)) +
  geom_boxplot() +
  scale_color_manual(values = method_colors) +
  labs(x = "N cells", y="Time fold", col = "") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "left") +
  ggtitle("") +
  geom_hline(yintercept = 1, color = "darkslategray", linetype = 'dashed') +
  coord_flip()
timing_plot

#ggsave(filename = "main_fig/timing.svg", dpi=300, width = 80, height = 100, plot = last_plot(), units = 'mm')
#ggsave(filename = "main_fig/timing.pdf", dpi=300, width = 80, height = 100, plot = last_plot(), units = 'mm')

layout <- "
KKKKAAAA
KKKKAAAA
EEEEFFFF
EEEEFFFF
HHHHLLLL
HHHHLLLL
GGGGGGGG
"

layout <- "
AABB
AABB
AABB
AACC
AACC
AACC
DDFF
DDFF
DDFF
DDFF
DDGG
DDGG
EEGG
EEGG
EEHH
EEHH
EEHH
EEHH
"

final_plot <- patchwork::wrap_plots(
  A = ggplot(),
  B = free(pB),
  C = free(pC),
  D = free(pfalse),
  E = free(ptrue),
  F = free(ks_false),
  G=free(ks_true),
  H=free(timing_plot),
  design = layout) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 12))
final_plot
ggsave(filename = "main_fig/main_v0.pdf", plot = final_plot, dpi = 400, width = 11.7, height = 11.7, units = "in")
ggsave(filename = "~/Dropbox/2023. DEVIL/figures/main2_v0.svg", plot = final_plot, dpi = 400, width = 11.7, height = 11.7, units = "in")
