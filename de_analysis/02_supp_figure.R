
rm(list = ls())
require(tidyverse)
require(ggh4x)
require(patchwork)
# source("utils.R")
source("utils_plots.R")

author_name <- "yazar"
for (author_name in c("bca", "yazar", "kumar", "hsc")) {
  author <- AUTHOR <- author_name
  method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
  method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

  is.pb <- F

  method_wise <- ifelse(is.pb==T, "Patient-wise", "Cell-wise")

  null_plots_cellwise <- all_null_plots(author, FALSE, algos = method_cellwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "p-value", xtype="Cell-wise")
  pow_plots_cellwise <- all_pow_plots(author, FALSE, algos = method_cellwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "-log10 p-value", xtype="Cell-wise")

  pA <- null_plots_cellwise %>%
    #sample_n(1000) %>%
    ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
    geom_point(size = .5) +
    scale_color_manual(values = method_colors) +
    #facet_grid(ytype+n.genes~xtype, scales = "free") +
    ggh4x::facet_nested(paste(n.samples, " P")+paste(n.genes, " G")~xtype, scales = "free") +
    theme_bw() +
    labs(x = "Uniform quantiles", y="p-value", col="Algorithm") +
    #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c(0,1)) +
    scale_y_continuous(breaks = c(0,1))

  pB <- pow_plots_cellwise %>%
    #sample_n(size = 1000, replace = T) %>%
    ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
    geom_point(size = .5) +
    geom_line() +
    scale_color_manual(values = method_colors) +
    #facet_grid(ytype+n.genes~xtype, scales = "free") +
    ggh4x::facet_nested(paste(n.samples, " P")+paste(n.genes, " G")~xtype, scales = "free") +
    theme_bw() +
    labs(x = "Uniform quantiles", y="-log10 p-value", col="Algorithm") +
    #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c(0,1)) +
    scale_y_continuous(trans = "log10", breaks = c(.1,100))

  null_plots_pwise <- all_null_plots(author, TRUE, method_patientwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "p-value", xtype="Patient-wise")
  pow_plots_pwise <- all_pow_plots(author, TRUE, method_patientwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "-log10 p-value", xtype="Patient-wise")

  pC <- null_plots_pwise %>%
    sample_n(1000) %>%
    ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
    geom_point(size = .5) +
    scale_color_manual(values = method_colors) +
    #facet_grid(ytype+n.genes~xtype, scales = "free") +
    ggh4x::facet_nested(paste(n.samples, " P")+paste(n.genes, " G")~xtype, scales = "free") +
    theme_bw() +
    labs(x = "Uniform quantiles", y="p-value", col="Algorithm") +
    #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c(0,1)) +
    scale_y_continuous(breaks = c(0,1))

  pD <- pow_plots_pwise %>%
    #sample_n(size = 1000, replace = T) %>%
    ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
    geom_point(size = .5) +
    geom_line() +
    scale_color_manual(values = method_colors) +
    #facet_grid(ytype+n.genes~xtype, scales = "free") +
    ggh4x::facet_nested(paste(n.samples, " P")+paste(n.genes, " G")~xtype, scales = "free") +
    theme_bw() +
    labs(x = "Uniform quantiles", y="-log10 p-value", col="Algorithm") +
    #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c(0,1)) +
    scale_y_continuous(trans = "log10", breaks = c(.1,100))

  pC
  pD

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
    #ggtitle(paste0("Cell-wise")) +
    labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
    #facet_wrap(~patients, ncol = 2) +
    ggh4x::facet_nested(~Npatients, scales = "free_y") +
    theme_bw() +
    #theme(text = element_text(size = 12)) +
    theme(legend.position = "none")

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
    #ggtitle(paste0("Patient-wise")) +
    labs(x = "N genes", y = "MCC", col="Algorithm", linetype = "N patients", shape = "N patients") +
    #facet_wrap(~patients, ncol = 2) +
    ggh4x::facet_nested(~Npatients, scales = "free_y") +
    theme_bw() +
    #theme(text = element_text(size = 12)) +
    theme(legend.position = "none")

  pfalse
  ptrue

  # kolmogorov smirnof plots - FALSE ####
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
    ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
    geom_line(linewidth = .8) +
    scale_color_manual(values = method_colors) +
    #ggtitle(paste0("Cell-wise")) +
    labs(x = "MCC", y = "Probability", col="Algorithm") +
    theme_bw() +
    theme(text = element_text(size = 12))

  L = ggplot_build(ks_false)$layout$panel_params[[1]]
  Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
  Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

  ks_false <- ks_false +
    annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
    annotate(geom='label', x=Lx, y=Ly - .12, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
    annotate(geom='label', x=Lx, y=Ly - .24, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
    theme(legend.position = "none")

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
    ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
    geom_line(linewidth = .8) +
    scale_color_manual(values = method_colors) +
    #ggtitle(paste0("Patient-wise")) +
    labs(x = "MCC", y = "Probability", col="Algorithm") +
    theme_bw() +
    theme(text = element_text(size = 12))

  L = ggplot_build(ks_true)$layout$panel_params[[1]]
  Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
  Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

  ks_true <- ks_true +
    annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]], fill='white') +
    annotate(geom='label', x=Lx, y=Ly - .12, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
    annotate(geom='label', x=Lx, y=Ly - .24, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
    theme(legend.position = "none")


  ks_true
  ks_false

  # timing plot ####
  timing_res <- readRDS(paste0("nullpower/timing_results/", a,".rds"))
  timing_plot <- timing_res %>%
    dplyr::filter(algo %in% c("Devil (base)", "glmGamPoi (cell)", "Nebula")) %>%
    dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
    dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>%
    ggplot(mapping = aes(x=cell_order, y=timings, col=algo)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "N cells", y="Time (s)", col = "Algorithm") +
    theme_bw() +
    theme(text = element_text(size = 12)) +
    theme(legend.position = "none") +
    ggtitle("") +
    scale_y_continuous(transform = 'log10')
  timing_plot


  # final figure ####
  layout <- "
  ABCD
  ABCD
  ABCD
  ABCD
  ABCD
  EEFF
  EEFF
  EEFF
  EEHH
  EEHH
  GGHH
  GGHH
  GGII
  GGII
  GGII
  LLLL
  "

  pA <- pA + guides(color = guide_legend(override.aes = list(size=2)))
  leg <- ggpubr::get_legend(pA)
  leg <- ggpubr::as_ggplot(leg)

  final_plot <- patchwork::wrap_plots(
    A = free(pA),
    B = free(pB),
    C = free(pC),
    D = free(pD),
    E = free(pfalse),
    F = free(ks_false),
    G = free(ptrue),
    H = free(ks_true),
    I=free(timing_plot + coord_flip()),
    L = leg + plot_layout(tag_level = "new"),
    design = layout) +
    plot_annotation(tag_levels = "A") &
    #plot_layout(guides = "collect") &
    theme(legend.position = "none", text = element_text(size = 12))
  # final_plot

  ggsave(paste0("figures/", AUTHOR, ".png"), dpi=400, width = 10, height = 13, plot = final_plot, units = 'in')
}
