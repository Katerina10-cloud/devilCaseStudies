setwd("~/GitHub/devilCaseStudies/de_analysis")
rm(list = ls())
require(tidyverse)
require(patchwork)
# source("utils.R")
source("utils_plots.R")

author_name <- "bca"

method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

# P-values Cell-wise ####
plots <- lapply(c("bca", "yazar", "kumar", "hsc"), function(author_name) {
  print(author_name)
  author <- AUTHOR <- author_name


  is.pb <- F
  method_wise <- ifelse(is.pb==T, "Patient-wise", "Cell-wise")

  nulls <- all_null_plots(author, FALSE, algos = method_cellwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "p-value", xtype="Cell-wise")

  pows <- all_pow_plots(author, FALSE, algos = method_cellwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "p-value", xtype="Cell-wise")


  pows$observed_p_value[pows$observed_p_value == 0] = min(pows$observed_p_value[pows$observed_p_value != 0])

  p1 = nulls %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = "p-value", y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    facet_grid(paste0(n.samples, " samples")~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme() +
    ggtitle(paste0(author_name, " dataset"))

  p2 = pows %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = "-log10 p-value", y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    facet_grid(paste0(n.samples, " samples")~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(transform = "log10") +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme() +
    ggtitle(paste0(author_name, " dataset"))

  p1
  p2

  list(p1, p2)
})

plots <- unlist(plots, recursive = FALSE)

cellwise_pvalues_plot <- patchwork::wrap_plots(plots, ncol = 2, nrow = 4) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
cellwise_pvalues_plot

# P-values Patient-wise ####
plots <- lapply(c("bca", "yazar", "kumar", "hsc"), function(author_name) {
  author <- AUTHOR <- author_name
  method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
  method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

  is.pb <- F
  method_wise <- ifelse(is.pb==T, "Patient-wise", "Cell-wise")

  nulls <- all_null_plots(author, FALSE, algos = method_patientwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "p-value", xtype="Cell-wise")

  pows <- all_pow_plots(author, FALSE, algos = method_patientwise, only_tibble = T) %>%
    do.call("rbind", .) %>% dplyr::mutate(ytype = "p-value", xtype="Cell-wise")


  pows$observed_p_value[pows$observed_p_value == 0] = min(pows$observed_p_value[pows$observed_p_value != 0])

  p1 = nulls %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = "p-value", y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    facet_grid(paste0(n.samples, " samples")~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme() +
    ggtitle(paste0(author_name, " dataset"))

  p2 = pows %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = "-log10 p-value", y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    facet_grid(paste0(n.samples, " samples")~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(transform = "log10") +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme() +
    ggtitle(paste0(author_name, " dataset"))

  p1
  p2

  list(p1, p2)
})

plots <- unlist(plots, recursive = FALSE)

pwise_pvalues_plot <- patchwork::wrap_plots(plots, ncol = 2, nrow = 4) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
pwise_pvalues_plot


# Boxplot andf Kolmogorov CellWise ####
res <- readRDS("nullpower/final_res/results.rds")

plots <- lapply(c("bca", "yazar", "kumar", "hsc"), function(author_name) {
  a <- author <- AUTHOR <- author_name

  boxplot <- res %>%
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
    theme(legend.position = "none") +
    ggtitle(paste0(author_name, " dataset"))

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

  list(boxplot, ks_false)
})

plots <- unlist(plots, recursive = FALSE)

boxplots_cwise <- patchwork::wrap_plots(plots, ncol = 2, nrow = 4, widths = c(1,.4)) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
boxplots_cwise

# Boxplot andf Kolmogorov PatientWise ####

plots <- lapply(c("bca", "yazar", "kumar", "hsc"), function(author_name) {
  a <- author <- AUTHOR <- author_name

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
    theme(legend.position = "none") +
    ggtitle(paste0(author_name, " dataset"))

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


  list(ptrue, ks_true)
})

plots <- unlist(plots, recursive = FALSE)

boxplots_pwise <- patchwork::wrap_plots(plots, ncol = 2, nrow = 4, widths = c(1,.4)) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
boxplots_pwise


# Timing - ALL ####
plots <- lapply(c("bca", "yazar", "kumar", "hsc"), function(author_name) {
  timing_res <- readRDS(paste0("nullpower/timing_results/", author_name,".rds"))
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
    scale_y_continuous(transform = 'log10') +
    ggtitle(paste0(author_name, " dataset"))
  timing_plot
})

timings <- patchwork::wrap_plots(plots, ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
timings

ggsave("../figures/supp_03_a.pdf", cellwise_pvalues_plot, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_b.pdf", boxplots_cwise, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_c.pdf", pwise_pvalues_plot, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_d.pdf", boxplots_pwise, width = 14, height = 12, units = "in", dpi = 600)
ggsave("../figures/supp_03_e.pdf", timings, width = 14, height = 12, units = "in", dpi = 600)
