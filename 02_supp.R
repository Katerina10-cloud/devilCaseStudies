
rm(list = ls())
require(tidyverse)
require(patchwork)
library(devil)
source("timing_scaling/utils_img.R")


# baronPancreas ####

lfc_cut <- 1
pval_cut <- .05

results_folder <- "timing_scaling/results/baronPancreas/"
fits_folder = "timing_scaling/results/baronPancreas/fits/"
results = get_results(results_folder)
results <- results %>% dplyr::filter(n_genes * n_cells != 100 * 500)

# Comparisons ####
plots <- plot_time_and_memory_comparison(results, n_extrapolation = 3)
plots$p_time <- plots$p_time + theme(legend.position = "right", legend.box = "vertical")
plots$p_time_ratio <- plots$p_time_ratio+ theme(legend.position = "right", legend.box = "vertical")
plots$p_mem <- plots$p_mem + theme(legend.position = "right", legend.box = "vertical")
plots$p_mem_ratio <- plots$p_mem_ratio + theme(legend.position = "right", legend.box = "vertical")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = lfc_cut, pval_cut = pval_cut)

# Plot devil volcano
de_fit <- readRDS("timing_scaling/results/baronPancreas/fits/cpu_devil_1000_ngene_4000_ncells_2_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))

volcano <- devil::plot_volcano(
  de_fit %>% dplyr::filter(adj_pval > 0),
  labels = F,
  lfc_cut = lfc_cut,
  pval_cut = pval_cut, title = "devil")


design <- "
AABB
AABB
AABB
CCDD
CCDD
CCDD
EEFF
EEFF
GGHH
GGHH
GGHH
"

main <- free(plots$p_time) +
  free(plots$p_time_ratio) +
  free(plots$p_mem) +
  free(plots$p_mem_ratio) +
  free(corr_plots$lfc) +
  free(corr_plots$theta) +
  free(volcano) +
  free(upset) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold')
  )
# main
ggsave("figures/supp_02_a.pdf", plot = main, width = 10, height = 12, units = 'in')


# Human blood ####
results_folder <- "timing_scaling/results/HumanBlood/"
fits_folder = "timing_scaling/results/HumanBlood/fits/"
results = get_results(results_folder)

# Comparisons ####
plots <- plot_time_and_memory_comparison(results, n_extrapolation = 3)
plots$p_time <- plots$p_time + theme(legend.position = "right", legend.box = "vertical")
plots$p_time_ratio <- plots$p_time_ratio+ theme(legend.position = "right", legend.box = "vertical")
plots$p_mem <- plots$p_mem + theme(legend.position = "right", legend.box = "vertical")
plots$p_mem_ratio <- plots$p_mem_ratio + theme(legend.position = "right", legend.box = "vertical")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = lfc_cut, pval_cut = pval_cut)

# Plot devil volcano
de_fit <- readRDS("timing_scaling/results/HumanBlood/fits/cpu_devil_1000_ngene_1e+06_ncells_3_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))

volcano <- devil::plot_volcano(
  de_fit %>% dplyr::filter(adj_pval > 0, abs(lfc) < 10),
  labels = F,
  lfc_cut = lfc_cut,
  pval_cut = pval_cut, title = "devil")


design <- "
AABB
AABB
AABB
CCDD
CCDD
CCDD
EEFF
EEFF
GGHH
GGHH
GGHH
"

main <- free(plots$p_time) +
  free(plots$p_time_ratio) +
  free(plots$p_mem) +
  free(plots$p_mem_ratio) +
  free(corr_plots$lfc) +
  free(corr_plots$theta) +
  free(volcano) +
  free(upset) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold')
  )
# main
ggsave("figures/supp_02_b.pdf", plot = main, width = 10, height = 12, units = 'in')
