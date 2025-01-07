
rm(list = ls())
require(tidyverse)
require(patchwork)
library(ggupset)
library(devil)
source("timing_scaling/utils_img.R")

lfc_cut <- 1
pval_cut <- .05

results_folder <- "timing_scaling/results/MacaqueBrain/"
fits_folder = "timing_scaling/results/MacaqueBrain/fits/"
results = get_results(results_folder)

# Comparisons ####
plots <- plot_time_and_memory_comparison(results, n_extrapolation = 3)
plots$p_time <- plots$p_time + theme(legend.position = "right", legend.box = "vertical")
plots$p_time_ratio <- plots$p_time_ratio+ theme(legend.position = "right", legend.box = "vertical")
plots$p_mem <- plots$p_mem + theme(legend.position = "right", legend.box = "vertical")
plots$p_mem_ratio <- plots$p_mem_ratio + theme(legend.position = "right", legend.box = "vertical")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)
corr_plots$lfc
corr_plots$theta

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = lfc_cut, pval_cut = pval_cut)

# Plot devil volcano
de_fit <- readRDS("timing_scaling/results/MacaqueBrain/fits/cpu_devil_1000_ngene_1e+06_ncells_2_celltypes.rds") %>%
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
ggsave("figures/main_02.pdf", plot = main, width = 10, height = 12, units = 'in')
