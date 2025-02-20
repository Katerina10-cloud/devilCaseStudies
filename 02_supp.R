
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
p <- plot_time_and_memory_comparison(results, n_extrapolation = 3, ncols = 4) + 
  theme(
    legend.position = "bottom", legend.box = "vertical", legend.spacing = unit(-10, "pt"),
    axis.title.y = element_blank()
  )
# plots$p_time <- plots$p_time + theme(legend.position = "right", legend.box = "vertical")
# plots$p_time_ratio <- plots$p_time_ratio+ theme(legend.position = "right", legend.box = "vertical")
# plots$p_mem <- plots$p_mem + theme(legend.position = "right", legend.box = "vertical")
# plots$p_mem_ratio <- plots$p_mem_ratio + theme(legend.position = "right", legend.box = "vertical")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = lfc_cut, pval_cut = pval_cut)

# Plot devil volcano
de_fit <- readRDS("timing_scaling/results/baronPancreas/fits/cpu_devil_1000_ngene_4000_ncells_2_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))

volcano <- plot_volc(
  de_fit,
  labels = F,
  lfc_cut = lfc_cut,
  pval_cut = pval_cut, title = "devil")


design <- "
AAAAA
AAAAA
BBCCD
"

# design <- "
# AB
# AC
# AD
# AD
# AE
# "

main <- free(p) +
  free(corr_plots$lfc + theme(title = element_blank())) +
  free(corr_plots$theta + theme(title = element_blank())) +
  #free(volcano) +
  free(upset) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold')
  )
main
ggsave("figures/supp_02_a.pdf", plot = main, width = 10, height = 12, units = 'in')


# Human blood ####
results_folder <- "timing_scaling/results/HumanBlood/"
fits_folder = "timing_scaling/results/HumanBlood/fits/"
results = get_results(results_folder)

# Comparisons ####
p <- plot_time_and_memory_comparison(results, n_extrapolation = 3, ncols = 4)
# plots$p_time <- plots$p_time + theme(legend.position = "right", legend.box = "vertical")
# plots$p_time_ratio <- plots$p_time_ratio+ theme(legend.position = "right", legend.box = "vertical")
# plots$p_mem <- plots$p_mem + theme(legend.position = "right", legend.box = "vertical")
# plots$p_mem_ratio <- plots$p_mem_ratio + theme(legend.position = "right", legend.box = "vertical")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = lfc_cut, pval_cut = pval_cut)

# Plot devil volcano
de_fit <- readRDS("timing_scaling/results/HumanBlood/fits/cpu_devil_1000_ngene_1e+06_ncells_3_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))

de_fit <- readRDS("timing_scaling/results/HumanBlood/fits/cpu_devilS_1000_ngene_1e+06_ncells_3_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))

volcano <- plot_volc(
  de_fit,
  labels = F,
  lfc_cut = lfc_cut,
  pval_cut = pval_cut, title = "devil")


design <- "
AAAAA
AAAAA
BBCCD
"

# design <- "
# AB
# AC
# AD
# AD
# AE
# "

main <- free(p) +
  free(corr_plots$lfc + theme(title = element_blank())) +
  free(corr_plots$theta + theme(title = element_blank())) +
  #free(volcano) +
  free(upset) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold')
  )
main
ggsave("figures/supp_02_b.pdf", plot = main, width = 10, height = 12, units = 'in')
