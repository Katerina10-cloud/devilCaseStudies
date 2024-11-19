rm(list = ls())
source("scripts/utils.R")
require(tidyverse)

# SMALL ####
results <- get_results("results/baronPancreas/")
plot_comps = plot_time_and_memory_comparison(results)

patchwork::wrap_plots(plot_comps, ncol = 3, guides = "collect")

# Correlations
corr_plots <- plot_correlations("results/baronPancreas/fits/")
patchwork::wrap_plots(corr_plots, nrow = 3, ncol = 1, guides = "collect")

# BIG ####
results <- get_results("results/macaque_brain/")
plot_comps = plot_time_and_memory_comparison(results)

plot_comps$time
patchwork::wrap_plots(plot_comps, ncol = 3, guides = "collect")
ggsave("img/macaque_brain/time_and_memory_comparison.pdf", plot = last_plot(), width = 11, height = 6, units = 'in')

# Correlations
corr_plots <- plot_correlations("results/macaque_brain/fits/")
patchwork::wrap_plots(corr_plots, nrow = 2, ncol = 1, guides = "collect")
ggsave("img/macaque_brain/correlations.pdf", plot = last_plot(), width = 11, height = 6, units = 'in')

# Volcano plots
v1 <- devil:::plot_volcano(
  readRDS("results/macaque_brain/fits/cpu_glmGamPoi_0.25_pgene_0.5_pcells_2_celltypes.rds") %>% 
    dplyr::mutate(name = row_number()) %>% 
    dplyr::filter(pval > 0), labels = FALSE, title = "glmGamPoi volcano plot", lfc_cut = .25)

v2 <- devil:::plot_volcano(
  readRDS("results/macaque_brain/fits/gpu_devil_0.25_pgene_0.5_pcells_2_celltypes.rds") %>% 
    dplyr::mutate(name = row_number()) %>% 
    dplyr::filter(pval > 0), labels = FALSE, title = "GPU devil volcano plot", lfc_cut = .25)

patchwork::wrap_plots(list(v1, v2), ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")
ggsave("img/macaque_brain/volcanos.pdf", plot = last_plot(), width = 11, height = 6, units = 'in')
