
rm(list = ls())
require(tidyverse)
library(patchwork)
source("utils_img.R")

dataset_name = "MacaqueBrain"

add_predicttions = TRUE
results_folder <- paste0("results/", dataset_name)
fits_folder = paste0("results/", dataset_name, "/fits/")
time_results = get_time_results(results_folder)
mem_results = get_memory_results(results_folder)
if (add_predicttions) {
  time_results = plyr::rbind.fill(time_results, predict_time_results(time_results))
  mem_results = plyr::rbind.fill(mem_results, predict_memory_results(mem_results))
  time_results$type[is.na(time_results$type)] = "observed"
  mem_results$type[is.na(mem_results$type)] = "observed"
} else {
  time_results$type = "observed"
  mem_results$type = "observed"
}

# Comparisons ####
pA <- time_comparison(time_results) + guides(linetype = "none")
pB <- memory_comparison(mem_results) + guides(linetype = "none")

pC <- time_comparison(time_results, ratio = "glmGamPoi / cpu") + guides(linetype = "none")
pD <- memory_comparison(mem_results, ratio = "glmGamPoi / cpu") + guides(linetype = "none")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)
corr_plots$lfc
corr_plots$theta

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = 1, pval_cut = .05)
upset = ggplotify::as.ggplot(upset)

design <- "
AAABBB
AAABBB
CCCDDD
CCCDDD
GHHHLL
GHHHLL
"

free(pA) + free(pB) + free(pC) + free(pD) + free(corr_plots$lfc) + free(corr_plots$theta) + free(upset) +
  plot_layout(design = design)

ggsave(paste0("img/",dataset_name,".pdf"), plot = last_plot(), width = 15, height = 12, units = 'in')
