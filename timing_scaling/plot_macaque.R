
rm(list = ls())
require(tidyverse)
library(patchwork)
source("utils_img.R")

results_folder <- "results/MacaqueBrain/"
fits_folder = "results/MacaqueBrain/fits/"
results = get_results(results_folder)

# Comparisons ####
pA <- time_comparison(results) + ggtitle("Time comparison")
pB <- memory_comparison(results) + ggtitle("Memory comparison")
pC <- time_per_gene_comparison(results) + ggtitle("Time per gene comparison")

pD <- time_comparison(results, ratio = "devil (GPU)") + ggtitle("Time ratio comparison")
pE <- memory_comparison(results, ratio = "devil (GPU)") + ggtitle("Memory ratio comparison")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)
corr_plots$lfc
corr_plots$theta

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = 1, pval_cut = .05)
upset = ggplotify::as.ggplot(upset)

design <- "
AABBCC
AABBCC
DDEEFF
DDEEFF
GHHHLL
GHHHLL
"


free(pA) + free(pB) + free(pC) + free(pD) + free(pE) + free(corr_plots$lfc) + guide_area() + free(upset) + free(corr_plots$theta) +
  plot_layout(design = design, guides = "collect")


ggsave("img/report_MacaqueBrain.pdf", plot = last_plot(), width = 15, height = 12, units = 'in')
