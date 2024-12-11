
rm(list = ls())
require(tidyverse)
require(patchwork)
library(devil)
source("utils_img.R")

results_folder <- "results/baronPancreas/"
fits_folder = "results/baronPancreas/fits/"
results = get_results(results_folder)


results <- results %>% dplyr::filter(n_genes * n_cells != 100 * 500)

# Comparisons ####
plots <- plot_time_and_memory_comparison(results, n_extrapolation = 3)
plots$p_time
plots$p_time_ratio
plots$p_mem
plots$p_mem_ratio

# pA <- time_comparison(results) + ggtitle("Time comparison")
# pB <- memory_comparison(results) + ggtitle("Memory comparison")
# pC <- time_per_gene_comparison(results) + ggtitle("Time per gene comparison")
#
# pD <- time_comparison(results, ratio = "devil (GPU)") + ggtitle("Time ratio comparison")
# pE <- memory_comparison(results, ratio = "devil (GPU)") + ggtitle("Memory ratio comparison")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)
corr_plots$lfc
corr_plots$theta

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = 1, pval_cut = .05)
upset = ggplotify::as.ggplot(upset)

# Plot devil volcano
de_fit <- readRDS("results/baronPancreas/fits/cpu_devil_1000_ngene_4000_ncells_2_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))
devil::plot_volcano(de_fit, labels = F)

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


ggsave("img/report_baronPancreas.pdf", plot = last_plot(), width = 15, height = 12, units = 'in')
