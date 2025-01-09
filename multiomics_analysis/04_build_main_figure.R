
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ComplexHeatmap")
sapply(pkgs, require, character.only = TRUE)

hm <- readRDS("plot/hm.rds")
venn <- readRDS("plot/venn_plot.rds")
volcanos <- readRDS("plot/volcanos.rds")
go_plots <- readRDS("plot/enrichment_dotplot.RDS")

design <- "
#AAA
BBBC
DDDD"

free(hm) +
  free(volcanos) +
  free(venn) +
  free(go_plots) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 9),
    plot.tag = element_text(face = 'bold')
  )
