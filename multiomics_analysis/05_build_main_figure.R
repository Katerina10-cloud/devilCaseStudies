
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "patchwork", "ComplexHeatmap")
sapply(pkgs, require, character.only = TRUE)

venn <- readRDS("plot/venn_plot.rds")
volcanos <- readRDS("plot/volcanos.rds")
go_plots <- readRDS("plot/enrichment_dotplot.RDS")

design <- "
BBBC
BBBC
#DDD
#DDD
#DDD"

design <- "
BBBBBBCC
BBBBBBCC
#DDDDDDD
#DDDDDDD
#DDDDDDD"

final_plot <- free(volcanos) +
  free(venn) +
  free(go_plots) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("C", 'D', "E"))) &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold')
  )
final_plot

x = 9
ggsave("plot/main_figure_bottom.pdf", width = 1.43 * x, height = x, dpi = 600, units = "in", plot = final_plot)
