
rm(list = ls())
require(patchwork)
require(ggplot2)

umap = readRDS("multiomics_analysis/plot/umap_both.rds")
umap = ggplotify::as.ggplot(umap)

venn <- readRDS("multiomics_analysis/plot/venn_plot_v2.rds") + 
  theme(legend.position = "none")
volcanos <- readRDS("multiomics_analysis/plot/volcanos.rds") +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
go_plots <- readRDS("multiomics_analysis/plot/enrichment_dotplot.RDS") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
simp_plot <- readRDS("multiomics_analysis/plot/simp_plot.rds") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
hm = readRDS("multiomics_analysis/plot/hm.rds")
hm = ggplotify::as.ggplot(hm)



design <- "
AAABBBB
AAABBBB
AAABBBB
CCCCCCC
CCCCCCC
DDFFFFF
DDFFFFF
EEFFFFF
EEFFFFF
EEFFFFF
"

final_plot <- free(umap) + free(hm) +
  free(volcanos) + free(venn) +
  free(simp_plot) + free(go_plots) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

l = 15
x_ratio = 1
y_ratio = 1.2
ggsave("figures/main_muscle.pdf", width = x_ratio * l, height = y_ratio * l, units = "in", dpi = 400, plot = final_plot)
