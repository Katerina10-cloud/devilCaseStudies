rm(list = ls())
require(tidyverse)

# MAIN ####
# Use MacaqueBrain and HSC
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/speedup.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/memory_ratio.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pCD = readRDS("de_analysis/img/RDS/hsc/pvalues.RDS")
pEF = readRDS("de_analysis/img/RDS/hsc/MCCs.RDS")

design = "
AABB
AABB
AABB
AABB
AABB
##CC
##CC
##DD
##DD
EEFF
EEFF
EEFF
EEFF"

final_plot = free(pA) + free(pB) +
  free(pCD$null_pvalue) + free(pCD$de_pvalue) +
  free(pEF$cellwise) + free(pEF$patientwise) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )
final_plot
ggsave(filename = "figures/main_2_v0.pdf", plot = final_plot, dpi = 600, width = 11.7, height = 11.7, units = "in")
rm(pA, pB, pCD, pEF, final_plot, design)

# EXTENDED ####
# Use MacaqueBrain and HSC
pAB = readRDS("timing_scaling/img/RDS/MacaqueBrain/correlation.RDS")
pC = readRDS("timing_scaling/img/RDS/MacaqueBrain/large.RDS")
pD = readRDS("timing_scaling/img/RDS/MacaqueBrain/disp_speedup.RDS") +
  theme(legend.direction='horizontal',
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(1, "pt"),
        legend.box.margin = margin(0, 0, 0, 0))
pE = readRDS("timing_scaling/img/RDS/MacaqueBrain/disp_runtime.RDS") +
  theme(legend.direction='horizontal',
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(1, "pt"),
        legend.box.margin = margin(0, 0, 0, 0))
pFG = readRDS("de_analysis/img/RDS/hsc/ks_test.RDS")
pFG$cellwise = pFG$cellwise + theme(legend.direction='horizontal',
                     legend.position = "bottom",
                     legend.box = "vertical",
                     legend.spacing.y = unit(0, "pt"),
                     legend.spacing.x = unit(1, "pt"),
                     legend.box.margin = margin(0, 0, 0, 0)) +
  guides(color = guide_legend(ncol = 2))
pFG$patientwise = pFG$patientwise + theme(legend.direction='horizontal',
                                    legend.position = "bottom",
                                    legend.box = "vertical",
                                    legend.spacing.y = unit(0, "pt"),
                                    legend.spacing.x = unit(1, "pt"),
                                    legend.box.margin = margin(0, 0, 0, 0)) +
  guides(color = guide_legend(ncol = 2))

pH = readRDS("de_analysis/img/RDS/hsc/timing.RDS") +
  theme(legend.direction='horizontal',
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(1, "pt"),
        legend.box.margin = margin(0, 0, 0, 0))

design = "
AABBCC
AABBCC
DDDEEE
DDDEEE
DDDEEE
FFGGHH
FFGGHH
FFGGHH"

final_plot = free(pAB$lfc) + free(pAB$theta) + free(pC) +
  free(pD) + free(pE) +
  free(pFG$cellwise) + free(pFG$patientwise) + free(pH) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )
final_plot
ggsave(filename = "figures/ext_2.pdf", plot = final_plot, dpi = 600, width = 11.7, height = 11.7, units = "in")
rm(pAB, pC, pD, pE, pFG, pH, final_plot, design)

# SUPP SCALING ####
dataset_name = "baronPancreas"
for (dataset_name in c("baronPancreas", "MacaqueBrain")) {
  path_to_rds = file.path("timing_scaling/img/RDS/", dataset_name)
  pA = readRDS(file.path(path_to_rds, "runtime.RDS"))
  pB = readRDS(file.path(path_to_rds, "speedup.RDS"))
  pC = readRDS(file.path(path_to_rds, "memory.RDS"))
  pD = readRDS(file.path(path_to_rds, "memory_ratio.RDS"))
  pEF = readRDS(file.path(path_to_rds, "correlation.RDS"))
  pG = readRDS(file.path(path_to_rds, "upset.RDS"))

  design = "
AAABBB
AAABBB
CCCDDD
CCCDDD
EEFFGG"

  supp_fig = free(pA) + free(pC) +
    free(pB) + free(pD) +
    free(pEF$lfc) + free(pEF$theta) + free(pG) +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "A") &
    theme(
      text = element_text(size = 12),
      plot.tag = element_text(face = 'bold'),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      panel.grid.minor = element_blank(),
      legend.direction='horizontal',
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(0, "pt"),
      legend.spacing.x = unit(1, "pt"),
      legend.box.margin = margin(0, 0, 0, 0)
    )

  ggsave(filename = paste0("figures/supp_scaling_",dataset_name,".pdf"), plot = supp_fig, dpi = 600, width = 11.7, height = 11.7, units = "in")
}

# SUPP DE ANALYSIS ####
# NULL
datasets = c("hsc", "bca", "yazar", "kumar")
lapply(datasets, function(p) {
  readRDS(paste0("de_analysis/img/RDS/",p,"/pvalues.RDS"))
})





