rm(list = ls())
require(tidyverse)
require(patchwork)

dir.create("all_figures/scaling_and_sim/", recursive = T)

# MAIN ####
# Use MacaqueBrain and HSC
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/speedup.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/memory_ratio.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pAB = (pA + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pB + theme(legend.position = "none")) +
  plot_layout(guides = "collect")

pCD = readRDS("de_analysis/img/RDS/hsc/pvalues.RDS")
pEF = readRDS("de_analysis/img/RDS/hsc/MCCs_boxplots.RDS") + theme(legend.position = "none")

design = "
AAAA
AAAA
AAAA
AAAA
AAAA
##CC
##CC
##DD
##DD
EEEE
EEEE
EEEE
EEEE"

final_plot = free(pAB) +
  free(pCD$null_pvalue) + free(pCD$de_pvalue) +
  free(pEF) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A", "", "C", "D", "E"))) &
  theme(
    text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = 'bold'),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )
final_plot
ggsave(filename = "all_figures/scaling_and_sim/main_2_v0.pdf", plot = final_plot, dpi = 600, width = 11.7, height = 11.7, units = "in")
rm(pA, pB, pCD, pEF, final_plot, design)

# EXTENDED ####
# Use MacaqueBrain and HSC
pAB = readRDS("timing_scaling/img/RDS/MacaqueBrain/correlation.RDS")
pC = readRDS("timing_scaling/img/RDS/MacaqueBrain/large.RDS")
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
FFGGHH
FFGGHH"

final_plot = free(pAB$lfc) + free(pAB$theta) + free(pC) +
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
ggsave(filename = "all_figures/scaling_and_sim/ext_2.pdf", plot = final_plot, dpi = 600, width = 11.7, height = 7, units = "in")
rm(pAB, pC, pD, pFG, pH, final_plot, design)

# SUPP SCALING ####
## Times and Memory ####
pA = readRDS("timing_scaling/img/RDS/MacaqueBrain/runtime.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pB = readRDS("timing_scaling/img/RDS/MacaqueBrain/memory.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pAB = (pA + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pB + theme(legend.position = "none")) +
  plot_layout(guides = "collect")

pC = readRDS("timing_scaling/img/RDS/MacaqueBrain/disp_speedup.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pD = readRDS("timing_scaling/img/RDS/MacaqueBrain/disp_runtime.RDS") +
  theme(legend.direction='horizontal', legend.position = "bottom", legend.box = "vertical")
pCD = (pC + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pD + theme(legend.position = "none")) +
  plot_layout(guides = "collect")

final_plot = free(pAB) + free(pCD) +
  plot_layout(nrow = 2) +
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
ggsave(filename = "all_figures/scaling_and_sim/supp_scaling_times_and_memory.pdf", plot = final_plot, dpi = 600, width = 11.7, height = 7, units = "in")
rm(pAB, pCD, final_plot, pA, pB, pC, pD)


## Small dataset ####
dataset_name = "baronPancreas"

path_to_rds = file.path("timing_scaling/img/RDS/", dataset_name)
pA = readRDS(file.path(path_to_rds, "runtime.RDS"))
pB = readRDS(file.path(path_to_rds, "speedup.RDS"))
pAB = (pA + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pB + theme(legend.position = "none")) +
  plot_layout(guides = "collect")
pC = readRDS(file.path(path_to_rds, "memory.RDS"))
pD = readRDS(file.path(path_to_rds, "memory_ratio.RDS"))
pCD = (pC + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
  (pD + theme(legend.position = "none")) +
  plot_layout(guides = "collect")
pEF = readRDS(file.path(path_to_rds, "correlation.RDS"))
pG = readRDS(file.path(path_to_rds, "upset.RDS"))

design = "
AAAAAA
AAAAAA
AAAAAA
CCCCCC
CCCCCC
CCCCCC
EEFFGG
EEFFGG"

supp_fig = free(pAB) + free(pCD) +
  free(pEF$lfc) + free(pEF$theta) + free(pG) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "", "D"))) &
  theme(
    text = element_text(size = 12),
    plot.tag = element_text(face = 'bold'),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(1, "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  )
ggsave(filename = paste0("all_figures/scaling_and_sim/supp_scaling_",dataset_name,".pdf"), plot = supp_fig, dpi = 600, width = 11.7, height = 11.7, units = "in")

# SUPP DE ANALYSIS ####
## All methods ####
all_models_plots = readRDS("de_analysis/img/RDS/all_models.RDS")
pA = all_models_plots$MCC + theme(legend.position = "bottom")
pB = all_models_plots$timing + theme(legend.position = "bottom")
pC = all_models_plots$failure_rate + theme(legend.position = "bottom")

des = "
AAAAA
AAAAA
BBBCC
BBBCC"

all_models_plot = free(pA) + free(pB) + free(pC) +
  plot_layout(design = des) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))
ggsave("all_figures/scaling_and_sim/supp_sim_allmodels.pdf", all_models_plot, width = 10, height = 10, dpi = 600, units = "in")
rm(all_models_plots, pA, pB, pC, all_models_plot, des)


## All datasets ####
for (author in c("hsc", "kumar", "yazar", "bca")) {
  pA = readRDS(file.path("de_analysis/img/RDS/", author, "pvalues.RDS"))
  pA = (pA$null_pvalue + theme(legend.direction='vertical', legend.position = "none", legend.box = "horizontal")) +
    (pA$de_pvalue + theme(legend.position = "none"))

  pBC = readRDS(file.path("de_analysis/img/RDS/", author, "MCCs.RDS"))
  pB = pBC$cellwise + ggtitle("Cell-wise")
  pC = pBC$patientwise + ggtitle("Patient-wise")
  pBC = (pB + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
    (pC + theme(legend.position = "none")) +
    plot_layout(guides = "collect")

  pDE = readRDS(file.path("de_analysis/img/RDS/", author, "ks_test.RDS"))
  pD = pDE$cellwise + ggtitle("Cell-wise")
  pE = pDE$patientwise + ggtitle("Patient-wise")
  pDE = (pD + theme(legend.direction='vertical', legend.position = "bottom", legend.box = "horizontal")) +
    (pE + theme(legend.position = "none")) +
    plot_layout(guides = "collect")

  pF = readRDS(file.path("de_analysis/img/RDS/", author, "timing.RDS"))

  design = "
AAAA
BBBB
CCCC
#DD#"

  final_plot = pA / pBC / pDE / pF +
    plot_layout(design = design) +
    plot_annotation(tag_levels = c("A", "B", "C", "D")) &
    theme(plot.tag = element_text(face = 'bold'))

  ggsave(paste0("all_figures/scaling_and_sim/supp_sim_",author,".pdf"), final_plot, width = 10, height = 10, dpi = 600, units = "in")
}

