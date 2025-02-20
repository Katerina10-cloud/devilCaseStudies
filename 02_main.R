
rm(list = ls())
require(tidyverse)
require(patchwork)
library(ggupset)
library(devil)
source("timing_scaling/utils_img.R")

lfc_cut <- 1
pval_cut <- .05

results_folder <- "timing_scaling/results/MacaqueBrain/"
fits_folder = "timing_scaling/results/MacaqueBrain/fits/"
results = get_results(results_folder)

# Comparisons ####
p <- plot_time_and_memory_comparison(results, n_extrapolation = 3, ncols = 4) + 
  theme(
    legend.position = "bottom", legend.box = "vertical", legend.spacing = unit(-10, "pt"),
    axis.title.y = element_blank()
  )
p
# plots$p_time <- plots$p_time + theme(legend.position = "right", legend.box = "vertical")
# plots$p_time_ratio <- plots$p_time_ratio+ theme(legend.position = "right", legend.box = "vertical")
# plots$p_mem <- plots$p_mem + theme(legend.position = "right", legend.box = "vertical")
# plots$p_mem_ratio <- plots$p_mem_ratio + theme(legend.position = "right", legend.box = "vertical")

# Correlations ####
corr_plots <- plot_correlations(fits_folder)
corr_plots$lfc
corr_plots$theta

# Plot UpSet ####
upset <- plot_upset(fits_folder, lfc_cut = lfc_cut, pval_cut = pval_cut)

# Plot devil volcano
de_fit <- readRDS("timing_scaling/results/MacaqueBrain/fits/cpu_devil_1000_ngene_1e+06_ncells_2_celltypes.rds") %>%
  dplyr::mutate(name = paste0("Gene ", row_number()))

library(biomaRt)
# Connect to human dataset
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Connect to macaque dataset
macaque = useMart("ensembl", dataset="mmulatta_gene_ensembl")

# Get list of human IDs
human_genes <- rownames(de_fit) # add your list of genes here

# Get human symbols and macaque orthologs
df_macaque_genes = getBM(attributes=c("ensembl_gene_id", "external_gene_name", 
                   "mmulatta_homolog_ensembl_gene", "mmulatta_homolog_associated_gene_name"),
      filters="ensembl_gene_id", 
      values=human_genes,
      mart=human)

# GABAergic neuron markers (most reliable)
gabaergic_core <- c("GABRB2", "ERBB4")

# Glutamatergic neuron markers (most reliable)
glutamatergic_core <- c("GRIA2", "SNAP25", "CAMK2A")

genes_of_interest = c(gabaergic_core, glutamatergic_core)

de_fit$ensembl_gene_id = rownames(de_fit)
volcano = de_fit %>% 
  dplyr::left_join(df_macaque_genes, by="ensembl_gene_id") %>% 
  dplyr::select(pval, adj_pval, lfc, mmulatta_homolog_associated_gene_name) %>% 
  dplyr::mutate(name = mmulatta_homolog_associated_gene_name) %>% 
  #dplyr::mutate(adj_pval = ifelse(adj_pval != 0, adj_pval, min(adj_pval[adj_pval != 0]))) %>% 
  #dplyr::filter(adj_pval <= .05, abs(lfc) > .5) %>% 
  dplyr::mutate(name = ifelse(name %in% genes_of_interest, name, "")) %>% 
  plot_volc(
    labels = TRUE,
    lfc_cut = .5,
    pval_cut = pval_cut, 
    title = "GABAergic v Glutamatergic neurons"
  )
volcano

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
ggsave("figures/main_02.pdf", plot = main, width = 12, height = 8, units = 'in')
#ggsave("figures/main_02.pdf", plot = main, width = 10, height = 12, units = 'in')

# design <- "
# AABB
# AABB
# AABB
# CCDD
# CCDD
# CCDD
# EEFF
# EEFF
# GGHH
# GGHH
# GGHH
# "
# 
# main <- free(plots$p_time) +
#   free(plots$p_time_ratio) +
#   free(plots$p_mem) +
#   free(plots$p_mem_ratio) +
#   free(corr_plots$lfc) +
#   free(corr_plots$theta) +
#   free(volcano) +
#   free(upset) +
#   plot_layout(design = design) +
#   plot_annotation(tag_levels = "A") &
#   theme(
#     text = element_text(size = 12),
#     plot.tag = element_text(face = 'bold')
#   )
# # main
# ggsave("figures/main_02.pdf", plot = main, width = 10, height = 12, units = 'in')
