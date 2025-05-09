
rm(list = ls())
require(tidyverse)
library(patchwork)
source("utils_img.R")


# MacaqueBrain ####
results_folder <- "results/MacaqueBrain/"
fits_folder = "results/MacaqueBrain/fits/"
results = get_results(results_folder)

table <- results %>%
  dplyr::mutate(memory = as.numeric(memory * 1e-9)) %>%
  dplyr::group_by(n_genes, n_cells, model_name) %>%
  dplyr::summarise(time = mean(time), memory=mean(memory)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_genes, n_cells) %>%
  dplyr::mutate(time_ratio = time / time[model_name == 'devil (GPU)']) %>%
  dplyr::mutate(memory_ratio = memory / memory[model_name == 'devil (GPU)']) %>%
  dplyr::select(model_name, n_genes, time, time_ratio, memory, memory_ratio) %>%
  tidyr::pivot_wider(names_from = model_name, values_from = c(time, time_ratio, memory, memory_ratio))
table <- table[,!grepl("ratio_devil (GPU)", colnames(table), fixed = T)]
table <- table[,!grepl("ratio", colnames(table), fixed = T)]
xtable::xtable(table, caption = "Macauqe Brain", label = "tab:MacaqueBrain", digits = 2)

# HumanBlood ####
results_folder <- "results/HumanBlood/"
fits_folder = "results/HumanBlood/fits/"
results = get_results(results_folder)

table <- results %>%
  dplyr::mutate(memory = as.numeric(memory * 1e-9)) %>%
  dplyr::group_by(n_genes, n_cells, model_name) %>%
  dplyr::summarise(time = mean(time), memory=mean(memory)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_genes, n_cells) %>%
  dplyr::mutate(time_ratio = time / time[model_name == 'devil (GPU)']) %>%
  dplyr::mutate(memory_ratio = memory / memory[model_name == 'devil (GPU)']) %>%
  dplyr::select(model_name, n_genes, time, time_ratio, memory, memory_ratio) %>%
  tidyr::pivot_wider(names_from = model_name, values_from = c(time, time_ratio, memory, memory_ratio))
table <- table[,!grepl("ratio_devil (GPU)", colnames(table), fixed = T)]
table <- table[,!grepl("ratio", colnames(table), fixed = T)]
xtable::xtable(table, caption = "HumanBlood", label = "tab:HumanBlood", digits = 2)

# baronPancreas ####
results_folder <- "results/baronPancreas/"
fits_folder = "results/baronPancreas/fits/"
results = get_results(results_folder)
results <- results %>% dplyr::filter(n_genes * n_cells != 100 * 500)

table <- results %>%
  dplyr::mutate(memory = as.numeric(memory * 1e-9)) %>%
  dplyr::group_by(n_genes, n_cells, model_name) %>%
  dplyr::summarise(time = mean(time), memory=mean(memory)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_genes, n_cells) %>%
  dplyr::mutate(time_ratio = time / time[model_name == 'devil (GPU)']) %>%
  dplyr::mutate(memory_ratio = memory / memory[model_name == 'devil (GPU)']) %>%
  dplyr::select(model_name, n_genes, time, time_ratio, memory, memory_ratio) %>%
  tidyr::pivot_wider(names_from = model_name, values_from = c(time, time_ratio, memory, memory_ratio))
table <- table[,!grepl("ratio_devil (GPU)", colnames(table), fixed = T)]
table <- table[,!grepl("ratio", colnames(table), fixed = T)]
xtable::xtable(table, caption = "baronPancreas", label = "tab:baronPancreas", digits = 2)
