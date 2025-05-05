
rm(list = ls())
require(tidyverse)

method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

res <- readRDS("nullpower/final_res/results.rds")
# res <- res %>%
#   dplyr::select(name, MCC, F1, FPR, TPR, is.pb, author, patients, ngenes)

## Cell-wise ####
cw_table <- res %>%
  dplyr::filter(is.pb == FALSE) %>%
  dplyr::filter(name %in% method_cellwise) %>%
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>%
  dplyr::group_by(name, author, patients, is.pb) %>%
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", median(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(name, author, patients, MCC) %>% 
  dplyr::mutate(name = factor(name, levels=c("devil", "Nebula", "glmGamPoi", "limma"))) %>% 
  dplyr::arrange(name) %>% 
  tidyr::pivot_wider(names_from = name, values_from = MCC)

cw_table <- res %>%
  dplyr::filter(is.pb == FALSE) %>%
  dplyr::filter(name %in% method_cellwise) %>%
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>%
  dplyr::group_by(name, author, patients, is.pb) %>%
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", quantile(., .25), quantile(., .75)),
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(name, author, patients, MCC) %>% 
  dplyr::mutate(name = factor(name, levels=c("devil", "Nebula", "glmGamPoi", "limma"))) %>% 
  dplyr::arrange(name) %>% 
  tidyr::pivot_wider(names_from = name, values_from = MCC)

cw_table

## PatientWise ####
pw_table <- res %>%
  dplyr::filter(is.pb == TRUE) %>%
  dplyr::filter(name %in% method_patientwise) %>%
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>%
  dplyr::group_by(name, author, patients, is.pb) %>%
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", median(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(name, author, patients, MCC) %>% 
  dplyr::mutate(name = factor(name, levels=c("devil", "Nebula", "glmGamPoi", "limma"))) %>% 
  dplyr::arrange(name) %>% 
  tidyr::pivot_wider(names_from = name, values_from = MCC)
pw_table


# Print tables ####
cw_table %>% xtable::xtable()
pw_table %>% xtable::xtable()


# Prep timing tables
folder_path = "nullpower/timing_results/"
lf <- list.files(folder_path)
df <- lapply(1:length(lf), function(i) {
  readRDS(paste0(folder_path, lf[i]))   %>% 
    dplyr::filter(algo %in% c("Devil (base)", "glmGamPoi (cell)", "Nebula")) %>%
    dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
    dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>% 
    dplyr::select(algo, timings, author, cell_order) %>% 
    dplyr::group_by(algo, author, cell_order) %>% 
    dplyr::summarise(m = median(timings), s = sd(timings))  
}) %>% do.call("bind_rows", .)

df %>% 
  dplyr::rename(name=algo) %>%  
  dplyr::mutate(time = sprintf("%.3f ± %.3f", m, s)) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>%
  dplyr::select(name, author, cell_order, time) %>% 
  tidyr::pivot_wider(names_from = name, values_from = time) %>% 
  xtable::xtable()
  



