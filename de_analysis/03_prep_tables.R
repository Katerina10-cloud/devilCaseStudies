
rm(list = ls())
require(tidyverse)

method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

res <- readRDS("nullpower/final_res/results.rds")
res <- res %>% 
  dplyr::select(name, MCC, F1, FPR, TPR, is.pb, author, patients, ngenes)

# BCA ####
A <- "bca"
## Cell-wise ####
cw_table <- res %>% 
  dplyr::filter(is.pb == FALSE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_cellwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

## PatientWise ####
pw_table <- res %>% 
  dplyr::filter(is.pb == TRUE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_patientwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

full_table <- dplyr::bind_rows(cw_table, pw_table) %>% 
  dplyr::select(!author) %>% 
  dplyr::mutate(is.pb = ifelse(is.pb, "PW", "CW"))
xtable::xtable(full_table, caption = A, label = paste0("tab:", A))


# KUMAR ####
A <- "kumar"
## Cell-wise ####
cw_table <- res %>% 
  dplyr::filter(is.pb == FALSE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_cellwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

## PatientWise ####
pw_table <- res %>% 
  dplyr::filter(is.pb == TRUE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_patientwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

full_table <- dplyr::bind_rows(cw_table, pw_table) %>% 
  dplyr::select(!author) %>% 
  dplyr::mutate(is.pb = ifelse(is.pb, "PW", "CW"))
xtable::xtable(full_table, caption = A, label = paste0("tab:", A))

# YAZAR ####
A <- "yazar"
## Cell-wise ####
cw_table <- res %>% 
  dplyr::filter(is.pb == FALSE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_cellwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

## PatientWise ####
pw_table <- res %>% 
  dplyr::filter(is.pb == TRUE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_patientwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

full_table <- dplyr::bind_rows(cw_table, pw_table) %>% 
  dplyr::select(!author) %>% 
  dplyr::mutate(is.pb = ifelse(is.pb, "PW", "CW"))
xtable::xtable(full_table, caption = A, label = paste0("tab:", A))

# HSC ####
A <- "hsc"
## Cell-wise ####
cw_table <- res %>% 
  dplyr::filter(is.pb == FALSE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_cellwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

## PatientWise ####
pw_table <- res %>% 
  dplyr::filter(is.pb == TRUE) %>% 
  dplyr::filter(author == A) %>% 
  dplyr::filter(name %in% method_patientwise) %>% 
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>% 
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>% 
  dplyr::group_by(name, author, patients, ngenes, is.pb) %>% 
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", mean(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = c(MCC))

full_table <- dplyr::bind_rows(cw_table, pw_table) %>% 
  dplyr::select(!author) %>% 
  dplyr::mutate(is.pb = ifelse(is.pb, "PW", "CW"))
xtable::xtable(full_table, caption = A, label = paste0("tab:", A))
