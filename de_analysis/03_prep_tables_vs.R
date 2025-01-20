
rm(list = ls())
require(tidyverse)

method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula")
method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)")

res <- readRDS("nullpower/final_res/results.rds")
res <- res %>%
  dplyr::select(name, MCC, F1, FPR, TPR, is.pb, author, patients, ngenes)

## Cell-wise ####
cw_table <- res %>%
  dplyr::filter(is.pb == FALSE) %>%
  dplyr::filter(name %in% method_cellwise) %>%
  dplyr::mutate(name = ifelse(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
  dplyr::mutate(name = ifelse(grepl("Devil", name), "devil", name)) %>%
  dplyr::group_by(name, author, is.pb) %>%
  dplyr::summarise(
    across(
      c(MCC), # Explicitly list the columns
      ~sprintf("%.3f ± %.3f", median(.), sd(.)), # Compute mean ± sd with 3 significant digits
      .names = "{.col}" # Column naming
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = author, values_from = c(MCC))

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
  pivot_wider(names_from = c(author, patients), values_from = c(MCC))
pw_table
