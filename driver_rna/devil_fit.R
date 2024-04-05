#!/usr/bin/env Rscript

### Devil testing ###

setwd("D:/PhD_AI/sc_devil/")

library(devil)
library(tidyverse)

load("data/metadata.R")
load("data/rna_counts.R")

#Data filtering for donor ID and cell type
metadata_filtered <- metadata %>%
  filter(donor_id %in% c("Donor_1"),
         cell_type %in% c("retinal progenitor cell", "retinal ganglion cell"),
         sequencing_platform %in% c("Illumina NovaSeq 6000"))

metadata <- metadata_filtered %>% select("cell_type") %>%
  mutate(cell_clusters = case_when(
    cell_type == "retinal progenitor cell"  ~ '0',
    #cell_type == "retinal rod cell"  ~ '1',
    cell_type == "retinal ganglion cell"  ~ '1')) %>%
  select("cell_clusters")

#Create design matrix
design_matrix <- model.matrix(~ seurat_clusters, data = metadata)

### Parameters inference ###

devil_fit_retina <- devil::fit_devil(input_matrix = rna_counts,
                             design_matrix = design_matrix,
                             overdispersion = TRUE,
                             offset=0,
                             size_factors=TRUE,
                             verbose=TRUE,
                             max_iter=500,
                             eps=1e-4,
                             parallel = TRUE)

### Statistical test ###
contrast_vector <- c(0,1)

stat_test_res <- devil::test_de(devil.fit = devil_fit_retina,
                                contrast = contrast_vector,
                                pval_adjust_method = "BH",
                                max_lfc = 10,
                                clusters = NULL)
