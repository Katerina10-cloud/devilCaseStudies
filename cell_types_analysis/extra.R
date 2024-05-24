#setwd("~/GitHub/cell_types_analysis")
rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "AnnotationDbi", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
library(scMayoMap)
library(Seurat)
source("utils.R")

d_umap <- readRDS("results/BigBloodData/umap_tibble.rds")
d_umap %>%
  ggplot(mapping = aes(x=x, y=y, col=cluster)) +
  geom_point()
