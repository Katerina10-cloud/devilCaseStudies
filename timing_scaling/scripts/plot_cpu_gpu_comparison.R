rm(list = ls())
require(tidyverse)

# SMALL ####

results_folder <- "results/baronPancreas/"
results_paths <- list.files(results_folder)
results_paths <- results_paths[results_paths != "fits"]
p <- results_paths[1]
results <- lapply(results_paths, function(p) {
  print(p)
  info <- unlist(strsplit(p, "_"))
  res = readRDS(paste0(results_folder, p))
  
  
  dplyr::tibble(
    model_name = paste(info[1], info[2]),
    p_gene = info[3],
    p_cells = info[5],
    n_cell_types = info[7],
    time = as.numeric(unlist(res$time), units = "secs"),
    memory = res$mem_alloc
  )
}) %>% do.call("bind_rows", .) %>%
  dplyr::mutate(p_cells = as.numeric(p_cells))

results %>%
  ggplot(mapping = aes(x = p_cells, y = time, col = model_name)) +
  geom_point() +
  geom_smooth() +
  ggh4x::facet_nested(~"Gene percentage"+p_gene) +
  theme_bw() +
  labs(x = "Percentage of cells", y = "Time (s)", col = "Model")

results %>%
  ggplot(mapping = aes(x = factor(p_cells), y = memory * 1e-9, fill=model_name)) +
  geom_col(color="black", position=position_dodge()) +
  ggh4x::facet_nested(~"Gene percentage"+p_gene) +
  theme_bw() +
  labs(x = "Cells percentage", y = "Memory (GB)", fill = "Model")

results %>%
  dplyr::group_by(p_gene, p_cells, n_cell_types) %>%
  # dplyr::mutate(n = n()) %>%
  # dplyr::filter(n == 3) %>%
  dplyr::mutate(ratio_time = time / time[model_name == "gpu devil"]) %>%
  ggplot(mapping = aes(x = p_cells, y = ratio_time, col = model_name)) +
  geom_point() +
  geom_smooth() +
  ggh4x::facet_nested(~"Gene percentage"+p_gene) +
  theme_bw() +
  labs(x = "Percentage of cells", y = "Time (s)", col = "Model")

# Correlations
gpu.devil.res <- readRDS("results/baronPancreas/fits/gpu_devil_1_pgene_1_pcells_2_celltypes.rds")
devil.res <- readRDS("results/baronPancreas/fits/cpu_devil_1_pgene_1_pcells_2_celltypes.rds")
glm.res <- readRDS("results/baronPancreas/fits/cpu_glmGamPoi_1_pgene_1_pcells_2_celltypes.rds")

my_theme = function(p) {
  p +
    labs(x="devil", y="glmGamPoi") +
    theme_bw()
}


hist(devil.res$lfc)
hist(glm.res$lfc)
hist(gpu.devil.res$lfc)

dplyr::tibble(x = devil.res$lfc, y = glm.res$lfc) %>% 
  dplyr::filter(abs(x) < 10 & abs(y) < 10) %>% 
  ggplot(mapping = aes(x=x, y=y)) +
  geom_point() +
  ggtitle("LFC")
 
dplyr::tibble(x = gpu.devil.res$lfc, y = glm.res$lfc) %>% 
  dplyr::filter(abs(x) < 10 & abs(y) < 10) %>% 
  ggplot(mapping = aes(x=x, y=y)) +
  geom_point() +
  ggtitle("LFC")
my_theme(p)

dplyr::tibble(x = devil.res$theta, y = glm.res$theta) %>% 
  ggplot(mapping = aes(x=x, y=y)) +
  geom_point() +
  theme_bw() +
  ggtitle(bquote(Overdispersion ~ theta))

dplyr::tibble(x = gpu.devil.res$theta, y = glm.res$theta) %>% 
  ggplot(mapping = aes(x=x, y=y)) +
  geom_point() +
  theme_bw() +
  ggtitle(bquote(Overdispersion ~ theta))

dplyr::tibble(x = devil.res$pval, y = glm.res$pval) %>% 
  ggplot(mapping = aes(x=x, y=y)) +
  geom_point()

dplyr::tibble(x = devil.res$adj_pval, y = glm.res$adj_pval) %>% 
  ggplot(mapping = aes(x=x, y=y)) +
  geom_point()


# BIG ####

results_folder <- "results/macaque_brain/"
results_paths <- list.files(results_folder)
results_paths <- results_paths[results_paths != "fits"]

results <- lapply(results_paths, function(p) {
  print(p)
  info <- unlist(strsplit(p, "_"))
  res = readRDS(paste0(results_folder, p))

  dplyr::tibble(
    model_name = paste(info[1], info[2]),
    p_gene = info[3],
    p_cells = info[5],
    n_cell_types = info[7],
    time = as.numeric(res$median, units = "secs"),
    memory = res$mem_alloc
  )
}) %>% do.call("bind_rows", .) %>%
  dplyr::mutate(p_cells = as.numeric(p_cells))

results %>%
  ggplot(mapping = aes(x = p_cells, y = time, col = model_name)) +
  geom_point() +
  geom_line() +
  ggh4x::facet_nested(~"Gene percentage"+p_gene) +
  theme_bw() +
  labs(x = "Percentage of cells", y = "Time (s)", col = "Model")

results %>%
  ggplot(mapping = aes(x = factor(p_cells), y = memory * 1e-9, fill=model_name)) +
  geom_col(color="black", position=position_dodge()) +
  ggh4x::facet_nested(~"Gene percentage"+p_gene) +
  theme_bw() +
  labs(x = "Cells percentage", y = "Memory (GB)", fill = "Model")

results %>%
  dplyr::group_by(p_gene, p_cells, n_cell_types) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n == 3) %>%
  dplyr::mutate(ratio_time = time / time[model_name == "gpu devil"]) %>%
  ggplot(mapping = aes(x = p_cells, y = ratio_time, col = model_name)) +
  geom_point() +
  geom_line() +
  ggh4x::facet_nested(~"Gene percentage"+p_gene) +
  theme_bw() +
  labs(x = "Percentage of cells", y = "Time (s)", col = "Model")

# Correlations
devil.beta <- readRDS("results/macaque_brain/fits/beta_cpu_devil_0.25_pgene_1_pcells_2_celltypes.rds")
glm.beta <- readRDS("results/macaque_brain/fits/beta_cpu_glmGamPoi_0.25_pgene_1_pcells_2_celltypes.rds")

devil.beta %>% is.na() %>% sum()
