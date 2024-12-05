
method_colors = c(
  "glmGamPoi" = "#EAB578",
  "nebula" =  "steelblue",
  "NEBULA" =  "steelblue",
  "devil" = "#099668",
  "devil (GPU)" = "indianred"
)

get_results <- function(results_folder) {
  results_paths <- list.files(results_folder)
  results_paths <- results_paths[results_paths != "fits"]

  results <- lapply(results_paths, function(p) {
    print(p)
    info <- unlist(strsplit(p, "_"))
    res = readRDS(paste0(results_folder, p))


    dplyr::tibble(
      model_name = paste(info[1], info[2]),
      n_genes = info[3],
      n_cells = info[5],
      n_cell_types = info[7],
      time = as.numeric(unlist(res$time), units = "secs"),
      memory = res$mem_alloc
    )
  }) %>% do.call("bind_rows", .) %>%
    dplyr::mutate(n_cells = as.numeric(n_cells), n_genes = as.numeric(n_genes))

  results$model_name %>% unique()
  results %>%
    dplyr::mutate(model_name = ifelse(model_name == "cpu devil", "devil", model_name)) %>%
    dplyr::mutate(model_name = ifelse(model_name == "cpu glmGamPoi", "glmGamPoi", model_name)) %>%
    dplyr::mutate(model_name = ifelse(model_name == "gpu devil", "devil (GPU)", model_name))
}

time_comparison = function(results, ratio = FALSE) {
  if (!isFALSE(ratio)) {
    results %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
      dplyr::summarise(y = mean(time)) %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types) %>%
      dplyr::mutate(y = y / y[model_name == ratio]) %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name)) +
      geom_bar(position = "dodge", stat = "identity") +
      #scale_y_continuous(transform = "log10") +
      #geom_point() +
      #geom_line() +
      #ggh4x::facet_nested(~"N cells"+n_cells, scales = "free", space = "free_y") +
      #facet_grid(~n_cells, scales = "free_y", shrink = T) +
      facet_wrap(~paste0(n_cells, " cells"), scales = "free_y", shrink = T) +
      theme_bw() +
      labs(x = "N genes", y = "Time ratio (s)", col = "Model", fill="Model") +
      scale_fill_manual(values = method_colors) +
      scale_color_manual(values = method_colors)
  } else {
    results %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
      dplyr::summarise(y = mean(time), sd =sd(time)) %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, col = model_name, fill=model_name)) +
      geom_bar(position = "dodge", stat = "identity") +
      #scale_y_continuous(transform = "log10") +
      #geom_point() +
      #geom_line() +
      #ggh4x::facet_nested(~"N cells"+n_cells, scales = "free", space = "free_y") +
      #facet_grid(~n_cells, scales = "free_y", shrink = T) +
      facet_wrap(~paste0(n_cells, " cells"), scales = "free_y", shrink = T) +
      theme_bw() +
      labs(x = "N genes", y = "Time (s)", col = "Model", fill="Model") +
      scale_fill_manual(values = method_colors) +
      scale_color_manual(values = method_colors)
  }
}

memory_comparison = function(results, ratio=FALSE) {
  if (!isFALSE(ratio)) {
    results %>%
      dplyr::mutate(memory = as.numeric(memory) * 1e-9) %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
      dplyr::summarise(y = mean(memory)) %>%
      dplyr::group_by(n_genes, n_cells, n_cell_types) %>%
      dplyr::mutate(y = y / y[model_name == ratio]) %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = y, fill = model_name)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      facet_wrap(~paste0(n_cells, " cells"), scales = "free_y", shrink = T) +
      theme_bw() +
      labs(x = "N genes", y = "Memory ratio (GB)", col = "Model", fill="Model") +
      scale_fill_manual(values = method_colors) +
      scale_color_manual(values = method_colors)
  } else {
    results %>%
      dplyr::mutate(memory = as.numeric(memory) * 1e-9) %>%
      ggplot(mapping = aes(x = as.factor(n_genes), y = memory, fill = model_name)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() +
      facet_wrap(~paste0(n_cells, " cells"), scales = "free_y", shrink = T) +
      theme_bw() +
      labs(x = "N genes", y = "Memory used (GB)", col = "Model", fill="Model") +
      scale_fill_manual(values = method_colors) +
      scale_color_manual(values = method_colors)
  }
}

time_per_gene_comparison = function(results) {

  results %>%
    dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(time / n_genes)) %>%
    ggplot(mapping = aes(x = as.factor(n_genes), y = y, fill = model_name)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    facet_wrap(~paste0(n_cells, " cells"), scales = "free_y", shrink = T) +
    theme_bw() +
    labs(x = "N genes", y = "Time per gene (s)", col = "Model", fill="Model") +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors)
}

plot_time_and_memory_comparison = function(results) {
  p1 <- results %>%
    dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(time), sd =sd(time)) %>%
    ggplot(mapping = aes(x = n_cells, y = y, ymin=y-sd, ymax=y+sd, col = model_name)) +
    #geom_bar(position = "dodge", stat = "identity") +
    geom_pointrange() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n_genes, scales = "free_x") +
    theme_bw() +
    labs(x = "N cells", y = "Time (s)", col = "Model")
  p1

  results %>%
    dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(time), sd =sd(time)) %>%
    ggplot(mapping = aes(x = n_genes, y = y / n_genes, col = model_name)) +
    #geom_bar(position = "dodge", stat = "identity") +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n_cells, scales = "free_x") +
    theme_bw() +
    labs(x = "N cells", y = "Time (s)", col = "Model")

  p1

  p2 <- results %>%
    dplyr::mutate(memory = as.numeric(memory) * 1e-9) %>%
    ggplot(mapping = aes(x = as.factor(n_cells), y = memory, fill = model_name)) +
    geom_bar(stat = "identity", position = "dodge", col = "black") +
    ggh4x::facet_nested(~"N genes"+n_genes, scales = "free_x") +
    theme_bw() +
    labs(x = "N cells", y = "Memory (GB)", fill = "Model")
  p2

  p3 <- results %>%
    dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
    dplyr::summarise(y = mean(time / n_genes)) %>%
    #dplyr::group_by(n_genes, n_cells, n_cell_types) %>%
    #dplyr::mutate(ratio_time = time / time[model_name == "devil (GPU)"]) %>%
    #dplyr::group_by(n_genes, n_cells, n_cell_types, model_name) %>%
    #dplyr::summarise(y = mean(ratio_time), sd =sd(ratio_time)) %>%
    ggplot(mapping = aes(x = n_cells, y = y, col = model_name)) +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n_genes, scales = "free_x") +
    theme_bw() +
    labs(x = "N cells", y = "Time per gene (s)", col = "Model")

  list(time=p1, memory=p2, ratio_time=p3)
}


plot_correlations <- function(fits_folder) {
  fits <- list.files(fits_folder, full.names = T)

  devil.res <- readRDS(fits[grepl("/cpu_devil_", fits)])
  gpu.devil.res <- readRDS(fits[grepl("/gpu_devil_", fits)])
  glm.res <- readRDS(fits[grepl("/cpu_glmGam", fits)])

  p1 <- dplyr::bind_rows(
    dplyr::tibble(
      x = gpu.devil.res$lfc,
      x_name = "gpu devil",
      y = devil.res$lfc,
      y_name = "cpu devil"
    ),
    dplyr::tibble(
      x = gpu.devil.res$lfc,
      x_name = "gpu devil",
      y = glm.res$lfc,
      y_name = "cpu glmGamPoi"
    )
  ) %>%
    dplyr::group_by(x_name, y_name) %>%
    dplyr::filter(abs(x) < 10 & abs(y) < 10) %>%
    ggplot(mapping = aes(x=x, y=y)) +
    geom_point() +
    ggpubr::stat_cor() +
    facet_grid(x_name~y_name, switch = "both") +
    theme_bw() +
    labs(x = bquote(LFC[1]), y=bquote(LFC[2])) +
    ggtitle("Log fold change correlation")
  p1

  p2 <- dplyr::bind_rows(
    dplyr::tibble(
      x = gpu.devil.res$theta,
      x_name = "gpu devil",
      y = devil.res$theta,
      y_name = "cpu devil"
    ),
    dplyr::tibble(
      x = gpu.devil.res$theta,
      x_name = "gpu devil",
      y = glm.res$theta,
      y_name = "cpu glmGamPoi"
    )
  ) %>%
    dplyr::group_by(x_name, y_name) %>%
    ggplot(mapping = aes(x=x, y=y)) +
    geom_point() +
    ggpubr::stat_cor() +
    facet_grid(x_name~y_name, switch = "both") +
    theme_bw() +
    labs(x = bquote(theta[1]), y=bquote(theta[2])) +
    ggtitle("Overdisperions correlation")
  p2

  list(lfc=p1, theta=p2)
}

plot_upset <- function(fits_folder, lfc_cut, pval_cut) {
  fits <- list.files(fits_folder, full.names = T)

  l = list()
  for (f in fits) {
    info = unlist(strsplit(unlist(strsplit(f, "fits//"))[2], "_"))
    name = paste(info[1], info[2])
    de_genes = readRDS(f) %>%
      dplyr::mutate(name = paste0("Gene ", row_number())) %>%
      dplyr::filter(abs(lfc) >= lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::pull(name)
    l[[name]] = de_genes
  }

  UpSetR::upset(UpSetR::fromList(l), order.by = "freq")
}
