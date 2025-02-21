
require(ggupset)

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

plot_time_and_memory_comparison = function(results, n_extrapolation = 3, ncols=2) {
  d <- results %>%
    dplyr::group_by(model_name, n_genes, n_cells) %>%
    dplyr::summarise(memory = mean(memory), time = mean(time)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(size = n_genes * n_cells)

  # Extraolate unkwonw sizes
  sizes = unique(d$size)
  models <- unique(d$model_name)

  d <- lapply(models, function(m) {
    dd <- d %>% dplyr::filter(model_name == m) %>% dplyr::arrange(-size)

    x = c(log10(dd$size[1:n_extrapolation]))
    y_time = c(log10(dd$time[1:n_extrapolation]))
    y_memory = c(log10(dd$memory[1:n_extrapolation]))

    time_lm <- lm(y_time ~ x)
    memory_lm <- lm(y_memory ~ x)

    ddd = dd %>%
      ungroup() %>%
      dplyr::filter(size %in% sizes) %>%
      dplyr::select(model_name, time, memory, size) %>%
      dplyr::mutate(is_extrapolated = "FALSE")

    sizes_to_extrapolate <- sizes[!(sizes %in% ddd$size)]
    if (length(sizes_to_extrapolate) > 0) {
      s <- sizes_to_extrapolate[1]

      dd_extr <- lapply(sizes_to_extrapolate, function(s) {
        time_pred = 10^predict(time_lm, newdata = data.frame(x = c(log10(s))))
        mem_pred = 10^predict(memory_lm, newdata = data.frame(x = c(log10(s))))
        dplyr::tibble(model_name = m, time = time_pred, memory = mem_pred, size = s, is_extrapolated = "TRUE")
      }) %>% do.call("bind_rows", .)

      ddd <- dplyr::bind_rows(ddd, dd_extr)
    }

    if (any(ddd$is_extrapolated == "TRUE")) {
      fake_row = ddd %>% dplyr::filter(is_extrapolated == "FALSE") %>% dplyr::filter(size == max(size))
      fake_row$is_extrapolated <- "TRUE"
      ddd <- dplyr::bind_rows(ddd, fake_row)
    }

    ddd
  }) %>% do.call("bind_rows", .)

  # Compute Ratios
  d <- d %>%
    dplyr::group_by(size) %>%
    dplyr::mutate(time_ratio = time / time[model_name == "devil (GPU)"]) %>%
    dplyr::mutate(memory_ratio = memory / memory[model_name == "devil (GPU)"]) %>%
    dplyr::mutate(is_extrapolated = ifelse(is_extrapolated == "FALSE", "Observed", "Extrapolated")) %>%
    dplyr::mutate(is_extrapolated = factor(is_extrapolated, levels = c("Observed", "Extrapolated")))


  p = d %>%
    `colnames<-`(c("model_name", "Time (s)", "Memory(GB)", "size", "is_extrapolated", "Speed ratio", "Memory consumption ratio")) %>%
    tidyr::pivot_longer(!c(model_name, is_extrapolated, size)) %>%
    dplyr::mutate(name = factor(name, levels = c("Time (s)", "Memory(GB)", "Speed ratio", "Memory consumption ratio"))) %>%
    ggplot(mapping = aes(x = size, y= value, col = model_name, linetype = is_extrapolated)) +
    geom_point(aes(shape = is_extrapolated), size = 3) +
    geom_line() +
    scale_x_continuous(transform = "log10") +
    scale_y_continuous(transform = "log10") +
    labs(x = "Dataset size", y = "", col = "Model", shape="Measurement", linetype="Measurement") +
    scale_size_continuous(guide = "none") +
    theme_bw() +
    scale_color_manual(values = method_colors) +
    facet_wrap(name~., scales = "free", ncol = ncols, strip.position = "left")

  return(p)

  # p_time = d %>%
  #   ggplot(mapping = aes(x = size, y= time, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Time (s)", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # p_time_ratio = d %>%
  #   ggplot(mapping = aes(x = size, y= time_ratio, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Time ratio", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # p_mem <- d %>%
  #   ggplot(mapping = aes(x = size, y= memory * 1e-9, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Memory (GB)", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # p_mem_ratio <- d %>%
  #   ggplot(mapping = aes(x=size, y= memory_ratio, col = model_name, linetype = is_extrapolated)) +
  #   geom_point(aes(shape = is_extrapolated), size = 3) +
  #   geom_line() +
  #   scale_x_continuous(transform = "log10") +
  #   scale_y_continuous(transform = "log10") +
  #   labs(x = "Dataset size", y = "Memory ratio", col = "Model", shape="Measurement", linetype="Measurement") +
  #   scale_size_continuous(guide = "none") +
  #   theme_bw() +
  #   scale_color_manual(values = method_colors)
  #
  # list(
  #   p_time=p_time,
  #   p_time_ratio=p_time_ratio,
  #   p_mem=p_mem,
  #   p_mem_ratio=p_mem_ratio
  # )
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

fits_folder <- "results/baronPancreas/fits/"
pval_cut <- .05
lfc_cut <- 1
plot_upset <- function(fits_folder, lfc_cut, pval_cut) {
  fits <- list.files(fits_folder, full.names = T)

  l = list()
  fits = fits[grepl("cpu_devil_", fits) | grepl("gpu_devil_", fits) | grepl("glmGamPoi", fits)]
  res <- lapply(fits, function(f) {
    info = unlist(strsplit(unlist(strsplit(f, "fits//"))[2], "_"))
    name = paste(info[1], info[2])
    de_genes = readRDS(f) %>%
      dplyr::mutate(name = paste0("Gene ", row_number())) %>%
      dplyr::filter(abs(lfc) >= lfc_cut, adj_pval <= pval_cut) %>%
      dplyr::pull(name)
    dplyr::tibble(algorithm = name, de_genes = de_genes)
  }) %>% do.call("bind_rows", .)

  res %>%
    group_by(de_genes) %>%
    summarize(algorithm = list(algorithm)) %>%
    ggplot(aes(x = algorithm)) +
    geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
    scale_y_continuous(breaks = NULL, lim = c(0, nrow(res) / length(fits)), name = "") +
    scale_x_upset() +
    theme_bw()
}

plot_volc = function (
    devil.res,
    lfc_cut = 1,
    pval_cut = 0.05,
    labels = TRUE,
    colors = c("gray", "forestgreen", "steelblue", "indianred"),
    color_alpha = 0.7,
    point_size = 1,
    center = TRUE,
    title = "Volcano plot") {

  if (sum(is.na(devil.res))) {
    message("Warning: some of the reults are unrealiable (i.e. contains NaN)\n Those genes will not be displayed")
    devil.res <- stats::na.omit(devil.res)
  }

  d <- devil.res %>%
    dplyr::mutate(pval_filter = .data$adj_pval <= pval_cut, lfc_filter = abs(.data$lfc) >= lfc_cut) %>%
    dplyr::mutate(class = dplyr::if_else(
      .data$pval_filter & .data$lfc_filter,
      "p-value and lfc",
      dplyr::if_else(
        .data$pval_filter,
        "p-value",
        dplyr::if_else(
          .data$lfc_filter,
          "lfc",
          "non-significant"
          )
        )
      )
      ) %>% dplyr::mutate(class = factor(class, levels = c("non-significant", "lfc", "p-value", "p-value and lfc"))) %>%
    dplyr::mutate(label = dplyr::if_else(class == "p-value and lfc", .data$name, NA))
  d

  if (sum(d$adj_pval == 0) > 0) {
    message(paste0(sum(d$adj_pval == 0), " genes have adjusted p-value equal to 0, will be set to the minimum non-zero value"))
    d$adj_pval[d$adj_pval == 0] <- min(d$adj_pval[d$adj_pval != 0])
  }

  p <- d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$lfc, y = -log10(.data$adj_pval), col = .data$class, label = .data$label)) +
    ggplot2::geom_point(alpha = color_alpha, size = point_size) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col = "") + ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed") + ggplot2::geom_hline(yintercept = -log10(pval_cut),
                                                                   linetype = "dashed") + ggplot2::ggtitle(title) + ggplot2::theme(legend.position = "bottom")
  if (center) {
    p <- p + ggplot2::xlim(c(-max(abs(d$lfc)), max(abs(d$lfc))))
  }

  if (labels) {
    p = p + ggrepel::geom_label_repel(max.overlaps = Inf, col = "black")
  }
  p
}
