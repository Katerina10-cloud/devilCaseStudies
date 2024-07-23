
method_colors = c(
  "glmGamPoi (Pb)" = "#A22E29",
  "edgeR" = "#7D629E",
  "limma" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578",
  "Nebula" =  "#E4A6A7",
  "Devil" = "#ABD5D0",
  "DevilCL" = "#31482F"
)

plot_null <- function(n_samples, ct.index, author, is.pb, algos = c("glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil", "DevilCL")) {
  if (is.pb) {
    head_foler = "nullpower/null_subject/"
  } else {
    head_foler = "nullpower/null_cell/"
  }

  d <- read.delim(paste0(head_foler, author ,".n.", n_samples, ".ct.",ct.index,".fc.0.5.csv"), sep = ",")
  colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil", "DevilCL")

  mask <- colnames(d) %in% algos
  mask[1] <- TRUE
  d <- d[,mask]

  cols <- colnames(d)
  d <- lapply(2:ncol(d), function(c) {
    values = d[,c] %>% sort()
    x = seq(0,1,length = length(values))
    dplyr::tibble(x = x, observed_p_value = values, name = colnames(d)[c])
  }) %>% do.call("bind_rows", .)

  d$name %>% unique()

  p_null <- d %>%
    #dplyr::filter(grepl("devil", name) | grepl("Nebula", name)) %>%
    ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
    geom_point(size = .8) +
    #geom_line() +
    geom_abline(slope = 1, intercept = 0, col = "black", linetype="dashed") +
    ggtitle(paste0("Cell type ", ct.index, " - ", n_samples, " patients"), subtitle = paste0("Patient hierarchy ", is.pb)) +
    theme_bw() +
    labs(x = "Expected p-value", y="Observed p-value", col="Algorithm") +
    #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
    ggsci::scale_color_locuszoom() +
    theme(legend.position = "bottom")

  p_null
}


plot_pow <- function(n_samples, ct.index, author, is.pb, algos = c("glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil", "DevilCL")) {
  if (is.pb) {
    head_foler = "nullpower/pow_subject/"
  } else {
    head_foler = "nullpower/pow_cell/"
  }

  d <- read.delim(paste0(head_foler, author ,".n.", n_samples, ".ct.",ct.index,".fc.0.5.csv"), sep = ",")

  colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil", "DevilCL")

  mask <- colnames(d) %in% algos
  mask[1] <- TRUE

  d <- d[,mask]
  cols <- colnames(d)
  d <- lapply(2:ncol(d), function(c) {
    values = d[,c] %>% sort(decreasing = TRUE)
    values[values <= 1e-300] = 1e-300
    x = seq(0,1,length = length(values))
    dplyr::tibble(x = x, observed_p_value = -log10(values), name = colnames(d)[c])
  }) %>% do.call("bind_rows", .)

  d$name %>% unique()

  p_pow <- d %>%
    #dplyr::filter(grepl("devil", name) | grepl("Nebula", name)) %>%
    ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
    #geom_line() +
    geom_point(size = .8) +
    ggtitle(paste0("Cell type ", ct.index, " - ", n_samples, " patients"), subtitle = paste0("Patient hierarchy ", is.pb)) +
    theme_bw() +
    scale_y_continuous(trans = "log10") +
    labs(x = "Expected p-value", y="Observed p-value", col="Algorithm") +
    theme(legend.position = "bottom") +
    #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
    ggsci::scale_color_locuszoom()

  p_pow
}
