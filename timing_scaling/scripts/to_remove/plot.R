rm(list = ls())
library(tidyverse)
library(ggplot2)
library(patchwork)

# Beta ONLY ####
n_cond <- 1
for (n_cond in c(1,2,4)) {
  t <- readRDS(paste0("results/",n_cond,"_timing_beta_only.rds"))

  t <- lapply(1:nrow(t), function(i) {
    times <- t[i,]$time %>% unlist()
    lapply(1:length(times), function(j) {
      t[i,] %>% dplyr::mutate(timing = times[j], iter=j)
    }) %>% do.call('bind_rows', .)

  }) %>% do.call('bind_rows', .)

  t %>%
    ggplot(mapping = aes(x=n.cell, y=timing, col=algo)) +
    geom_point() +
    geom_line() +
    facet_wrap(~n.genes) +
    scale_y_continuous(trans = 'log10')

  p1 <- t %>%
    ggplot(mapping = aes(x=n.cell, y=timing, col=algo)) +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n.genes) +
    theme_bw() +
    scale_y_continuous(trans = 'log10') +
    labs(y = "Time(s)") +
    theme(legend.position = 'bottom') +
    ggtitle(paste0("Timing"), paste0("Dataset with ", max(t$n.cell), " cells and ", n_cond, " features"))

  p2 <- t %>%
    dplyr::group_by(iter, p.genes, p.cell) %>%
    dplyr::mutate(ratio_time = timing / timing[algo == 'devil_gpu']) %>%
    ggplot(mapping = aes(x=n.cell, y=ratio_time, col=algo)) +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n.genes) +
    theme_bw() +
    theme(legend.position = 'bottom')


  p <- p1 / p2
  ggsave(paste0("img/", n_cond, "_features_beta_only.pdf"), plot = p, dpi=600, width = 10, height = 7, units = 'in')
}

# beta AND overdisp ####
for (n_cond in c(1,2,4)) {
  t <- readRDS(paste0("results/",n_cond,"_timing.rds"))

  t <- lapply(1:nrow(t), function(i) {
    times <- t[i,]$time %>% unlist()
    lapply(1:length(times), function(j) {
      t[i,] %>% dplyr::mutate(timing = times[j], iter=j)
    }) %>% do.call('bind_rows', .)

  }) %>% do.call('bind_rows', .)

  t %>%
    ggplot(mapping = aes(x=n.cell, y=timing, col=algo)) +
    geom_point() +
    geom_line() +
    facet_wrap(~n.genes) +
    scale_y_continuous(trans = 'log10')

  p1 <- t %>%
    ggplot(mapping = aes(x=n.cell, y=timing, col=algo)) +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n.genes) +
    theme_bw() +
    scale_y_continuous(trans = 'log10') +
    labs(y = "Time(s)") +
    theme(legend.position = 'bottom') +
    ggtitle(paste0("Timing"), paste0("Dataset with ", max(t$n.cell), " cells and ", n_cond, " features"))

  p2 <- t %>%
    dplyr::group_by(iter, p.genes, p.cell) %>%
    dplyr::mutate(ratio_time = timing / timing[algo == 'devil_gpu']) %>%
    ggplot(mapping = aes(x=n.cell, y=ratio_time, col=algo)) +
    geom_point() +
    geom_line() +
    ggh4x::facet_nested(~"N genes"+n.genes) +
    theme_bw() +
    theme(legend.position = 'bottom')


  p <- p1 / p2
  ggsave(paste0("img/", n_cond, "_features.pdf"), plot = p, dpi=600, width = 10, height = 7, units = 'in')
}

