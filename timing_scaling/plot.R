rm(list = ls())
require(tidyverse)
require(ggplot2)

t <- readRDS("timings.rds")
t <- lapply(1:nrow(t), function(i) {
  times <- t[i,]$time %>% unlist()
  lapply(1:length(times), function(j) {
    t[i,] %>% dplyr::mutate(timing = times[j], iter=j)
  }) %>% do.call('bind_rows', .)

}) %>% do.call('bind_rows', .)

t %>%
  ggplot(mapping = aes(x=n.cell, y=median, col=algo)) +
  geom_point() +
  geom_line()

t %>%
  ggplot(mapping = aes(x=n.cell, y=timing, col=algo)) +
  geom_point() +
  geom_smooth() +
  scale_y_continuous(trans = 'log10')
