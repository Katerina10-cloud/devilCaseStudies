#!/usr/bin/env Rscript


### UMAP plot ###
library(Seurat)
library(ggplot2)
library(ggmin)
#saving plot cluster
options(bitmapType='cairo')
png(file="plots/umap.png", width = 900, height = 600)

clustering <- DimPlot(seurat, 
                      dims = c(1, 2),
                      reduction = "umap",
                      group.by = 'cell_type',
                      repel = TRUE, 
                      label = FALSE
)
clustering + ggmin::theme_powerpoint() +
  labs(title = "Human Fetal retina, 223 288 (Human Cell Atlas)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
dev.off()

### Volcano Plot ###

library(ggplot2)
library(gridExtra)
library(tidyverse)

stat_test_retina$diffexpressed <- "NO"
stat_test_retina$diffexpressed[stat_test_retina$lfc >= 1 & stat_test_retina$adj_pval < 0.05] <- "UP"
stat_test_retina$diffexpressed[stat_test_retina$lfc <= -1 & stat_test_retina$adj_pval < 0.05] <- "DOWN"

p1 <- ggplot(data = stat_test_retina, aes(x = lfc, y = -log10(adj_pval), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-8, 8, 1))+
  labs(title="Devil:scRNA, multipatient (80 000 cells)",x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p1
