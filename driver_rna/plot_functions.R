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
library(ggpubr)
library(tidyverse)

stat_test_res$diffexpressed <- "NO"
stat_test_res$diffexpressed[stat_test_res$lfc >= 1 & stat_test_res$adj_pval < 0.05] <- "UP"
stat_test_res$diffexpressed[stat_test_res$lfc <= -1 & stat_test_res$adj_pval < 0.05] <- "DOWN"

p1 <- ggplot(data = stat_test_res, aes(x = lfc, y = -log10(adj_pval), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-8, 8, 2))+
  labs(title="Devil:scRNA, single patients (7400 cells)",x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p1

### Enhanced volcano ###

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(ggrepel)

res <- stat_test_res %>% remove_rownames %>% column_to_rownames(var="name")

p1 <- EnhancedVolcano(res,
                lab = rownames(res),
                x = 'lfc',
                y = 'adj_pval',
                selectLab = c("PRUNE2", "MTUS2", "KLHL1", "CDH18", 
                              "FGFR2", "FGF19", "CASP7"),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                col=c('black', 'black', 'black', 'purple'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                #title = 'Devil: single patient (7400 cells)',
                subtitle = "",
                #legend=c("NS", "Log2 FC", "padj", "padj & Log2 FC"),
                gridlines.major = F, gridlines.minor = F,
                border = 'full', borderWidth = 0.5, borderColour = 'black',
                titleLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
plot <- p1 + ggplot2::labs(title="Human neural retina DGE: single patient (7400 cells)") +
                     theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot                     
