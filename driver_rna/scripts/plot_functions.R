#!/usr/bin/env Rscript


### UMAP plot ###
library(Seurat)
library(ggplot2)
library(ggmin)
#saving plot cluster
options(bitmapType='cairo')
png(file="plots/umap.png", width = 900, height = 600)

clustering <- DimPlot(sc_retina_atac, 
                      dims = c(1, 2),
                      reduction = "umap",
                      group.by = 'tissue',
                      repel = TRUE, 
                      label = FALSE
)
clustering + ggmin::theme_powerpoint() +
  labs(title = "Human fetal retina snATAC tissue specific clustering (223 288 cells)") +
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

res <- stat_test_retina %>% remove_rownames %>% column_to_rownames(var="name")

p1 <- EnhancedVolcano::EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'lfc',
                y = 'adj_pval',
                selectLab = c("RGS16", "ATOH7", "GADD45G", "HES1", "OTX2", "SFRP2", "HKDC1", "CNGB1", "NFIC"),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 3.0,
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
plot <- p1 + ggplot2::labs(title="Differentiated neurons vs RPCs, (145 000 cells)") +
                     theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot         

res <- res[!rownames(res) %in% c(), ]

### Visualize fgsea enrichment results ###
ggplot(data=screencast, aes(x=Reason,y=Percentage,fill=factor(Type))) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  ggtitle("Strategies for Using Homework Solution and Mini-Lecture Screencasts")



