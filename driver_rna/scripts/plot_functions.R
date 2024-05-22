###--------------------------------------------------------------###
### UMAP plot ###
###--------------------------------------------------------------###

library(Seurat)
library(ggplot2)
library(ggmin)

#saving plot cluster 
options(bitmapType='cairo')
png(file="umap_blood.png", width = 1000, height = 1000)

clustering1 <- DimPlot(seurat_rna, 
                      dims = c(1, 2),
                      reduction = "umap",
                      group.by = 'celltype.l2',
                      repel = TRUE, 
                      label = TRUE)
p1 <- clustering1 + ggmin::theme_powerpoint() +
  labs(title = "Custom clusters(muscle dataset, ~290 000 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

clustering2 <- DimPlot(seurat_retina, 
                       dims = c(1, 2),
                       reduction = "umap",
                       group.by = 'majorclass',
                       repel = TRUE, 
                       label = TRUE
                       
)

p2 <- clustering2 + ggmin::theme_powerpoint() +
  labs(title = "Custom labels (cell types)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

p1+p2

dev.off()


###-----------------------------------------------------------------###
### Volcano Plot ###
###-----------------------------------------------------------------###

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

plot1 <- p1 + ggplot2::labs(title="Macula vs Peripheral retina, RPCs (66 000 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))


### Visualize fgsea enrichment results ###
# making ggdotchart
library(ggpubr)

res_gse$log_padjust <- -log10(res_gse$p.adjust)

plot1 <- ggdotchart(res_gse, x = "Description", y = "log_padjust",
                    color = "gene_clusters",
                    palette = c("#00AFBB", "#FC4E07"),
                    sorting = "descending",
                    rotate = TRUE,
                    group = "gene_clusters",
                    dot.size = 6,
                    add = "segments",
                    title = "Late RPCs + Muller cells vs Early RPCs",
                    xlab = "Pathways",
                    ylab = "-log10(padjust)",
                    ggtheme = theme_pubr()
)
plot1 + plot2


### Heatmap of the most significant genes ###

gene_list1 <- c("FGF19", "CYP1B1", "CYP26A1", "DIO2", "ANXA2", "FRZB", "CRYAB", "HES1", "PTGDS", 
                "GPX3", "FOXG1", "TBX20")

gene_list2 <- c("FGF19", "HAS2", "ZNF676", "FOS", "NFIB", "PID1", "NETO1", "PTGDS", "ATOH7", 
                "PTF1A", "OTX2", "LUM")

heatmap1 <- Seurat::DoHeatmap(subset(sc_retina, downsample = 10000), 
                              features = gene_list1, 
                              group.by = "tissue", 
                              size = 3,
                              group.label.rot = TRUE)

heatmap1 <- heatmap1 + scale_fill_gradientn(colors = c("darkblue", "white", "yellow")) +
  labs(title = "Macula lutea vs Peripheral retina") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

heatmap2 <- Seurat::DoHeatmap(subset(sc_retina, downsample = 10000), 
                              features = gene_list2, 
                              group.by = "dev_stage", 
                              size = 3)

heatmap2 <- heatmap2 + scale_fill_gradientn(colors = c("darkblue", "white", "yellow")) +
  labs(title = "Late RPCs vs Early RPCs") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))


### Using dittoSeq ###
#BiocManager::install("dittoSeq")
library(dittoSeq)



