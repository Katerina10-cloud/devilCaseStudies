
### UMAP plot ###

library(Seurat)
library(ggplot2)
library(ggmin)
library(tidyverse)

#saving plot cluster
options(bitmapType='cairo')
png(file="plots/umap_atac.png", width = 1800, height = 600)

clustering1 <- DimPlot(sc_retina, 
                      dims = c(1, 2),
                      reduction = "umap",
                      group.by = 'cell_type',
                      repel = TRUE, 
                      label = FALSE
)
p1 <- clustering1 + ggmin::theme_powerpoint() +
  labs(title = "HFR snATAC cell specific (223 288 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

clustering2 <- DimPlot(sc_retina, 
                      dims = c(1, 2),
                      reduction = "umap",
                      group.by = 'tissue',
                      repel = TRUE, 
                      label = FALSE
)
p2 <- clustering2 + ggmin::theme_powerpoint() +
  labs(title = "HFR snATAC tissue specific (223 288 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p1 + p2
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
library(patchwork)

res <- subset(stat_test_retina, lfc >= -5)
res <- res %>% remove_rownames %>% column_to_rownames(var="name")

# Tissue specific markers RPC: FGF19", "CYP1B1", "CYP26A1", "DIO2", "CDKN1A", "ANXA2", "FRZB", "CRYAB", "HES1"
# Diff. neurons vs RPCs : "HKDC1", "HES1", "SFRP2", "ATOH7", "RGS1G", "GADD45G", "NFIC", "OTX2", "CNGB1"
# RGCs vs RPCs : "HES1", "SFRP2", "FRZB", "ATOH7", "NFIA", "PTF1A", "MYC", "ISL1", "POU4F1", "POU6F2", "EBF3"
# Dev.stage: "FGF19", "HAS2", "ZNF676", "FOS", "NFIB", "PID1", "NETO1", "PTGDS", "ATOH7", "PTF1A", "OTX2"

p1 <- EnhancedVolcano::EnhancedVolcano(res,
                lab = rownames(res),
                x = 'lfc',
                y = 'adj_pval',
                selectLab = c("FGF19", "HAS2", "ZNF676", "FOS", "NFIB", "PID1", "NETO1", "PTGDS", "ATOH7", "PTF1A", "OTX2"),
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
                subtitle = "",
                gridlines.major = F, gridlines.minor = F,
                border = 'full', borderWidth = 0.5, borderColour = 'black',
                titleLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

plot1 <- p1 + ggplot2::labs(title="Macula vs Peripheral retina, RPCs (66 000 cells)") +
                     theme(plot.title=element_text(hjust=0.5, vjust=0.5))

plot2 <- p1 + ggplot2::labs(title="Differentiated neurons vs RPCs (145 000 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

plot3 <- p1 + ggplot2::labs(title="RGCs vs RPCs (107 000 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

plot4 <- p1 + ggplot2::labs(title="Early RPCs vs Late RPCs + Muller cells (75 000 cells)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

#wrap_plots(plot1, plot2, plot3, plot4)


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
library("ComplexHeatmap")
library("RcolorBrewer")
library("circlize")
library("grid")

# Extract DEG #

# Macula vs Peripheral retina
res_deg1 <- res_deg %>%
  filter(name %in% c("CYP1B1", "FGF19", "HES1", "FOXG1", "TBX20", "DIO2", "ANXA2", "FRZB", "CRYAB", "CYP26A1",
                     "CDKN1A", "PTGDS", "GPX3", "APOE"))

# RGCs vs RPCs
res_deg2 <- res_deg %>%
  filter(name %in% c("HES1", "SFRP2", "NFIA", "FRZB", "PTFA1", "ATOH7", "MYC", "POU4F1", "ISL2", "EBF3",
                     "POUGF2"))

# Differentiated neurons vs RPCs
res_deg3 <- res_deg %>%
  filter(name %in% c("HKDC1", "SFRP2", "HES1", "ATOH7", "NFIC", "GADD456", "OTX2", "CNGB1", "PCA3", "MYO3B"))

# Early RPCs vs Late RPCs
res_deg4 <- res_deg %>%
  filter(name %in% c("HAS2", "FGF19", "ZNF676", "ATOH7", "PTF1A", "OTX2", "FOS", "NFIB", "PID1", 
                     "NETO1", "PTGDS", "FGF19"))

ann_col_info <- as.data.frame(metadata)
anno_info_colors <- list(cell_type = c(retinal progenitor cell = "lightblue",
                                       retinal ganglion cell = "red"))

plot1 <- pheatmap(rna_exp,
                  cluster_rows = FALSE,
                  cluster_cols = TRUE,
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  scale = "row",
                  color =colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  annotation_col = ann_col_info,
                  annotation_colors = anno_info_colors,
                  main = "DGE ")










