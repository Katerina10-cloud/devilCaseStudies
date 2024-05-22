# Plot results #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel")
sapply(pkgs, require, character.only = TRUE)

set.seed(12345)

# UMAP #

metadata_rna <- metadata_rna[ (metadata_rna$tech_id %in% c("1")), ]
umap_rna <- metadata_rna %>% dplyr::select(umap_1, umap_2, cell_type)
umap_rna$cell_type <- as.factor(umap_rna$cell_type)

#Seurat::DiscretePalette()

p_umap_rna <- ggplot(umap_rna, aes(x = umap_1, y = umap_2, 
                                       color = cell_type))+
  geom_point()+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snRNA (~ 212 000 cells)",
       color = "cell type") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_minimal() +
  scale_fill_manual("cell_type")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)
p_umap_rna <- Seurat::LabelClusters(plot = p_umap_rna, id = 'cell_type', color="black")


umap_atac <- metadata_atac %>% dplyr::select(UMAP_1, UMAP_2, cell_type)
umap_atac$cell_type <- as.factor(umap_atac$cell_type)

p_umap_atac <-ggplot(umap_atac, aes(x = UMAP_1, y = UMAP_2, 
                                            color = cell_type))+
  geom_point()+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snATAC (~ 90 000 cells)")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_minimal() +
  scale_fill_manual("cell type")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)

p_umap_atac <- Seurat::LabelClusters(plot = p_umap_atac, id = 'cell_type', color="black")

p_umap_rna + p_umap_atac


# Vienn diagram #

# Select DE genes #
sum(res_atac_rna$adj_pval_snATAC < 0.05 & res_atac_rna$lfc_snATAC >= 1.0, na.rm=TRUE) #up_regulated
sum(res_atac_rna$adj_pval_snATAC < 0.05 & res_atac_rna$lfc_snATAC <= -1.0, na.rm=TRUE) #down-reg

sum(res_atac_rna$adj_pval_snRNA < 0.05 & res_atac_rna$lfc_snRNA >= 1.0, na.rm=TRUE) #up_regulated
sum(res_atac_rna$adj_pval_snRNA < 0.05 & res_atac_rna$lfc_snRNA <= -1.0, na.rm=TRUE) #down-reg

overlap <- draw.pairwise.venn(area1=4822, area2=1024,cross.area=401, 
                              category=c("snATAC","snRNA"),fill=c("Red","Blue"))


# Correlation plot #

res_atac_rna <- res_atac_rna[-c(101,290),]

cor_coefs <- cor.test(res_atac_rna$lfc_snRNA, res_atac_rna$lfc_snATAC)
corr_plot <- ggplot(data = res_atac_rna, aes(x = lfc_snRNA, y = lfc_snATAC)) + 
  geom_point(shape = 21, fill = 'black', size = 2) +
  geom_smooth(method=lm , color="red", fill="black", se=TRUE) +
  labs(title = "Correlation of significantly DAG (snATAC-seq) and GE (snRNA-seq)")+
  xlab("LogFC2 snRNA") +
  ylab ("LogFC2 snATAC") +
  theme_minimal()+
  annotate("text", x = 10, y = 4, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 10, y = 3.5, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))
corr_plot

corr_plot <- ggplot(mapping = aes(x = res_atac_rna$lfc_snRNA, y = res_atac_rna$lfc_snATAC)) +
  geom_point(shape = 21, fill = 'black', size = 2) +
  labs(title = "Correlation of significantly DAG (snATAC-seq) and GE (snRNA-seq)")+
  xlab("snRNA log2FC") +
  ylab ("snATAC log2FC") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  sm_statCorr()+
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_minimal()
corr_plot


# Volcano plot #

# top genes
top_genes <- filter(res, name %in% c("PPARA","PER2","MYH1","MYH2","MYH4","PDE7B", 
                     "TNNT2","ID1","SAA2-SAA4","JUN","JUND","FOS","EGR1")) 

p2 <- EnhancedVolcano::EnhancedVolcano(res_atac_nodup,
                                       lab = res_atac_nodup$geneID,
                                       x = 'lfc',
                                       y = 'adj_pval',
                                       selectLab = c("PPARA","PER2","MYH1","MYH2","MYH4","PDE7B", 
                                                     "TNNT2","ID1","SAA2-SAA4","JUN","JUND","FOS","EGR1"),
                                       xlab = bquote(~Log[2]~ 'fold change'),
                                       pCutoff = 0.05,
                                       FCcutoff = 1.0,
                                       pointSize = 1.0,
                                       labSize = 2.0,
                                       labCol = 'black',
                                       labFace = 'bold',
                                       boxedLabels = TRUE,
                                       col=c('gray', 'gray', 'gray', 'purple'),
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

plot1 <- p1 + ggplot2::labs(title="snRNA (old cohort myonuclei sv young cohort)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

plot2 <- p2 + ggplot2::labs(title="snATAC (old cohort myonuclei sv young cohort)") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

plot1 + plot2
grid.arrange(plot1, plot2, nrow = 1)

### Visualize Gsea enrichment results ###

res_gse$log_padjust <- -log10(res_gse$p.adjust)

plot1 <- ggdotchart(res_gse, x = "Description", y = "log_padjust",
                    color = "gene_clusters",
                    palette = c("blue", "#FC4E07"),
                    sorting = "descending",
                    rotate = TRUE,
                    group = "gene_clusters",
                    dot.size = 6,
                    add = "segments",
                    title = "Enrichment of key identified genes",
                    xlab = "Pathways",
                    ylab = "-log10(padjust)",
                    ggtheme = theme_pubr()
)
plot1
