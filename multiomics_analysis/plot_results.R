# Plot results #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis/results/")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

load("~/metadata_atac_umap.Rdata")
load("~/metadata_rna_umap.Rdata")

### UMAP ###
metadata_rna <- metadata_rna[ (metadata_rna$tech_id %in% c("1")), ]
umap_rna <- metadata_rna %>% dplyr::select(umap_1, umap_2, cell_type)
umap_rna$cell_type <- as.factor(umap_rna$cell_type)

p_umap_rna <- ggplot(umap_rna, aes(x = umap_1, y = umap_2, color = cell_type))+
  geom_point(size=0.2)+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snRNA (~ 212 000 nuclei)",
       color = "cell_type") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_classic() +
  scale_fill_manual("cell_type")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)
  
p_umap_rna <- Seurat::LabelClusters(plot = p_umap_rna, id = 'cell_type', color="black")+
  theme(legend.position="none")


umap_atac <- metadata_atac %>% dplyr::select(UMAP_1, UMAP_2, cell_type)
umap_atac$cell_type <- as.factor(umap_atac$cell_type)

p_umap_atac <-ggplot(umap_atac, aes(x = UMAP_1, y = UMAP_2, color = cell_type))+
  geom_point(size=0.2)+
  labs(x = "UMAP_1",
       y = "UMAP_2",
       subtitle = "snATAC (~ 90 000 nuclei)")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))+
  theme_classic() +
  scale_fill_manual("cell type")+
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE)

p_umap_atac <- Seurat::LabelClusters(plot = p_umap_atac, id = 'cell_type', color="black")+
  theme(legend.position="none")

p_umap_rna + p_umap_atac


### Volcano plot ###
  
# Top cell type & muscle aging specific genes #
top_genes <- filter(res, name %in% c("PPARA","PER2","MYH1","MYH2","MYH4","PDE7B", 
                     "TNNT2","ID1","SAA2-SAA4","JUN","JUND","FOS","EGR1")) 

data_path <- "MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(data_path)
p <- EnhancedVolcano::EnhancedVolcano(rna_devil,
                                       lab = rna_devil$name,
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

plot1 <- p + ggplot2::labs(title="Devil: snRNA") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot1






