-------------------------------------------------------
# Plot results #
-------------------------------------------------------  

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis/results/devil/")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

set.seed(12345)

-------------------------------------------------------
### UMAP ###
-------------------------------------------------------
metadata_rna <- metadata_rna[ (metadata_rna$tech_id %in% c("1")), ]
umap_rna <- metadata_rna %>% dplyr::select(umap_1, umap_2, cell_type)
umap_rna$cell_type <- as.factor(umap_rna$cell_type)

#Seurat::DiscretePalette()

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


--------------------------------------------------------
# Vienn diagram #
--------------------------------------------------------  

devil <- list(snATAC=atac_devil$geneID, snRNA=rna_devil$geneID)
glm <- list(snATAC=atac_glm$geneID, snRNA=rna_glm$geneID)

p1 <- ggvenn(devil, c("snATAC", "snRNA"), show_percentage = FALSE,
             set_name_size = 4)
p1 <- p1 + ggplot2::labs(title="Devil") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

p2 <- ggvenn(glm, c("snATAC", "snRNA"), show_percentage = FALSE,
             set_name_size = 4)
p2 <- p2 + ggplot2::labs(title="glnGamPoi") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))

grid.arrange(p1, p2, nrow = 2)


----------------------------------------------------------
# Correlation plot #
----------------------------------------------------------
  
corr_plot2 <- ggplot(mapping = aes(x = overl_glm$lfc_snRNA, y = overl_glm$lfc_snATAC)) +
  geom_point(shape = 21, fill = 'black', size = 2) +
  labs(title = "glmGamPoi")+
  xlab("snRNA log2FC") +
  ylab ("snATAC log2FC") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  sm_statCorr()+
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()
corr_plot1+corr_plot2


---------------------------------------------------------
# Volcano plot #
---------------------------------------------------------  

# top genes #
top_genes <- filter(res, name %in% c("PPARA","PER2","MYH1","MYH2","MYH4","PDE7B", 
                     "TNNT2","ID1","SAA2-SAA4","JUN","JUND","FOS","EGR1")) 

p1 <- EnhancedVolcano::EnhancedVolcano(res_rna_neb,
                                       lab = res_rna_neb$name,
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

plot1 <- p1 + ggplot2::labs(title="Nebula: snRNA") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot1

plot2 <- p2 + ggplot2::labs(title="Devil: snATAC") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot2

plot3 <- p3 + ggplot2::labs(title="glmGamPoi: snRNA") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot3

plot4 <- p4 + ggplot2::labs(title="glmGamPoi: snATAC") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot4

grid.arrange(plot3, plot4, nrow = 1)


------------------------------------------------------------
### Gsea enrichment results ###
------------------------------------------------------------

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


---------------------------------------------------------
### P-values and lfc comparison ###
---------------------------------------------------------
  
devil <- devil[devil$geneID %in% glm$geneID,]
glm <- glm[glm$geneID %in% devil$geneID,]

p_atac_devil <- devil %>% dplyr::select(geneID,adj_pval_snATAC)
p_atac_glm <- glm %>% dplyr::select(geneID,adj_pval_snATAC)
p_atac <- cbind(p_atac_devil, p_atac_glm)
colnames(p_atac) <- c("geneID1", "adjpval_atac_devil", "geneID2", "adjpval_atac_glm" )

pval_atac <- ggplot(p_atac, aes(x=adjpval_atac_glm, y=adjpval_atac_devil)) + 
  #geom_point(shape=15, color="blue")+
  geom_pointdensity(shape=20) +
  #scale_color_viridis()+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "p-value snATAC")+
  xlab("adjpval snATAC from glmGamPoi") +
  ylab ("adjpval snATAC from Devil") +
  theme_classic()+
  theme(legend.position="none")

p_rna_devil <- devil %>% dplyr::select(geneID,adj_pval_snRNA)
p_rna_glm <- glm %>% dplyr::select(geneID,adj_pval_snRNA)
p_rna <- cbind(p_rna_devil, p_rna_glm)
colnames(p_rna) <- c("geneID1", "adjpval_rna_devil", "geneID2", "adjpval_rna_glm" )

pval_rna<- ggplot(p_rna, aes(x=adjpval_rna_glm, y=adjpval_rna_devil)) + 
  #geom_point(shape=15, color="blue")+
  geom_pointdensity(shape=20) +
  #scale_color_viridis()+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "p-value snRNA")+
  xlab("adjpval snRNA from glmGamPoi") +
  ylab ("adjpval snRNA from Devil") +
  theme_classic()+
  theme(legend.position="none")

lfc_atac_devil <- devil %>% dplyr::select(geneID,lfc_snATAC)
lfc_atac_glm <- glm %>% dplyr::select(geneID,lfc_snATAC)
lfc_atac <- cbind(lfc_atac_devil, lfc_atac_glm)
colnames(lfc_atac) <- c("geneID1", "lfc_atac_devil", "geneID2", "lfc_atac_glm" )

p_lfc_atac <- ggplot(lfc_atac, aes(x=lfc_atac_glm, y=lfc_atac_devil)) + 
  #geom_point(shape=15, color="blue")+
  geom_pointdensity(shape=20) +
  #scale_color_viridis()+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "log2FC snATAC")+
  xlab("log2FC snATAC from glmGamPoi") +
  ylab ("log2FC snATAC from Devil") +
  theme_classic()+
  theme(legend.position="none")


lfc_rna_devil <- devil %>% dplyr::select(geneID,lfc_snRNA)
lfc_rna_glm <- glm %>% dplyr::select(geneID,lfc_snRNA)
lfc_rna <- cbind(lfc_rna_devil, lfc_rna_glm)
colnames(lfc_rna) <- c("geneID1", "lfc_rna_devil", "geneID2", "lfc_rna_glm" )

p_lfc_rna<- ggplot(lfc_rna, aes(x=lfc_rna_glm, y=lfc_rna_devil)) + 
  #geom_point(shape=15, color="blue")+
  geom_pointdensity(shape=20) +
  #scale_color_viridis()+
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "log2FC snRNA")+
  xlab("log2FC snRNA from glmGamPoi") +
  ylab ("log2FC snRNA from Devil") +
  theme_classic()+
  theme(legend.position="none")

grid.arrange(pval_atac,p_lfc_atac,pval_rna,p_lfc_rna,  nrow = 2)


### ROC curve ###

rna_devil <- rna_devil[ rna_devil$name %in% rna_glm$name, ]
rna_glm <- rna_glm[ rna_glm$name %in% rna_devil$name, ]

rna_devil <- rna_devil %>% dplyr::select(name, lfc, adj_pval)
colnames(rna_devil) <- c("name", "lfc_devil", "adjpval_devil")

rna_glm <- rna_glm %>% dplyr::select(name, lfc, adj_pval)
colnames(rna_glm) <- c("name", "lfc_glm", "adjpval_glm")
rna <- merge(rna_devil, rna_glm, by = "name")

library(metaseqR2)
library(tidyverse)

p1 <- matrix(rna$adjpval_devil)
colnames(p1) <- "Devil"
p2 <- matrix(rna$adjpval_glm)
colnames(p2) <- "glmGamPoi"
p <- cbind(p1,p2)

res <- rna %>% 
  mutate(truth = case_when(
    lfc_devil >= 1.0 & adjpval_devil < 0.01 ~ "1",
    lfc_devil <= -1.0 & adjpval_devil < 0.01 ~ "1",
    lfc_devil < 1.0 | lfc_devil > -1.0 & adjpval_devil >= 0.01 ~ "0")) 

truth <- as.vector(as.numeric(res$truth))
names(truth) <- res$name

rna_roc <- metaseqR2::diagplotRoc(truth = truth, p = p1, sig = 0.01, x = "fpr",
                                   y = "tpr", path = NULL, draw = TRUE)
