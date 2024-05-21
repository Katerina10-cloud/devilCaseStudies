# Plot results #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram")
sapply(pkgs, require, character.only = TRUE)

source("utils.R")

set.seed(12345)

# UMAP #

metadata_rna <- metadata_rna[ (metadata_rna$tech_id %in% c("1")), ]
metadata_rna$cell_cluster <- as.numeric(as.factor(metadata_rna$cell_cluster))
#Seurat::DiscretePalette()
umap_rna <- metadata_rna %>% dplyr::select(umap_1, umap_2, cell_type)
umap_rna$cell_type <- as.factor(umap_rna$cell_type)

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

overlap <- draw.pairwise.venn(area1=1507, area2=1034,cross.area=340, 
                              category=c("snATAC","snRNA"),fill=c("Red","Blue"))


# Correlation plot #
cor_coefs <- cor.test(res_atac_rna$lfc_snRNA, res_atac_rna$lfc_snATAC)
corr_plot <- ggplot(data = res_atac_rna, aes(x = lfc_snRNA, y = lfc_snATAC)) + 
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  annotate("text", x = 10, y = 4, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 10, y = 3.5, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))
corr_plot

corr_plot <- ggplot(mapping = aes(x = res_atac_rna$lfc_snRNA, y = res_atac_rna$lfc_snATAC)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  labs(title = "Correlation of significantly DAG (snATAC-seq) and GE (snRNA-seq)")+
  xlab("LogFC2 snRNA") +
  ylab ("LogFC2 snATAC") +
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  sm_statCorr()

