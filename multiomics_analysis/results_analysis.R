### Results downstream analysis ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR")
sapply(pkgs, require, character.only = TRUE)

grange_path <- "multiomics_analysis/results/grange_annot.RDS"
grange <- readRDS(grange_path)

atac_path <- "multiomics_analysis/results/MuscleATAC/devil_atac.RDS"
atac_res <- readRDS(atac_path) %>% rename(ranges=name)

rna_path <- "multiomics_analysis/results/MuscleRNA/devil_rna.RDS"
rna_res <- readRDS("multiomics_analysis/results/MuscleRNA/devil_rna.RDS") %>% rename(geneID=name)


### Select non duplicated atac genes ###

grange <- grange %>% dplyr::select(SYMBOL, ranges, annotation) 
grange <- grange[ (grange$annotation %in% c("Promoter (<=1kb)")) , ]
atac_res <- dplyr::full_join(atac_res, grange, by = "ranges") %>% 
  dplyr::rename(geneID=SYMBOL) %>% na.omit()

atac_up <- atac_res %>% dplyr::filter(lfc > 0) %>% 
  dplyr::group_by(geneID) %>% 
  dplyr::arrange(lfc)

atac_up <- atac_up[order(atac_up$lfc, decreasing = TRUE), ] 
atac_up_nodup <- atac_up[!duplicated(atac_up$geneID), ]

atac_down <- atac_res %>% filter(lfc < 0) %>% 
  dplyr::group_by(geneID) %>% 
  dplyr::arrange(lfc)

atac_down <- atac_down[order(atac_down$lfc), ] 
atac_down_nodup <- atac_down[!duplicated(atac_down$geneID), ]

atac_nodup <- rbind(atac_up_nodup, atac_down_nodup)
atac_nodup <- atac_nodup[!duplicated(atac_nodup$geneID), ]

atac_deg <- atac_nodup %>% 
  dplyr::filter(adj_pval < 0.05 & abs(lfc) > 0.5)

rna_deg <- rna_res %>% 
  dplyr::filter(adj_pval < 0.05 & abs(lfc) > 0.5)


### scaDA results analysis ###

grange_scaDA_path <- "multiomics_analysis/results/grange_annot_scADA.RDS"
grange_scaDA <- readRDS(grange_scaDA_path)

grange <- "multiomics_analysis/results/grange_annot.RDS"
grange <- readRDS(grange)

atac_scaDA_path <- "multiomics_analysis/results/scADA_res.RDS"
atac_scaDA <- readRDS(atac_scaDA_path) 

rna_devil <- "multiomics_analysis/results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "multiomics_analysis/results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_edge <- "multiomics_analysis/results/MuscleRNA/edge_rna.RDS"
rna_edge <- readRDS(rna_edge)
#rna_edge <- topTags(rna_edge, n=15563)$table
rna_edge$geneID <- rownames(rna_edge)
colnames(rna_edge) <- c("lfc", "stError", "tval", "adj_pval", "geneID")

grange_scaDA <- grange_scaDA %>% dplyr::select(SYMBOL) %>% 
  dplyr::rename(geneID=SYMBOL)

atac_scaDA <- cbind(grange_scaDA, atac_scaDA)
atac_scaDA$log2fc <- -1 * atac_scaDA$log2fc

grange <- grange[ (grange$annotation %in% c("Promoter (<=1kb)")) , ]
atac_scaDA <- atac_scaDA[ (atac_scaDA$geneID %in% grange$SYMBOL) , ]

atac_up <- atac_scaDA %>% dplyr::filter(log2fc > 0) %>% 
  dplyr::group_by(geneID) %>% 
  dplyr::arrange(log2fc)

atac_up <- atac_up[order(atac_up$log2fc, decreasing = TRUE), ] 
atac_up_nodup <- atac_up[!duplicated(atac_up$geneID), ]

atac_down <- atac_scaDA %>% dplyr::filter(log2fc < 0) %>% 
  dplyr::group_by(geneID) %>% 
  dplyr::arrange(log2fc)

atac_down <- atac_down[order(atac_down$log2fc), ] 
atac_down_nodup <- atac_down[!duplicated(atac_down$geneID), ]

atac_nodup <- rbind(atac_up_nodup, atac_down_nodup)
atac_nodup <- atac_nodup[!duplicated(atac_nodup$geneID), ]

saveRDS(atac_nodup, file = "multiomics_analysis/results/atac_nodup_scaDA.RDS")

# Gene selection based on LFC & pvalue cutoff #
atac_deg <- atac_nodup %>% 
  dplyr::filter(FDR < 0.05 & abs(log2fc) >= 1.0)

rna_deg_devil <- rna_devil %>% 
  dplyr::filter(adj_pval < 0.05 & abs(lfc) >= 1.0)

rna_deg_glm <- rna_glm %>% 
  dplyr::filter(adj_pval < 0.05 & abs(lfc) >= 1.0)

rna_deg_edge <- rna_edge %>% 
  dplyr::filter(adj_pval < 0.05 & abs(lfc) > 1.0)


# Gene selection based on pvalue cutoff #
atac_deg <- atac_nodup %>% 
  dplyr::filter(FDR < 0.05)

rna_deg_devil <- rna_devil %>% 
  dplyr::filter(adj_pval < 0.05)

rna_deg_glm <- rna_glm %>% 
  dplyr::filter(adj_pval < 0.05)

rna_deg_edge <- rna_edge %>% 
  dplyr::filter(adj_pval < 0.05)

# Vienn diagram #
devil <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_devil$geneID)
glm <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_glm$geneID)
edge <- list(snATAC=atac_deg$geneID, snRNA=rna_deg_edge$geneID)

p1 <- ggvenn::ggvenn(devil, c("snATAC", "snRNA"), show_percentage = FALSE,
             set_name_size = 4)
p1 <- p1 + ggplot2::labs(title="Devil") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p1

p2 <- ggvenn::ggvenn(glm, c("snATAC", "snRNA"), show_percentage = FALSE,
             set_name_size = 4)
p2 <- p2 + ggplot2::labs(title="glmGamPoi") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p2

p3 <- ggvenn::ggvenn(edge, c("snATAC", "snRNA"), show_percentage = FALSE,
                     set_name_size = 4)
p3 <- p3 + ggplot2::labs(title="edgeR") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p3

gridExtra::grid.arrange(p1, p2, p3, nrow = 1)


### Filter overlapped genes ###
atac_deg_devil <- atac_deg[ atac_deg$geneID %in% rna_deg_devil$geneID,]
rna_deg_devil <- rna_deg_devil[ rna_deg_devil$geneID %in% atac_deg$geneID,]
atac_deg_devil <- atac_deg_devil %>% dplyr::select(geneID, FDR, log2fc)
colnames(atac_deg_devil) <- c("geneID", "adj_pval_snATAC", "lfc_snATAC")
rna_deg_devil <- rna_deg_devil %>% dplyr::select(geneID, adj_pval, lfc)
colnames(rna_deg_devil) <- c("geneID", "adj_pval_snRNA", "lfc_snRNA")
overlap_devil <- dplyr::full_join(atac_deg_devil, rna_deg_devil, by = "geneID")

saveRDS(overlap_devil, file = "multiomics_analysis/results/overlap_devil.RDS")

atac_deg_glm <- atac_deg[ atac_deg$geneID %in% rna_deg_glm$geneID,]
rna_deg_glm <- rna_deg_glm[ rna_deg_glm$geneID %in% atac_deg_glm$geneID,]
atac_deg_glm<- atac_deg_glm %>% dplyr::select(geneID, FDR, log2fc)
colnames(atac_deg_glm) <- c("geneID", "adj_pval_snATAC", "lfc_snATAC")
rna_deg_glm <- rna_deg_glm %>% dplyr::select(geneID, adj_pval, lfc)
colnames(rna_deg_glm) <- c("geneID", "adj_pval_snRNA", "lfc_snRNA")
overlap_glm <- dplyr::full_join(atac_deg_glm, rna_deg_glm, by = "geneID")

saveRDS(overlap_glm, file = "multiomics_analysis/results/overlap_glm.RDS")

atac_deg_edge <- atac_deg[ atac_deg$geneID %in% rna_deg_edge$geneID,]
rna_deg_edge <- rna_deg_edge[ rna_deg_edge$geneID %in% atac_deg_edge$geneID,]
atac_deg_edge<- atac_deg_edge %>% dplyr::select(geneID, log2fc, FDR)
colnames(atac_deg_edge) <- c("geneID", "lfc_snATAC", "adj_pval_snATAC")
rna_deg_edge <- rna_deg_edge %>% dplyr::select(geneID, lfc, adj_pval)
colnames(rna_deg_edge) <- c("geneID", "lfc_snRNA", "adj_pval_snRNA")
overlap_edge <- dplyr::full_join(atac_deg_edge, rna_deg_edge, by = "geneID")


### Correlation plot ###

corr1 <- ggplot2::ggplot(mapping = aes(x = -log10(overlap_devil$adj_pval_snRNA), y = -log10(overlap_devil$adj_pval_snATAC))) +
  geom_point(shape = 21, fill = 'black', size = 2) +
  labs(title = "Devil")+
  xlab("snRNA -log10(adj_pval)") +
  ylab ("snATAC -log10(adj_pval)") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  sm_statCorr()+
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()

corr2 <- ggplot2::ggplot(mapping = aes(x = -log10(overlap_glm$adj_pval_snRNA), y = -log10(overlap_glm$adj_pval_snATAC))) +
  geom_point(shape = 21, fill = 'black', size = 2) +
  labs(title = "glmGamPoi")+
  xlab("snRNA -log10(adj_pval)") +
  ylab ("snATAC -log10(adj_pval)") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  sm_statCorr()+
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()

corr3 <- ggplot2::ggplot(mapping = aes(x = -log10(overlap_edge$adj_pval_snRNA), y = -log10(overlap_edge$adj_pval_snATAC))) +
  geom_point(shape = 21, fill = 'black', size = 2) +
  labs(title = "edgeR")+
  xlab("snRNA -log10(adj_pval)") +
  ylab ("snATAC -log10(adj_pval)") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  sm_statCorr()+
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()

corr1+corr2+corr3



### Gene set Enrichement analysis ###

#rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

# Convert Gene symbols to EntrezID #
hs <- org.Hs.eg.db
my_symbols <- overlap_glm$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")

gene_list <- na.omit(gene_list)
gene_list <- gene_list[!duplicated(gene_list$SYMBOL),]
gene_list <- gene_list[gene_list$SYMBOL %in% overlap_glm$geneID,]
gene_list_rank <- as.vector(overlap_glm$lfc_snRNA)
names(gene_list_rank) <- gene_list$ENTREZID
gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

gseGO <- clusterProfiler::gseGO(gene_list_rank, ont = "All", OrgDb = org.Hs.eg.db,
                                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)

# Select enriched pathways #
res_gse <- gseGO@result

res_gse <- res_gse %>%
  filter(Description %in% c("cellular modified amino acid metabolic process","actin filament binding", "sarcoplasm",
                            "actin binding", "response to BMP","response to fibroblast growth factor", 
                             "negative regulation of mitotic cell cycle", "cellular response to oxygen levels",
                            "MAPK cascade", "response to calcium ion", "ERK1 and ERK2 cascade",
                            "CD4-positive, alpha-beta T cell activation")) %>% 
  mutate(gene_clusters = case_when(
  NES > 0  ~ 'up-regulated',
  NES < 0  ~ 'down-regulated'))

#res_gse <- res_gse %>%
  #filter(Description %in% c("positive regulation of cellular component biogenesis","MAPK cascade", 
                            #"regulation of apoptotic process","actin filament-based movement", 
                            #"actin-mediated cell contraction", "muscle contraction", "muscle system process", 
                            #"response to oxygen-containing compound")) %>% 
  #mutate(gene_clusters = case_when(
    #NES > 0  ~ 'up-regulated',
    #NES < 0  ~ 'down-regulated'))

### Visualize enrichment results ###

res_gse$log_padjust <- -log10(res_gse$p.adjust)

plot1 <- ggpubr::ggdotchart(res_gse, x = "Description", y = "log_padjust",
                    color = "gene_clusters",
                    palette = c("blue", "#FC4E07"),
                    sorting = "descending",
                    rotate = TRUE,
                    group = "gene_clusters",
                    dot.size = "setSize",
                    add = "segments",
                    title = "GO enrichment - Devil",
                    xlab = "GO Pathways",
                    ylab = "-log10(padjust)",
                    ggtheme = theme_pubr()
)
plot1 + theme(legend.position = "right")+
  theme_bw()+
  font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 10, color = "black", face = "bold")+
  font("xlab", size = 10)+
  font("ylab", size = 10)

#res_gse_devil <- res_gse
#res_gse_glm <- res_gse
res_gse_devil <- res_gse_devil %>% mutate(method = "Devil")
res_gse_glm <- res_gse_glm %>% mutate(method = "glmGamPoi")
res_gse <- rbind(res_gse_devil, res_gse_glm)
res_gse$method <- as.factor(res_gse$method)

# Barplot #

p1 <- ggplot(res_gse, aes(x = log_padjust, y = Description, fill=method)) + 
  geom_bar(stat="identity", position="dodge", width = 0.50)+
  labs(x = "-Log10(pvalue)",
       y = "GO Pathways",
       title = "GO analysis")+
  theme(strip.text.x = element_text(size=12, color="black", face="bold.italic"))+
  theme_bw()
p1 <- p1 + scale_fill_manual(values=c("#0072B2", "#CC79A7"))

p1 <- p1 + font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 10, color = "black", face = "bold")+
  font("xlab", size = 10)+
  font("ylab", size = 10)
p1



# Load 10X count matrices from multiple directories #
library(Seurat)
library(scCustomize)

base_path <- "/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/data/multiomics/rna/counts_data/"

counts_mtx <- scCustomize::Read10X_Multi_Directory(base_path = base_path, 
                                                   default_10X = FALSE, 
                                                   secondary_path = "gex_matrices",
                                                   merge = TRUE)



#metadata_1$`_index` %>% stringr::str_remove(pattern = "._[0-9]_[0-9]_[0-9]_[0-9]$")
#metadata <- metadata[ (metadata$age_pop %in% c("young_pop", "old_pop") & metadata$Annotation %in% c("Type I", "Type II")),]

# Formatting cell indexes #

metadata <- metadata[ (metadata$tech %in% c("snRNA")) ,]
#metadata <- metadata %>% 
  #mutate(cell_index = paste(metadata$orig.ident, metadata$`_index`, sep = "_"))
metadata_1 <- metadata[ (metadata$sample %in% c("OM1", "OM2", "OM3")) ,]
metadata_1$`_index` <- substr(metadata_1$`_index`, 1, 21)

metadata_2 <- metadata[ (metadata$sample %in% c("OM4", "OM5", "OM6", "OM7", "OM8", "OM9")) ,]
metadata_2$`_index` <- substr(metadata_2$`_index`, 1, 11)
metadata_2$`_index` <- gsub("._+[0-9]+_$", "", metadata_2$`_index`)
metadata_2$`_index` <- metadata_2$`_index` %>% stringr::str_remove(pattern = "_")
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL32N32"] <- "CELL32N3"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10080N"] <- "CELL10080N1"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL77N31"] <- "CELL77N3"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10209N"] <- "CELL10209N1"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10077N"] <- "CELL10077N1"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10226N"] <- "CELL10226N1"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10396N"] <- "CELL10396N1"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10210N"] <- "CELL10210N1"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL2N"] <- "CELL2N2"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL77N22"] <- "CELL77N2"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL38N32"] <- "CELL38N3"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL61N22"] <- "CELL61N2"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL83N22"] <- "CELL83N2"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL88N32"] <- "CELL88N3"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL47N22"] <- "CELL47N2"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL6N"] <- "CELL6N2"
names(metadata_2$`_index`)[names(metadata_2$`_index`) == "CELL10317N"] <- "CELL10317N1"


metadata_3 <- metadata[ (metadata$sample %in% c("P17", "P21", "P23", "P27", "P29", "P3", "P5")) ,]
metadata_3$`_index` <- substr(metadata_3$`_index`, 1, 11)
metadata_3$`_index` <- gsub("._+[0-9]+_$", "", metadata_3$`_index`)
metadata_3$`_index` <- sub("_2", "", metadata_3$`_index`)
metadata_3$`_index` <- metadata_3$`_index` %>% stringr::str_remove(pattern = "_")

metadata_4 <- metadata[ (metadata$sample %in% c("YM1", "YM2", "YM3", "YM4")) ,]
metadata_4$`_index` <- substr(metadata_4$`_index`, 1, 11)
metadata_4$`_index` <- gsub("._+[0-9]+_$", "", metadata_4$`_index`)
metadata_4$`_index` <- sub("_2", "", metadata_4$`_index`)
metadata_4$`_index` <- metadata_4$`_index` %>% stringr::str_remove(pattern = "_")

metadata <- rbind(metadata_1, metadata_2, metadata_3, metadata_4)

metadata <- metadata %>% 
  mutate(cell_index = paste(metadata$orig.ident, metadata$`_index`, sep = "_"))

metadata_adj_rna <- metadata_adj_rna[ (metadata_adj_rna$Annotation %in% c("Type I", "Type II")) ,]
metadata_adj_rna$cell_index <- metadata_adj_rna$cell_index %>% stringr::str_remove(pattern = "_")
rownames(metadata_adj_rna) <- metadata_adj_rna$cell_index


colnames(counts_mtx) <- colnames(counts_mtx) %>% stringr::str_remove(pattern = "_")
counts_mtx <- counts_mtx[ , (colnames(counts_mtx) %in% rownames(metadata)) ]
metadata <- metadata[ (rownames(metadata) %in% colnames(counts_mtx)), ]
counts_mtx <- counts_mtx[!duplicated(rownames(counts_mtx)), ]
metadata_adj_rna <- metadata_adj_rna[!duplicated(metadata_adj_rna$cell_index), ]
rownames(metadata_adj_rna) <- metadata_adj_rna$cell_index

saveRDS(metadata_adj_rna, file = "multiomics_analysis/results/metadata_adj_rna.RDS")
