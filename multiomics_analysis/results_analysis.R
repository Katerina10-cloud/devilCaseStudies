### Results downstream analysis ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

grange_path <- "multiomics_analysis/results/grange_annot.RDS"
grange <- readRDS(grange_path)

atac_path <- "multiomics_analysis/results/MuscleATAC/devil_atac.RDS"
atac_res <- readRDS(atac_path) %>% rename(ranges=name)

rna_path <- "multiomics_analysis/results/MuscleRNA/devil_rna.RDS"
rna_res <- readRDS("multiomics_analysis/results/MuscleRNA/devil_rna.RDS") %>% rename(geneID=name)

### Select non duplicated atac genes ###

grange <- grange %>% select(SYMBOL, ranges) 
atac_res <- full_join(atac_res, grange, by = "ranges") %>% rename(geneID=SYMBOL)

atac_up <- subset(atac_res, atac_res$lfc > 0) %>% 
  group_by(geneID) %>% 
  arrange(lfc)

atac_up <- atac_up[order(atac_up$lfc, decreasing = TRUE), ] 
atac_up_nodup <- atac_up[!duplicated(atac_up$geneID), ]

atac_down <- subset(atac_res, atac_res$lfc < 0) %>% 
  group_by(geneID) %>% 
  arrange(lfc)

atac_down <- atac_down[order(atac_down$lfc), ] 
atac_down_nodup <- atac_down[!duplicated(atac_down$geneID), ]

atac_nodup <- rbind(atac_up_nodup, atac_down_nodup)
#atac_nodup <- atac_nodup[!duplicated(atac_nodup$geneID), ]

saveRDS(atac_nodup, file = "multiomics_analysis/results/MuscleATAC/atac_nodup.RDS")


### Find overlapped genes ###

atac_up <- subset(atac_nodup, atac_nodup$adj_pval < 0.05 & atac_nodup$lfc >= 1.0)
atac_down <- subset(atac_nodup, atac_nodup$adj_pval < 0.05 & atac_nodup$lfc <= -1.0)
atac_deg <- rbind(atac_up, atac_down)

rna_up <- subset(rna_res, rna_res$adj_pval < 0.05 & rna_res$lfc >= 1.0)
rna_down <- subset(rna_res, rna_res$adj_pval < 0.05 & rna_res$lfc <= -1.0)
rna_deg <- rbind(rna_up, rna_down)

# Vienn diagram #
devil <- list(snATAC=atac_deg$geneID, snRNA=rna_deg$geneID)

p1 <- ggvenn(devil, c("snATAC", "snRNA"), show_percentage = FALSE,
             set_name_size = 4)
p1 <- p1 + ggplot2::labs(title="Devil") +
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))


atac_deg <- atac_deg[ atac_deg$geneID %in% rna_deg$geneID,]
rna_deg <- rna_deg[ rna_deg$geneID %in% atac_deg$geneID,]

atac_deg <- atac_deg %>% select(geneID, adj_pval, lfc)
colnames(atac_deg) <- c("geneID", "adj_pval_snATAC", "lfc_snATAC")

rna_deg <- rna_deg %>% select(geneID, adj_pval, lfc)
colnames(rna_deg) <- c("geneID", "adj_pval_snRNA", "lfc_snRNA")

overlap_devil <- full_join(atac_deg, rna_deg, by = "geneID")

saveRDS(overlap_devil, file = "multiomics_analysis/results/overlap_devil.RDS")


### Correlation plot ###

corr_plot <- ggplot2::ggplot(mapping = aes(x = overlap_devil$lfc_snRNA, y = overl_devil$lfc_snATAC)) +
  geom_point(shape = 21, fill = 'black', size = 2) +
  labs(title = "devil")+
  xlab("snRNA log2FC") +
  ylab ("snATAC log2FC") +
  geom_smooth(method='lm',formula=y~x, color="red", fill="black", se=TRUE)+
  sm_statCorr()+
  geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  theme_classic()
corr_plot


### Gene set Enrichement analysis ###

rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

# Convert Gene symbols to EntrezID #
hs <- org.Hs.eg.db
my_symbols <- overlap_devil$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")

gene_list <- na.omit(gene_list)
gene_list <- gene_list[!duplicated(gene_list$SYMBOL),]
gene_list <- gene_list[gene_list$SYMBOL %in% overlap_devil$geneID,]
gene_list_rank <- as.vector(overlap_devil$lfc_snRNA)
names(gene_list_rank) <- gene_list$ENTREZID
gene_list_rank <- sort(gene_list_rank, decreasing = TRUE)

gseGO <- clusterProfiler::gseGO(geneList = gene_list_rank, ont = "BP", OrgDb = org.Hs.eg.db,
                                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)

# Select enriched pathways #
res_gse <- gseGO@result
res_gse <- res_gse %>%
  filter(Description %in% c("positive regulation of cellular component biogenesis","MAPK cascade", 
                            "regulation of apoptotic process","actin filament-based movement", 
                            "actin-mediated cell contraction", "muscle contraction", "muscle system process", 
                            "response to oxygen-containing compound")) %>% 
  mutate(gene_clusters = case_when(
    NES > 0  ~ 'up-regulated',
    NES < 0  ~ 'down-regulated'))

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


#res_gse_devil <- res_gse_devil %>% mutate(method = "Devil")
#res_gse_glm <- res_gse_glm %>% mutate(method = "glmGamPoi")
#res_gse <- rbind(res_gse_devil, res_gse_glm)
#res_gse$method <- as.factor(res_gse$method)

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



