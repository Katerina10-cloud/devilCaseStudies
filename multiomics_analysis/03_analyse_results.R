### Results downstream analysis ###

rm(list = ls())

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram')
sapply(pkgs, require, character.only = TRUE)

rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))


# Gene selection based on LFC & pvalue cutoff #
lfc_cut <- 0.5
pval_cut <- .05

rna_deg_devil <- rna_devil %>%
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>% 
  dplyr::mutate(method = "devil")

rna_deg_glm <- rna_glm %>%
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>% 
  dplyr::mutate(method = "glmGamPoi")

rna_deg_nebula <- rna_nebula %>%
  dplyr::filter(adj_pval < pval_cut, abs(lfc) > lfc_cut) %>% 
  dplyr::mutate(method = "nebula")



# Vienn diagram ####

x <- list(
  devil = rna_deg_devil$geneID,
  glmGamPoi = rna_deg_glm$geneID,
  nebula = rna_deg_nebula$geneID
)

venn_plot <- ggVennDiagram::ggVennDiagram(x, color = 1, lwd = 0.7) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
venn_plot

saveRDS(venn_plot, "plot/venn_plot.rds")


# Volcano plot
rna_deg_devil <- rna_deg_devil %>% dplyr::filter(adj_pval > 4.467475e-90 ,)
rna_deg_nebula <- rna_deg_nebula %>% dplyr::filter(adj_pval > 4.787713e-22 ,)

rna_join <- rbind(rna_deg_devil, rna_deg_glm, rna_deg_nebula)
rna_join <- rna_join %>% 
  dplyr::mutate(
    isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut),
    DEtype = if_else(!isDE, "n.s.", if_else(lfc > 0, "Up-reg", "Down-reg")))

de_colors <- c("Down-reg" = "steelblue", "Up-reg" = "indianred", "n.s." = "grey")

gene_markers <- c("TNNT1", "MYH7", "MYH7B", "TNNT2", "PDE4B", "JUN", "FOSB",
                  "ID1", "MDM2", "TNNT3", "MYH2", "MYH1", "ENOX1", "SAA2", "SAA1",
                  "DCLK1", "ADGRB3", "NCAM1", "COL22A1", "PHLDB2", "CHRNE")

p_volcanos <- rna_join %>%
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(col = DEtype), size = 2.0, alpha = 0.2) + 
  scale_color_manual(values = de_colors) + 
  geom_label_repel(
    data = rna_join %>% filter(geneID %in% gene_markers),  
    aes(label = geneID),
    size = 4.0,               
    fontface = "bold",        
    color = "black",          
    fill = "white",       
    box.padding = 0.8,        
    point.padding = 0.5,      
    max.overlaps = Inf,       
    segment.color = "black", 
    segment.size = 0.5,       
    label.padding = unit(0.15, "lines"), 
    label.r = unit(0.4, "lines"),       
    min.segment.length = 0    
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)), 
                                  ceiling(max(rna_join$lfc)), by = 2)) +
  facet_wrap(~factor(method, levels = c("devil", "glmGamPoi", "nebula")), nrow = 1, scales = "free") +
  labs(x = expression(Log[2] ~ FC), 
       y = expression(-log[10] ~ Pvalue), 
       col = "DE type") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(
    legend.position = '',  # Hide legend (can be adjusted if needed)
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, color = "black"),  
    strip.text = element_text(size = 16, face = "plain", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
p_volcanos

ggsave("plot/volcano_all_methods.pdf", dpi = 400, width = 16.0, height = 5.5, plot = p_volcanos)



# Gene set Enrichement analysis devil ####

pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

de_res <- rna_deg_devil
de_res$RankMetric <- -log10(de_res$adj_pval) * sign(de_res$lfc)
de_res <- de_res %>% arrange(-RankMetric)

genes <- de_res$RankMetric
names(genes) <- de_res$geneID

hs <- org.Hs.eg.db

gseGO <- clusterProfiler::gseGO(
  genes,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,
  maxGSSize = 350,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

gseGO_devil <- gseGO@result %>% as.data.frame()

saveRDS(gseGO_devil, file = "results/gsea_GO/gseGO_devil.RDS")



# Gene set Enrichement analysis glmGamPoi ####

de_res <- rna_deg_glmGamPoi
de_res$RankMetric <- -log10(de_res$adj_pval) * sign(de_res$lfc)
de_res <- de_res %>% arrange(-RankMetric)

genes <- de_res$RankMetric
names(genes) <- de_res$geneID

hs <- org.Hs.eg.db

gseGO <- clusterProfiler::gseGO(
  genes,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,
  maxGSSize = 350,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

gseGO_glmGamPoi <- gseGO@result %>% as.data.frame()

saveRDS(gseGO_glmGamPoi, file = "results/gsea_GO/gseGO_glmGamPoi.RDS")


# Gene set Enrichement analysis nebula ####
de_res <- rna_deg_nebula
de_res$RankMetric <- -log10(de_res$adj_pval) * sign(de_res$lfc)
de_res <- de_res %>% arrange(-RankMetric)

genes <- de_res$RankMetric
names(genes) <- de_res$geneID

hs <- org.Hs.eg.db

gseGO <- clusterProfiler::gseGO(
  genes,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  minGSSize = 10,
  maxGSSize = 350,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

gseGO_nebula <- gseGO@result %>% as.data.frame()

saveRDS(gseGO_nebula, file = "results/gsea_GO/gseGO_nebula.RDS")


# Enrichment results preprocessing #

filtered_devil <- gseGO_devil %>%
  dplyr::filter(enrichmentScore < -0.3 | enrichmentScore > 0.4)

filtered_glmGamPoi <- gseGO_glmGamPoi %>%
  dplyr::filter(enrichmentScore < -0.3 | enrichmentScore > 0.4)

filtered_nebula <- gseGO_nebula %>%
  dplyr::filter(enrichmentScore < -0.3 | enrichmentScore > 0.4)


# Remove redundant terms #

filtered_devil <- filtered_devil %>%
  mutate(genes = strsplit(core_enrichment, "/"))


n_terms <- nrow(filtered_devil)
overlap_matrix <- matrix(0, nrow = n_terms, ncol = n_terms,
                         dimnames = list(filtered_devil$Description, filtered_devil$Description))

for (i in 1:n_terms) {
  for (j in i:n_terms) {
    shared_genes <- length(intersect(filtered_devil$genes[[i]], filtered_devil$genes[[j]]))
    total_genes <- length(union(filtered_devil$genes[[i]], filtered_devil$genes[[j]]))
    jaccard_index <- shared_genes / total_genes
    
    # Fill overlap matrix with Jaccard index
    overlap_matrix[i, j] <- jaccard_index
    overlap_matrix[j, i] <- jaccard_index
  }
}


redundant_terms <- c()
threshold <- 0.5

for (i in 1:(n_terms - 1)) {
  for (j in (i + 1):n_terms) {
    if (overlap_matrix[i, j] > threshold) {
      redundant_terms <- c(redundant_terms, filtered_devil$Description[j])
    }
  }
}

non_redundant_terms <- filtered_devil %>%
  filter(!Description %in% redundant_terms)
