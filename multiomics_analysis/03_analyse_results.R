### Results downstream analysis ###

rm(list = ls())

#setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis")

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

venn_plot <- ggVennDiagram::ggVennDiagram(x, color = 1, lwd = 0.8) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
venn_plot

saveRDS(venn_plot, "plot/venn_plot.rds")

# Volcano plot
rna_deg_devil <- rna_deg_devil %>% dplyr::filter(adj_pval > 4.467475e-90 ,)
rna_deg_nebula <- rna_deg_nebula %>% dplyr::filter(adj_pval > 4.787713e-22 ,)
rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])

rna_join <- rbind(rna_deg_devil, rna_deg_glm, rna_deg_nebula)
rna_join <- rna_join %>%
  dplyr::mutate(
    isDE = (abs(lfc) >= lfc_cut) & (adj_pval <= pval_cut),
    DEtype = if_else(!isDE, "n.s.", if_else(lfc > 0, "Up-reg", "Down-reg")))

de_colors <- c("Down-reg" = "steelblue", "Up-reg" = "indianred", "n.s." = "grey")

gene_markers <- c("TNNT1", "MYH7", "MYH7B", "TNNT2", "PDE4B", "JUN", "FOSB",
                  "ID1", "MDM2", "TNNT3", "MYH2", "MYH1", "ENOX1", "SAA2", "SAA1",
                  "DCLK1", "ADGRB3", "NCAM1", "COL22A1", "PHLDB2", "CHRNE")

saveRDS(gene_markers, file = "results/gsea_GO/gene_markers_BRelevant.RDS")

p_volcanos <- rna_join %>%
  ggplot(mapping = aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(aes(col = DEtype), size = 2.0, alpha = 0.2) +
  scale_color_manual(values = de_colors) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  geom_label_repel(
    data = rna_join %>% filter(geneID %in% gene_markers),
    aes(label = geneID),
    #size = 4.0,
    #fontface = "bold",
    color = "black",
    fill = "white",
    # box.padding = 0.8,
    # point.padding = 0.5,
    max.overlaps = Inf
    # segment.color = "black",
    # segment.size = 0.5,
    # label.padding = unit(0.15, "lines"),
    # label.r = unit(0.4, "lines"),
    # min.segment.length = 0
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)),
                                  ceiling(max(rna_join$lfc)), by = 2)) +
  facet_wrap(~factor(method, levels = c("devil", "glmGamPoi", "nebula")), nrow = 1, scales = "free") +
  labs(x = expression(Log[2] ~ FC),
       y = expression(-log[10] ~ Pvalue),
       col = "DE type") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
p_volcanos

ggsave("plot/volcano_all_methods.pdf", dpi = 400, width = 16.0, height = 5.5, plot = p_volcanos)
saveRDS(p_volcanos, "plot/volcanos.rds")


### Gene set Enrichement analysis ###

pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr")
sapply(pkgs, require, character.only = TRUE)

enrichmentGO <- function(rna_deg_data) {
  # RankMetric Calculation
  rna_deg_data$RankMetric <- -log10(rna_deg_data$adj_pval) * sign(rna_deg_data$lfc)
  rna_deg_data <- rna_deg_data %>% arrange(-RankMetric)

  genes <- rna_deg_data$RankMetric
  names(genes) <- rna_deg_data$geneID

  gseGO <- clusterProfiler::gseGO(
    genes,
    ont = "BP",  # use BP for dotplot, 'All' to perfom gene classification
    OrgDb = org.Hs.eg.db,
    minGSSize = 10,
    maxGSSize = 350,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    verbose = TRUE,
    seed = 1234
  )

  return(gseGO@result %>% as.data.frame())
}


rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])

gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)


# Remove redundant terms #

remove_redundant_terms <- function(data, enrichment_col = "enrichmentScore", core_col = "core_enrichment", desc_col = "Description", threshold = 0.5) {
  filtered_data <- data %>%
    #dplyr::filter(!!sym(enrichment_col) < -0.3 | !!sym(enrichment_col) > 0.4) %>%
    mutate(genes = strsplit(!!sym(core_col), "/"))

  n_terms <- nrow(filtered_data)
  overlap_matrix <- matrix(0, nrow = n_terms, ncol = n_terms,
                           dimnames = list(filtered_data[[desc_col]], filtered_data[[desc_col]]))

  for (i in 1:n_terms) {
    for (j in i:n_terms) {
      shared_genes <- length(intersect(filtered_data$genes[[i]], filtered_data$genes[[j]]))
      total_genes <- length(union(filtered_data$genes[[i]], filtered_data$genes[[j]]))
      jaccard_index <- shared_genes / total_genes

      # Fill overlap matrix with Jaccard index
      overlap_matrix[i, j] <- jaccard_index
      overlap_matrix[j, i] <- jaccard_index
    }
  }
  redundant_terms <- c()
  for (i in 1:(n_terms - 1)) {
    for (j in (i + 1):n_terms) {
      if (overlap_matrix[i, j] > threshold) {
        redundant_terms <- c(redundant_terms, filtered_data[[desc_col]][j])
      }
    }
  }
  non_redundant_data <- filtered_data %>%
    filter(!(!!sym(desc_col) %in% redundant_terms))

  return(non_redundant_data)
}


filtered_devil_nonRed <- remove_redundant_terms(gseGO_devil,
                                                       enrichment_col = "enrichmentScore",
                                                       core_col = "core_enrichment",
                                                       desc_col = "Description")

filtered_glm_nonRed <- remove_redundant_terms(gseGO_glm,
                                                       enrichment_col = "enrichmentScore",
                                                       core_col = "core_enrichment",
                                                       desc_col = "Description")

filtered_nebula_nonRed <- remove_redundant_terms(gseGO_nebula,
                                                    enrichment_col = "enrichmentScore",
                                                    core_col = "core_enrichment",
                                                    desc_col = "Description")


saveRDS(filtered_devil_nonRed, file = "results/gsea_GO/gseGO_devil.RDS")
saveRDS(filtered_glm_nonRed, file = "results/gsea_GO/gseGO_glmGamPoi.RDS")
saveRDS(filtered_nebula_nonRed, file = "results/gsea_GO/gseGO_nebula.RDS")


# Enrichment plot #

filtered_devil_nonRed <- readRDS("results/gsea_GO/gseGO_devil.RDS")
filtered_glm_nonRed <- readRDS("results/gsea_GO/gseGO_glmGamPoi.RDS")
filtered_nebula_nonRed <- readRDS("results/gsea_GO/gseGO_nebula.RDS")

filtered_devil_nonRed$DE_type <- ifelse(filtered_devil_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")
filtered_glm_nonRed$DE_type <- ifelse(filtered_glm_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")
filtered_nebula_nonRed$DE_type <- ifelse(filtered_nebula_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")

filtered_devil_nonRed <- filtered_devil_nonRed %>%
  dplyr::mutate(method = "devil")

filtered_glm_nonRed <- filtered_glm_nonRed %>%
  dplyr::mutate(method = "glmGamPoi")

filtered_nebula_nonRed <- filtered_nebula_nonRed %>%
  dplyr::mutate(method = "nebula")


# Select pathways to plot

terms_devil <- c("myofibril assembly", "actin-mediated cell contraction", "muscle cell development",
                 "muscle contraction", "acute-phase response", "regulation of cell division", "cell surface pattern recognition receptor signaling pathway", "homophilic cell adhesion via plasma membrane adhesion molecules",
                 "regulation of G2/M transition of mitotic cell cycle", "positive regulation of ERK1 and ERK2 cascade",
                 "immune response-activating cell surface receptor signaling pathway",
                 "tRNA metabolic process", "positive regulation of cytokine production", "epithelial cell proliferation",
                 "regulation of leukocyte activation", "regulation of MAPK cascade", "defense response")

terms_glm <- c("myofibril assembly", "actin filament-based movement", "muscle cell development", "muscle contraction",
               "acute-phase response", "homophilic cell adhesion via plasma membrane adhesion molecules", "regulation of cell division",
               "positive regulation of ERK1 and ERK2 cascade", "immune response-activating cell surface receptor signaling pathway",
               "tRNA metabolic process", "positive regulation of cytokine production", "epithelial cell proliferation", "defense response")

filtered_devil_nonRed <- filtered_devil_nonRed %>%
  dplyr::filter(Description %in% terms_devil)

filtered_glm_nonRed <- filtered_glm_nonRed %>%
  dplyr::filter(Description %in% terms_glm)

data_join <- rbind(filtered_devil_nonRed, filtered_glm_nonRed, filtered_nebula_nonRed)

plot_GO <- ggplot(data_join, aes(x = method, y = Description)) +

  geom_point(aes(size = setSize, color = p.adjust)) +
  facet_wrap(~factor(DE_type, levels = c("Up-regulated", "Down-regulated")), scales = "free", ncol = 2) +
  scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p.adjust)") +
  theme_bw() +
  # theme(
  #   #panel.spacing = unit(1, "lines"),
  #   #axis.text.x = element_text(size = 16, color = "black", angle = 30, hjust = 1),
  #   #axis.text.y = element_text(size = 16, color = "black"),
  #   #legend.key.size = unit(0.6, "cm"),
  #   #legend.text = element_text(size = 12, color = "black"),
  #   #legend.title = element_text(size = 14, color = "black"),
  #   #legend.spacing.y = unit(0.2, 'cm'),
  #   #strip.text = element_text(size = 18, face = "plain", color = "black"),
  #   #axis.text = element_text(size = 14, color = "black"),
  #   #axis.title = element_text(size = 16)
  # ) +
  labs(
    title = "",
    x = "",
    y = "Biological Process GO term",
    size = "Gene Count"
  )
plot_GO

ggsave("plot/enrichment_dotplot.png", dpi = 400, width = 18.0, height = 8.0, plot = plot_GO)
saveRDS(plot_GO, "plot/enrichment_dotplot.RDS")


### Categorize the DE genes based on their presence in core_enrichment as
### Biologically significant and Less biologically significant

library(stringr)

classify_genes <- function(de_genes, enrichment_data, core_col = "core_enrichment") {
  core_enrichment_genes <- enrichment_data %>%
    pull(!!sym(core_col)) %>%
    strsplit("/") %>%
    unlist() %>%
    unique()

  # Classify DE genes
  classified_genes <- de_genes %>%
    mutate(
      BiolSignificance = ifelse(gene %in% core_enrichment_genes,
                                      "Biologically significant",
                                      "Less Biologically significant")
    )

  return(classified_genes)
}


perform_classification <- function(de_genes, filtered_genes, method_name) {
  classified <- classify_genes(de_genes, filtered_genes)

  # Separate biologically significant and less significant genes
  biol_sign <- classified %>%
    dplyr::filter(BiolSignificance == "Biologically significant") %>%
    pull(gene)

  less_biol_sign <- classified %>%
    dplyr::filter(BiolSignificance == "Less Biologically significant") %>%
    pull(gene)

  classification_list <- list(
    method = method_name,
    biol_significant_genes = biol_sign,
    less_biol_significant_genes = less_biol_sign
  )

  return(classification_list)
}

all_results <- list()

methods <- list(
  devil = list(de_genes = rna_deg_devil$geneID, filtered = filtered_devil_nonRed$core_enrichment),
  glmGamPoi = list(de_genes = rna_deg_glm$geneID, filtered = filtered_glm_nonRed$core_enrichment),
  nebula = list(de_genes = rna_deg_nebula$geneID, filtered = filtered_nebula_nonRed$core_enrichment)
)


for (method in names(methods)) {
  de_genes <- data.frame(gene = methods[[method]]$de_genes)
  filtered <- data.frame(core_enrichment = methods[[method]]$filtered)
  result <- perform_classification(de_genes, filtered, method)
  all_results[[method]] <- result
}

saveRDS(all_results, file = "results/gsea_GO/gene_classification_allRes.RDS")



### Pathways categorization ###

gseGO_devil <- readRDS("results/gsea_GO/gseGO_devil.RDS")
gseGO_glm <- readRDS("results/gsea_GO/gseGO_glmGamPoi.RDS")
gseGO_nebula <- readRDS("results/gsea_GO/gseGO_nebula.RDS")

immune_system <- c("defense response", "acute-phase response", "positive regulation of cytokine production",
                   "inflammatory response", "hemopoiesis", "cellular response to cytokine stimulus", "leukocyte migration",
                   "cytokine production")

cellular_processes <- c("homophilic cell adhesion via plasma membrane adhesion molecules", "epithelial cell proliferation",
                        "regulation of leukocyte activation", "cell-cell adhesion", "leukocyte cell-cell adhesion")

cell_cycle_related_processes <- c("positive regulation of cell division", "regulation of G2/M transition of mitotic cell cycle",
                                  "regulation of cell division", "cell division", "regulation of cell cycle process")

muscle_system_process <- c("myofibril assembly", "actin-mediated cell contraction", "muscle cell development",
                           "muscle contraction", "actin filament-based movement")

signal_transduction <- c("cell surface pattern recognition receptor signaling pathway", "positive regulation of ERK1 and ERK2 cascade",
                         "immune response-activating cell surface receptor signaling pathway", "regulation of MAPK cascade",
                         "regulation of response to external stimulus", "positive regulation of developmental process", "positive regulation of multicellular organismal process",
                         "positive regulation of intracellular signal transduction", "positive regulation of gene expression",
                         "apoptotic signaling pathway")

metabolism <- c("tRNA metabolic process", "ncRNA metabolic process")

system_pathways <- c("skeletal system development", "skeletal system morphogenesis")

genetic_information_processing <- c("ncRNA processing",
                                    "nucleosome organization", 
                                    "heterochromatin formation")

unrelated_pathways <- c("cellular response to radiation", 
                        "animal organ morphogenesis", 
                        "embryo development",
                        "cellular response to gamma radiation", 
                        "regulation of cardiac muscle contraction by calcium ion signaling",
                        "regulation of neurotransmitter secretion", 
                        "cardiac atrium development", 
                        "retina morphogenesis in camera-type eye", 
                        "pancreas development",
                        "response to estrogen", 
                        "neural nucleus development", 
                        "neutral amino acid transport", 
                        "regulation of viral process",
                        "cellular response to tumor necrosis factor") 

pathways_categories <- list(immune_system, cellular_processes, cell_cycle_related_processes, muscle_system_process, signal_transduction,
                            metabolism, system_pathways, genetic_information_processing, unrelated_pathways)

names(pathways_categories) <- c("immune system", "cellular processes", "cell cycle related processes", "muscle system process", "signal transduction",
                                "metabolism", "system pathways", "genetic information processing", "unrelated pathways")

saveRDS(pathways_categories, file = "results/gsea_GO/pathway_categories.RDS")

