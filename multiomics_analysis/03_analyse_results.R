### Results downstream analysis ###

rm(list = ls())

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/devilCaseStudies/multiomics_analysis")

pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork", 'ggVennDiagram', 'stringr')
sapply(pkgs, require, character.only = TRUE)

rna_devil <- "results/MuscleRNA/devil_rna.RDS"
rna_devil <- readRDS(rna_devil) %>% dplyr::rename(geneID=name)

rna_glm <- "results/MuscleRNA/glmGamPoi_rna.RDS"
rna_glm <- readRDS(rna_glm) %>% dplyr::rename(geneID=name)

rna_nebula <- "results/MuscleRNA/nebula_rna.RDS"
rna_nebula <- readRDS(rna_nebula) %>% dplyr::rename(geneID=name) %>% dplyr::mutate(lfc = lfc / log(2))


# Gene selection based on LFC & pvalue cutoff #
lfc_cut <- 1.0
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

saveRDS(venn_plot, "plot/venn_plot_v2.rds")


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
    size = 4.0,
    #fontface = "bold",
    color = "black",
    fill = "white",
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.5,
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.4, "lines"),
    min.segment.length = 1
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(floor(min(rna_join$lfc)),
                                  ceiling(max(rna_join$lfc)), by = 2)) +
  facet_wrap(~factor(method, levels = c("devil", "glmGamPoi", "nebula")), nrow = 1, scales = "free") +
  labs(x = expression(Log[2] ~ FC),
       y = expression(-log[10] ~ Pvalue),
       col = "DE type") +
  #ggplot2::theme(legend.position = 'right',
                 #legend.text = element_text(size = 16, color = "black"),
                 #legend.title = element_text(size = 18, color = "black"),  
                 #strip.text = element_text(size = 20, face = "plain", color = "black"),
                 #axis.text = element_text(size = 18, color = "black"),
                 #axis.title = element_text(size = 20, color = "black"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
p_volcanos

ggsave("plot/volcanos_v2.png", dpi = 400, width = 16.0, height = 5.0, plot = p_volcanos)
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
    #seed = 1234
  )

  return(gseGO@result %>% as.data.frame())
}


rna_deg_glm$adj_pval[rna_deg_glm$adj_pval == 0] <- min(rna_deg_glm$adj_pval[rna_deg_glm$adj_pval != 0])

gseGO_devil <- enrichmentGO(rna_deg_devil)
gseGO_glm <- enrichmentGO(rna_deg_glm)
gseGO_nebula <- enrichmentGO(rna_deg_nebula)


# Enrichment results preprocessing #

# Remove redundant terms 

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
  
  redundant_data <- filtered_data %>%
    filter(!!sym(desc_col) %in% redundant_terms)

  return(list(
    non_redundant_data = non_redundant_data,
    redundant_data = redundant_data
  ))
}


filtered_devil <- remove_redundant_terms(gseGO_devil,
                                         enrichment_col = "enrichmentScore",
                                         core_col = "core_enrichment",
                                         desc_col = "Description",
                                         threshold = 0.5)

redundant_devil <- as.data.frame(filtered_devil[["redundant_data"]])
nonRed_devil <- as.data.frame(filtered_devil[["non_redundant_data"]])


filtered_glm <- remove_redundant_terms(gseGO_glm,
                                       enrichment_col = "enrichmentScore",
                                       core_col = "core_enrichment",
                                       desc_col = "Description",
                                       threshold = 0.5)

redundant_glm <- as.data.frame(filtered_glm[["redundant_data"]])
nonRed_glm <- as.data.frame(filtered_glm[["non_redundant_data"]])


filtered_nebula <- remove_redundant_terms(gseGO_nebula,
                                          enrichment_col = "enrichmentScore",
                                          core_col = "core_enrichment",
                                          desc_col = "Description",
                                          threshold = 0.5)

redundant_nebula <- as.data.frame(filtered_nebula[["redundant_data"]])
nonRed_nebula <- as.data.frame(filtered_nebula[["non_redundant_data"]])


saveRDS(filtered_devil, file = "results/gsea_GO/gseGO_devil_list.RDS")
saveRDS(filtered_glm, file = "results/gsea_GO/gseGO_glmGamPoi_list.RDS")
saveRDS(filtered_nebula, file = "results/gsea_GO/gseGO_nebula_list.RDS")


# Enrichment results | Pathways selection #

filtered_devil <- readRDS("results/gsea_GO/gseGO_devil_list.RDS")
filtered_glm <- readRDS("results/gsea_GO/gseGO_glmGamPoi_list.RDS")
filtered_nebula <- readRDS("results/gsea_GO/gseGO_nebula_list.RDS")

filtered_devil_nonRed <- filtered_devil[["non_redundant_data"]]
filtered_glm_nonRed <- filtered_glm[["non_redundant_data"]]
filtered_nebula_nonRed <- filtered_nebula[["non_redundant_data"]]

filtered_devil_nonRed$DE_type <- ifelse(filtered_devil_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")
filtered_glm_nonRed$DE_type <- ifelse(filtered_glm_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")
filtered_nebula_nonRed$DE_type <- ifelse(filtered_nebula_nonRed$enrichmentScore > 0, "Up-regulated", "Down-regulated")

filtered_devil_nonRed <- filtered_devil_nonRed %>%
  dplyr::mutate(method = "devil")

filtered_glm_nonRed <- filtered_glm_nonRed %>%
  dplyr::mutate(method = "glmGamPoi")

filtered_nebula_nonRed <- filtered_nebula_nonRed %>%
  dplyr::mutate(method = "nebula")



# Select not biologically specific pathways 

terms_notSpecific_devil <- c("embryo development")

terms_notSpecific_nebula <- c("response to tumor necrosis factor", 
                              "positive regulation of gene expression",
                              "embryonic morphogenesis", 
                              "tube development", 
                              "gene expression",
                              "regulation of biological process")


terms_notSpecific_glm <- c("response to gamma radiation",
                           "artery morphogenesis",
                           "positive regulation of tumor necrosis factor superfamily cytokine production",
                           "regulation of viral process",
                           "cognition",
                           "inner ear development",
                           "cell fate commitment",
                           "response to steroid hormone",
                           "viral process",
                           "response to tumor necrosis factor",
                           "reproductive structure development",
                           "kidney development",
                           "positive regulation of gene expression",
                           "negative regulation of gene expression",
                           "embryonic morphogenesis",
                           "nervous system development",
                           "central nervous system development",
                           "cardiac muscle hypertrophy")

not_biol_specific_terms <- list(terms_notSpecific_devil, terms_notSpecific_glm, terms_notSpecific_nebula)
names(not_biol_specific_terms) <- c("devil", "glmGamPoi", "nebula")

saveRDS(not_biol_specific_terms, file = "results/gsea_GO/notBiolSpecific_terms.RDS")


# Check glmGamPoi private genes enrichment #

glm_private_genes <- setdiff(rna_deg_glm$geneID, union(rna_deg_devil$geneID, rna_deg_nebula$geneID))

cat("Number of private genes in GLM:", length(glm_private_genes), "\n")

glm_private_data <- rna_deg_glm %>%
  dplyr::filter(geneID %in% glm_private_genes)

glm_private_data$adj_pval[glm_private_data$adj_pval == 0] <- min(glm_private_data$adj_pval[glm_private_data$adj_pval != 0])

gseGO_glm_private <- enrichmentGO(glm_private_data)

saveRDS(gseGO_glm_private, file = "results/gsea_GO/gseGO_glm_private.RDS")


# Divide pathways into BP categories #

biological_processes <- list(
  immune_response_and_defense = c(
    "innate immune response activating cell surface receptor signaling pathway",
    "immune response-regulating cell surface receptor signaling pathway",
    "defense response",
    "inflammatory response",
    "positive regulation of immune system process",
    "regulation of response to external stimulus",
    "B cell activation",
    "cellular response to cytokine stimulus",
    "cytokine production",
    "negative regulation of immune system process",
    "positive regulation of response to external stimulus",
    "interleukin-8 production"
    
  ),
  cellular_processes_and_regulation = c(
    "regulation of cytokinesis",
    "positive regulation of cell cycle process",
    "positive regulation of cytokine production",
    "ncRNA processing",
    "cellular component assembly involved in morphogenesis",
    "regulation of protein secretion",
    "negative regulation of cell cycle phase transition",
    "extrinsic apoptotic signaling pathway",
    "extracellular matrix organization",
    "epithelial cell proliferation",
    "response to hydrogen peroxide",
    "regulation of receptor-mediated endocytosis",
    "regulation of RNA splicing",
    "activation of protein kinase activity"
  ),
  
  signaling = c("DNA damage checkpoint signaling",
                "hippo signaling"),
  
  cellular_adhesion_and_movement = c(
    "homophilic cell adhesion via plasma membrane adhesion molecules",
    "myeloid leukocyte migration",
    "regulation of cell adhesion",
    "positive regulation of multicellular organismal process",
    "leukocyte migration",
    "leukocyte cell-cell adhesion",
    "positive regulation of cell migration",
    "regulation of substrate adhesion-dependent cell spreading",
    
  ),
  
  metabolism = c("organic cyclic compound catabolic process",
                 "proteolysis",
                 "protein modification by small protein conjugation or removal",
                 "negative regulation of metabolic process",
                 "macromolecule modification"),
  
  transport = c("intracellular transport",
                "iron ion transport"
  ),
  
  development_and_differentiation = c(
    "osteoblast differentiation",
    "cartilage development"
   
  ),
  muscle_function = c(
    "actin filament-based movement",
    "muscle contraction",
    "actin filament-based process",
    "muscle system process"
  )
)

go_terms_levels <- unname(unlist(biological_processes))

data_join$Biological_process = lapply(data_join$Description, function(desc) {
  true_bp <- lapply(names(biological_processes), function(bp) {
    if (desc %in% biological_processes[[bp]]) {
      return(bp)
    }
  }) %>% unlist()
  true_bp
})


plot_GO = data_join %>%
  dplyr::mutate(Description = factor(Description, levels = go_terms_levels)) %>%
  dplyr::mutate(DE_type = factor(DE_type, levels = c("Up-regulated", "Down-regulated"))) %>%
  ggplot(aes(x = method, y = Description, size = setSize, color = p.adjust)) +
  geom_point() +
  facet_grid(~DE_type, space = "free", scales = "free") +
  #facet_wrap(~factor(DE_type, levels = c("Up-regulated", "Down-regulated")), scales = "free_y", ncol = 2, shrink = F) +
  scale_color_gradient(low = "cornflowerblue", high = "coral", name = "p.adjust)") +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Biological Process GO term",
    size = "Gene Count"
  )
plot_GO

ggsave("plot/enrichment_dotplot.png", dpi = 400, width = 18.0, height = 8.0, plot = plot_GO)
saveRDS(plot_GO, "plot/enrichment_dotplot.RDS")

plot_data <- as.data.frame(plot_GO[["data"]])


### Summary of results across methods ###

gse_devil <- readRDS("results/gsea_GO/gseGO_devil_list.RDS")
gse_glm <- readRDS("results/gsea_GO/gseGO_glmGamPoi_list.RDS")
gse_nebula <- readRDS("results/gsea_GO/gseGO_nebula_list.RDS")

notSpecific_terms <- readRDS("results/gsea_GO/notBiolSpecific_terms.RDS")

calculate_counts <- function(redundant_data, non_redundant_data, non_specific_terms) {
  redundant_terms_count <- length(redundant_data)
  non_redundant_terms_count <- length(non_redundant_data)
  non_specific_terms_count <- length(non_specific_terms)
  
  return(c(redundant_terms_count, non_redundant_terms_count, non_specific_terms_count))
}

devil_counts <- calculate_counts(
  gse_devil[["non_redundant_data"]][["Description"]],
  gse_devil[["redundant_data"]][["Description"]],
  notSpecific_terms[["devil"]]
)

glm_counts <- calculate_counts(
  gse_glm[["non_redundant_data"]][["Description"]],
  gse_glm[["redundant_data"]][["Description"]],
  notSpecific_terms[["glmGamPoi"]]
)

nebula_counts <- calculate_counts(
  gse_nebula[["non_redundant_data"]][["Description"]],
  gse_nebula[["redundant_data"]][["Description"]],
  notSpecific_terms[["nebula"]]
)

summary_table <- data.frame(
  method = c("devil", "glmGamPoi", "nebula"),
  redundant_terms = c(devil_counts[1], glm_counts[1], nebula_counts[1]),
  non_redundant_terms = c(devil_counts[2], glm_counts[2], nebula_counts[2]),
  non_specific_terms = c(devil_counts[3], glm_counts[3], nebula_counts[3])
)


### Categorize the DE genes based on their presence in core_enrichment as
### Biologically significant and Less biologically significant

#library(stringr)

#classify_genes <- function(de_genes, enrichment_data, core_col = "core_enrichment") {
  #core_enrichment_genes <- enrichment_data %>%
    #pull(!!sym(core_col)) %>%
    #strsplit("/") %>%
    #unlist() %>%
    #unique()

  # Classify DE genes
  #classified_genes <- de_genes %>%
    #mutate(
      #BiolSignificance = ifelse(gene %in% core_enrichment_genes,
                                      #"Biologically significant",
                                      #"Less Biologically significant")
    #)

  #return(classified_genes)
#}


#perform_classification <- function(de_genes, filtered_genes, method_name) {
  #classified <- classify_genes(de_genes, filtered_genes)

  # Separate biologically significant and less significant genes
  #biol_sign <- classified %>%
    #dplyr::filter(BiolSignificance == "Biologically significant") %>%
    #pull(gene)

  #less_biol_sign <- classified %>%
    #dplyr::filter(BiolSignificance == "Less Biologically significant") %>%
    #pull(gene)

  #classification_list <- list(
    #method = method_name,
    #biol_significant_genes = biol_sign,
    #less_biol_significant_genes = less_biol_sign
  #)

  #return(classification_list)
#}

#all_results <- list()

#methods <- list(
  #devil = list(de_genes = rna_deg_devil$geneID, filtered = filtered_devil_nonRed$core_enrichment),
  #glmGamPoi = list(de_genes = rna_deg_glm$geneID, filtered = filtered_glm_nonRed$core_enrichment),
  #nebula = list(de_genes = rna_deg_nebula$geneID, filtered = filtered_nebula_nonRed$core_enrichment)
#)


#for (method in names(methods)) {
  #de_genes <- data.frame(gene = methods[[method]]$de_genes)
  #filtered <- data.frame(core_enrichment = methods[[method]]$filtered)
  #result <- perform_classification(de_genes, filtered, method)
  #all_results[[method]] <- result
#}

#saveRDS(all_results, file = "results/gsea_GO/gene_classification_allRes.RDS")





