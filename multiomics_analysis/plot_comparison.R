rm(list = ls())

# Devil ####
res_devil_rna <- readRDS("results/devil/res_muscl_rna_age_devil.RDS")
res_devil_atac <- load("results/devil/res_atac_age_devil.Rdata")
res_devil_atac <- res_atac
rm(res_atac)

res_devil_rna$gene = res_devil_rna$name
res_devil_atac$gene = res_devil_atac$name

gg <- unique(c(res_devil_atac$gene, res_devil_rna$gene))

plot(res_devil_atac$lfc, res_devil_rna$lfc)
res_edgeR <- readRDS("results/edgeR_res.rds")
res_edgeR$adj_pval <- p.adjust(res_edgeR$`Pr(>|t|)`, "BH")
annot <- readRDS("results/grange_annot")

res_edgeR$gene <- annot$SYMBOL
res_edgeR$annot <- annot$annotation

d_join <- res_edgeR %>%
  na.omit() %>%
  as_tibble() %>%
  dplyr::left_join(res_devil, by='gene')

d_join %>%
  dplyr::filter(abs(lfc) > .5, abs(Estimate) > .5) %>%
  ggplot(mapping = aes(x=-Estimate, y=lfc)) +
  geom_point() +
  ggpubr::stat_cor(method="pearson") +
  labs(x = "scATC logFC", y="devil logFC")


# glmGamPoi ####
load("results/glm/res_rna_age_glm.Rdata")
res_rna_glm$gene = res_rna_glm$name

d_join <- res_edgeR %>%
  na.omit() %>%
  as_tibble() %>%
  dplyr::left_join(res_rna_glm, by='gene')

d_join %>%
  dplyr::filter(abs(lfc) > .5, abs(Estimate) > .5) %>%
  ggplot(mapping = aes(x=-Estimate, y=lfc)) +
  geom_point() +
  ggpubr::stat_cor(method="pearson") +
  labs(x = "scATC logFC", y="glmGamPoi logFC")

# nebula
load("results/nebula/res_rna_age_nebula.Rdata")
res_rna_neb$gene = res_rna_neb$name

d_join <- res_edgeR %>%
  na.omit() %>%
  as_tibble() %>%
  dplyr::left_join(res_rna_neb, by='gene')

d_join %>%
  dplyr::filter(abs(lfc) > .5, abs(Estimate) > .5) %>%
  ggplot(mapping = aes(x=-Estimate, y=lfc)) +
  geom_point() +
  ggpubr::stat_cor(method="pearson") +
  labs(x = "scATC logFC", y="nebula logFC")
