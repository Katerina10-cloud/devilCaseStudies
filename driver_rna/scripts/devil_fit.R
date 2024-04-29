###----------------------------------------------------------###
### Devil testing ###
###----------------------------------------------------------###

#devtools::install_github("caravagnalab/devil")

library(devil)
library(tidyverse)

#Create design matrix
design_matrix <- model.matrix(~ cell_clusters, data = metadata)

### Parameters inference ###

devil_fit_retina <- devil::fit_devil(input_matrix = rna_counts,
                             design_matrix = design_matrix,
                             overdispersion = TRUE,
                             offset=0,
                             size_factors=TRUE,
                             verbose=TRUE,
                             max_iter=500,
                             eps=1e-4,
                             parallel = TRUE)

### Statistical test ###
contrast_vector <- c(0,1)

stat_test_res <- devil::test_de(devil.fit = devil_fit_retina,
                                contrast = contrast_vector,
                                pval_adjust_method = "BH",
                                max_lfc = 10,
                                clusters = NULL)
