setwd("~/Desktop/dottorato/rdevil_project/de_analysis/nullpower")
rm(list = ls())
source("utils/glmGamPoi.R")
source("utils/devil.R")
source("utils/nebula.R")
require(magrittr)

data <- readRDS("test_data/pb.FALSE.bca.n.10.ct.1.fc.0.5.csv")

X <- data$count

meta <- data$meta

X <- X[1:400,]

glm.fit <- glmGamPoi::glm_gp(X, design=~1+tx_cell, col_data = meta, on_disk=FALSE, size_factors=FALSE)

design_matrix <- model.matrix(~1+tx_cell, data = meta)

batch_factor <- 4
count <- X %>% as.matrix()
devil.fit <- rdevil::fit_linear_model(
  input_matrix = count, model_matrix = design_matrix, group_matrix = NULL,
  variance = "Hessian",
  inference_method = "SVI",
  size_factors = FALSE,
  method_specific_args = list(
    optimizer_name = "ClippedAdam",
    steps = as.integer(1000),
    lr = 0.5,
    gamma_lr = 1e-6,
    cuda = TRUE,
    jit_compile = TRUE,
    batch_size = as.integer(dim(count)[2] / batch_factor),
    full_cov = TRUE,
    disp_loc = 3,
    gauss_loc = 5
    #prior_loc = 2,
    #theta_bounds = c(0, 1e5),
    #init_loc = .25,
    #threshold = 1e-9
  )
)

glm.fit$Beta[1,]
devil.fit$params$beta[,1]

plot(glm.fit$Beta[,1], devil.fit$params$beta[1,])
plot(glm.fit$Beta[,2], devil.fit$params$beta[2,])
plot(glm.fit$overdispersions, devil.fit$params$theta)

glm.res <- glmGamPoi::test_de(glm.fit, contrast = c(0,1))
devil.res <- rdevil::test_posterior_null(devil.fit, contrast = as.array(c(0,1)))

plot(glm.res$pval, devil.res$p_value)

max_idx <- as.integer(5000 * .05)

devil.null.p.values <- devil.res$p_value[(max_idx+1):nrow(devil.res)]
glm.null.p.values <- glm.res$pval[(max_idx+1):nrow(devil.res)]

plot(devil.null.p.values, glm.null.p.values)

gap::qqunif(devil.null.p.values)
gap::qqunif(glm.null.p.values)

