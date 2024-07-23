remotes::install_github("caravagnalab/rdevil", force = TRUE)
library("rdevil")

counts <- rnbinom(n=10, mu=5, size = 1/0.7)

fit <- rdevil::fit_linear_model(
	input_matrix = X,
	model_matrix = matrix(1, nrow = ncol(counts))
)
