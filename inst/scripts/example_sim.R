## library(matchfinder)
devtools::load_all()

##------------------------------------

## simulate some data

n_cols <- 8L
half_rows <- 500L
rows <- half_rows * 2L

a_mat_pre <- matrix(runif(n_cols * n_cols, 0, 1),
                    ncol = n_cols)
a_mat <- cov2cor(t(a_mat_pre) %*% a_mat_pre)

z_mat <- matrix(rnorm(n_cols * rows),
                ncol = n_cols)
x_mat <- z_mat %*% t(a_mat)

rank_cols <- c(4, 5)

cov_x <- covariance_with_ranks(x_mat, rank_cols)
x_ranked <- ranked_x(x_mat, rank_cols)

## cov(x_mat) - a_mat %*% t(a_mat)


prior_weights <- c(1, 1, 2, 2, 1, 1, 1, 3)
hier_list <- list(
    list(index = c(1, 2, 3), weight = 2),
    list(index = c(4, 5, 6), weight = 3),
    list(index = c(7, 8), weight = 2)
)
min_weights <- rep(1/(8 * 3), times = 8L)

treat_vec <- (1L:rows) %in% fixed_sample((1L:rows), size = half_rows - 10L)

coef_vec <- rnorm(n_cols, sd = 0.2)
y_vector <- x_mat %*% coef_vec + treat_vec * 0.4

##------------------------------------

random_weights <- generate_random_weights(prior_weights = prior_weights,
                                          number_vectors = 50,
                                          minimum_weights = min_weights,
                                          hierarchical_list = hier_list)

match_lists <- all_weighted_matches(x_mat = x_ranked,
                                    cov_x,
                                    weight_list = random_weights,
                                    n_sinks = 50L,
                                    treat_vec = treat_vec)
