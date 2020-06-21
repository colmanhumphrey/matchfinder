#' Generates a matrix where each row is generated
#' with covariance equal to \code{cov_mat}.
#'
#' @param cov_mat The covariance desired.
#' @param n_rows How many rows to produce.
#' @return A matrix such that the covariance would tend to
#'   \code{cov_mat} as \code{n_rows} grows.
#' @author Colman Humphrey
#'
#' @export
x_from_cov <- function(cov_mat,
                       n_rows) {
    ## recall: in R, chol gives U s.t. U' U = X, not U U'
    chol_cov <- chol(cov_mat)

    n_cols <- nrow(cov_mat)

    base_normal_mat <- matrix(rnorm(n_rows * n_cols),
        nrow = n_rows,
        ncol = n_cols
    )

    base_normal_mat %*% chol_cov
}


#' Shifts the linear component until the mean of the
#' expit is the target mean.
#'
#' @param target_mean Desired target mean.
#' @param linear_vector Basically \eqn{X \beta} (no shift).
#' @return \code{expit(linear_vector + a)} for some \code{a} such that
#'   \code{mean(expit(linear_vector + a))} is close to \code{target_mean}.
#' @author Colman Humphrey
#'
#' @export
target_mean_expit <- function(target_mean,
                              linear_vector) {
    mono_func <- function(alpha_value) {
        mean(expit(linear_vector + alpha_value))
    }

    alpha_value <- binary_search(
        target_value = target_mean,
        monotone_function = mono_func
    )

    expit(linear_vector + alpha_value)
}


#' An example function that generates treatment probabilities
#' given an input matrix.
#'
#' @param x_mat Numeric matrix.
#' @return A vector of probabilities, with length equal to \code{nrow(x_mat)}.
#' @author Colman Humphrey
#'
#' @export
example_treat_prob_generator <- function(x_mat) {
    flat_1 <- (x_mat[, 1] - mean(x_mat[, 1])) / sd(x_mat[, 1])
    numer <- (flat_1^3 - 2 * flat_1^2) / 10

    if (ncol(x_mat) > 1) {
        numer <- numer + sign(x_mat[, 2]) / 2 +
            x_mat %*% rnorm(
                n = ncol(x_mat),
                mean = 0, sd = 0.1
            )
    }

    expit(numer)
}


#' An example function that generates a mean vector
#' given an input matrix.
#'
#' @param x_mat Numeric matrix.
#' @return A vector of means, with length equal to \code{nrow(x_mat)}.
#' @author Colman Humphrey
#'
#' @export
example_mean_generator <- function(x_mat) {
    b_vec <- rnorm(
        n = ncol(x_mat),
        mean = 0,
        sd = 0.3
    )

    linear_part <- x_mat %*% b_vec

    square_part <- (x_mat[, 1] - mean(x_mat[, 1]))^2
    cross_part <- 0
    if (ncol(x_mat) > 1) {
        square_part <- square_part - (x_mat[, 2] - mean(x_mat[, 2]))^2
        cross_part <- (x_mat[, 1] - mean(x_mat[, 1])) *
            (x_mat[, ncol(x_mat)] - mean(x_mat[, ncol(x_mat)]))
    }

    linear_part + square_part + cross_part
}


#' A default function that generates an input data matrix.
#'
#' First generates a random covariance, then generates
#' normal data with that covariance (actually correlation).
#' @param n_rows How many rows to produce.
#' @param n_cols How many columns to produce.
#' @return A matrix of data
#' @author Colman Humphrey
#'
#' @export
default_x_generator <- function(n_rows,
                                n_cols) {
    sig_mat_pre <- matrix(
        stats::runif((n_cols - 1) * n_cols, 0, 0.1),
        n_cols, n_cols - 1
    )
    sig_mat_cov <- sig_mat_pre %*% t(sig_mat_pre) +
        diag(x = stats::runif(n_cols, 0, 0.3))
    sig_mat_cor <- stats::cov2cor(sig_mat_cov)

    x_from_cov(sig_mat_cor,
        n_rows = n_rows
    )
}


#' Default function to generate mean-zero noise. Very simple.
#'
#' Exported just so all functions have a default in
#' \code{generate_simulation_input}.
#' @param n_rows How many rows to produce.
#' @return Returns just mean 0 variance 1 normal noise.
#' @author Colman Humphrey
#'
#' @export
default_error_generator <- function(n_rows) {
    rnorm(n = n_rows)
}
