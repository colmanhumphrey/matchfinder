##' Solving GLS when variance is known up to a scale factor
##'
##' If we know that \eqn{y \sim \mathcal{N}(X \beta, \sigma^2 \Sigma)}
##' where \eqn{\Sigma} is known but \eqn{\beta} and \eqn{\sigma^2} aren't,
##' then we can solve using nearly straight-forward generalised least squares.
##' We return a list of results: the coefficients (\eqn{\beta}); the standard
##' error for those coefficients; the whole covariance matrix for \eqn{\beta}
##' and finally our estimate for \eqn{\sigma^2}.
##' @param y_vector Outcome vector (in future can potentially be matrix)
##' @param x_mat Matrix of explanatory variables. We **will not include the
##'   intercept for you** so add a column of 1s if you want this. You can
##'   also supply just a single vector and we'll upgrade it to a matrix
##'   (but again we will not add an intercept). Needs to match length
##'   of \code{y_vector}
##' @param rel_var_mat \eqn{\Sigma} in the equations above: the variance of the
##'   error up to some (unknown) constant factor.
##' @return List with the following elements:
##' * \code{beta_gls}: the coefficients, \eqn{\hat{\beta}}
##' * \code{beta_gls_stderr}: the standard errors for the coefficients
##' * \code{var_beta_gls}: the covariance matrix for \eqn{\beta}
##' * \code{sigma_sq_est}: our estimate for \eqn{\sigma^2}.
##' If \code{x_mat} has column names, they should hopefully pass through
##' to the resulting vectors/matrices.
##' @author Colman Humphrey
##'
##' @export
gls_known_relative_variance <- function(y_vector, x_mat, rel_var_mat) {
    if (is.matrix(y_vector)) {
        stop("matrix `y_vector` not yet implemented")
    }

    if (!is.matrix(rel_var_mat)) {
        stop("`rel_var_mat` needs to be a matrix")
    }

    if (nrow(rel_var_mat) != length(y_vector)) {
        err_msg <- paste0("`rel_var_mat` has ", nrow(rel_var_mat), " rows but ",
                          "y_vector has length ", length(y_vector))
        stop(err_msg)
    }

    if (nrow(rel_var_mat) != ncol(rel_var_mat)) {
        err_msg <-
            paste0("`rel_var_mat` should be square, instead has dimension (",
                   nrow(rel_var_mat), ", ", ncol(rel_var_mat), ")")
        stop(err_msg)
    }

    if (is.matrix(x_mat)) {
        if (nrow(x_mat) != length(y_vector)) {
            err_msg <- paste0("`x_mat` has ", nrow(x_mat), " rows but ",
                              "y_vector has length ", length(y_vector))
            stop(err_msg)
        }
    } else {
        if (length(x_mat) != length(y_vector)) {
            err_msg <-
                paste0("input vector (`x_mat`) not same length as `y_vector` ",
                       "- should `x_mat` be a matrix?")
            stop(err_msg)
        }
        x_mat <- matrix(x_mat, ncol = 1L)
    }

    gls_action(y_vector, x_mat, rel_var_mat)
}

##' Does the actual work of \code{gls_known_relative_variance}
##'
##' @keywords internal
gls_action <- function(y_vector, x_mat, rel_var_mat) {
    sig_inv_x <- solve(rel_var_mat, x_mat)
    x_sig_inv_x <- t(x_mat) %*% sig_inv_x

    beta_gls <-  solve(x_sig_inv_x, t(sig_inv_x) %*% y_vector)
    y_resid <- y_vector - x_mat %*% beta_gls

    ## this can be solved with e.g.
    ## \code{mahalanobis(c(y_resid), rep(0, 500), rel_var_mat)}
    ## (so, mahalanobis distance of `y` residuals to zero
    ##  with scale `rel_var_mat`) but on testing this was faster,
    ## because we don't actually need to invert the matrix
    sigma_sq_est <- c(t(y_resid) %*% solve(rel_var_mat, y_resid)) /
        (length(y_vector) - ncol(x_mat))

    var_beta_gls_pre <- solve(x_sig_inv_x)
    var_beta_gls <- var_beta_gls_pre * sigma_sq_est

    beta_stderr <- sqrt(diag(var_beta_gls))

    list(
        beta_gls = beta_gls[, 1, drop = TRUE],
        beta_gls_stderr = beta_stderr,
        var_beta_gls = var_beta_gls,
        sigma_sq_est = sigma_sq_est
    )
}

##' Not using this yet, but provides a useful alternative
##' to \code{gls_action}: using R's in-built \code{lm}
##' but is slower
##'
##' @keywords internal
gls_svd_lm <- function(y_vector, x_mat, rel_var_mat) {
    var_svd <- svd(rel_var_mat)
    sqrt_inv_diag <- diag(1 / sqrt(var_svd$d))
    half_inv <- sqrt_inv_diag %*% t(var_svd$u)

    y_star <- c(half_inv %*% y_vector)
    x_star <- half_inv %*% x_mat

    lm(y_star ~ x_star + 0)
}
