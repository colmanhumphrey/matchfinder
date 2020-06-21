#' computes the covariance of the input matrix, adjusting for ranks if given
#'
#' adjusts to ranks if given; fixes ties; rescales etc;
#' returns the covariance of x_mat. If no rank_cols are given, this
#' function is identical to just calling cov(x_mat)
#' @param x_mat matrix of variables (numeric)
#' @param rank_cols names of columns to be converted to ranks before analysis
#' @return covariance matrix of x_mat (potentially adjusted for ranks)
#'
#' @export
covariance_with_ranks <- function(x_mat,
                                  rank_cols = NULL) {
    if (!is.matrix(x_mat)) {
        stop("x_mat must be a matrix")
    }

    ## convert to ranks
    x_mat <- ranked_x(
        x_mat,
        rank_cols
    )

    cov_x <- cov(x_mat)

    if (is.null(rank_cols)) {
        return(cov_x)
    }

    rank_cols <- rank_integer_index(rank_cols, x_mat)

    ## fix the rank ties variances
    var_untied <- var((1:nrow(x_mat)) / nrow(x_mat))
    sd_ratio <- sqrt(var_untied / diag(cov_x))

    ## don't fix the ones that aren't ranked
    sd_ratio[!(1L:ncol(x_mat) %in% rank_cols)] <- 1

    diag(sd_ratio) %*% cov_x %*% diag(sd_ratio)
}
