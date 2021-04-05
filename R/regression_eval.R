##' Evaluate match using regression
##'
##' Generates a list with the estimated treatment effect,
##' and the standard error for that effect. If bipartite
##' (i.e. no tolerance used) just regresses differences in
##' \code{y_vector} against the intercept, basically difference in means.
##' If a tolerance vector is used, the differences are regressed
##' against the tolerance differences.
##' In both cases, if controls are repeated, we build up the appropriate
##' correlation matrix and use Generalised Linear Regression to solve.
##' @param match_list Usual match list
##' @param y_vector Usual outcome vector
##' @param tolerance_list Usual tolerance list
##' @return Returns a list:
##' * \code{estimate}: the estimated treatment effect. If we're in bipartite
##'     land, then this is the usual difference in means, else
##'     it's the effect of a "unit" change in tolerance
##' * \code{standard_error}: the standard deviation of the above coefficient
##' @author Colman Humphrey
##'
##' @export
regression_eval <- function(match_list,
                            y_vector,
                            tolerance_list = gen_tolerance_list()) {
    tol_validation <- tolerance_check(match_list, tolerance_list)
    if (tol_validation[["error"]]) {
        stop(tol_validation[["message"]])
    }

    tol_diffs <- if (is.null(tolerance_list)) {
                     rep(1.0, length(match_list[["treat_index"]]))
                 } else {
                     tol_vec <- tolerance_list[["tolerance_vec"]]
                     stopifnot(length(tol_vec) == length(y_vector))
                     tol_vec[match_list[["treat_index"]]] -
                         tol_vec[match_list[["control_index"]]]
                 }

    y_diffs <- y_vector[match_list[["treat_index"]]] -
        y_vector[match_list[["control_index"]]]

    if (anyNA(y_diffs)) {
        err_msg <- if (anyNA(y_vector)) {
                       paste0("`y_vector` has ", sum(is.na(y_vector)),
                              " missing values")
                   } else {
                       max_match_index <- max(
                           match_list[["treat_index"]],
                           match_list[["control_index"]])
                       if (max_match_index > length(y_vector)) {
                           paste0("`match_list` has max index ",
                                  max_match_index,
                                  " but `y_vector` is only of length ",
                                  length(y_vector))
                       } else {
                           ## shouldn't really happen
                           paste0(
                               "`y_diffs` has missing values, not sure why; ",
                               "check `match_list`"
                           )
                       }
                   }
        stop(err_msg)
    }

    ## TODO: treatment duplication, cross-over duplication
    all_unique_controls <- !any(duplicated(match_list[["control_index"]]))

    if (all_unique_controls) {
        diff_reg <- lm(y_diffs ~ tol_diffs + 0)
        coef_res <- summary(diff_reg)[["coefficients"]]

        list(
            estimate = coef_res[1L, "Estimate"],
            standard_error = coef_res[1L, "Std. Error"]
        )
    } else {
        rel_var_mat <- diag(rep(1.0, length(match_list[["treat_index"]])))
        half_list <- pair_indexing(build_index_vector(
            match_list[["control_index"]]))
        rel_var_mat[half_list] <- 0.5

        gls_result <- gls_known_relative_variance(
            y_vector = y_diffs,
            x_mat = matrix(tol_diffs, ncol = 1L),
            rel_var_mat = rel_var_mat
        )

        list(
            estimate = gls_result[["beta_gls"]],
            standard_error = gls_result[["beta_gls_stderr"]]
        )
    }
}
