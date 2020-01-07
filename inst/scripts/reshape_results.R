reshape_parallel_sim <- function(parallel_result,
                                 treat_model,
                                 mu_model,
                                 n_rows,
                                 n_cols,
                                 p_cut) {

    do.call(rbind, lapply(p_cut, function(p_cut_val) {
        do.call(rbind, lapply(parallel_result, function(par_res) {
            p_briers <- unlist(lapply(
                par_res[["weighted_results"]], function(x) {
                    x[["permutation_brier"]]
                }))
            if (all(p_briers > p_cut_val)) {
                given_cut_ind <- which.min(p_briers)
            } else {
                given_cut_ind <- which(p_briers <= p_cut_val)[1]
            }

            data.frame(
                treat_model = treat_model,
                mu_model = mu_model,
                p_cut = p_cut_val,
                p_brier =
                    par_res$weighted_results[[given_cut_ind]]$permutation_brier,
                raw_brier =
                    par_res$weighted_results[[given_cut_ind]]$raw_brier,
                n_rows = n_rows,
                n_cols = n_cols,
                n_sinks = par_res$weighted_results[[given_cut_ind]]$n_sinks,
                naive_est = par_res$naive_est,
                propensity_est =
                    par_res$propensity_results[[given_cut_ind]]$est,
                mahal_est = par_res$mahal_results[[given_cut_ind]]$est,
                weighted_est = par_res$weighted_results[[given_cut_ind]]$est
            )
        }))
    }))
}
