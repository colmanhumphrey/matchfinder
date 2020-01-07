library(dplyr)
library(matchfinder)
library(magrittr)
source("sim_funcs.R")
source("reshape_results.R")
source("plot_sim_results.R")

##------------------------------------

## last_time <- Sys.time()
## for(treat_model in c("a", "b", "c", "d")) {
##     for (mu_model in c("a", "b", "c", "d")) {
##         print(paste0(treat_model, mu_model))
##         print(Sys.time())
##         print(difftime(Sys.time(), last_time))
##         last_time <- Sys.time()
##         big_res <- parallel_sim(
##             treat_model = treat_model,
##             mu_model = mu_model,
##             n_sink_gen = n_sink_func_func(),
##             n_rows = 2000,
##             n_cols = 10,
##             num_weight_vectors = 100,
##             match_method = "with_replacement",
##             num_cores = parallel::detectCores(),
##             iterations = 100,
##             x_norm = TRUE,
##             silent = FALSE)
##         saveRDS(big_res,
##                 file = paste0("~/Downloads/sim_res/loop_", treat_model,
##                               "_", mu_model, ".rds"))
##     }
## }
## print(Sys.time())
## print(difftime(Sys.time(), last_time))

p_cut <- seq(from = 1, to = 0.5, length.out = 6)

stupid_name_fix <- function(listo) {
    lapply(listo, function(x) {
        if ("permutation_results" %in% names(x)) {
            return(list(
                naive_est = x$naive_est,
                propensity_results = x$permutation_results,
                mahal_results = x$mahal_results,
                weighted_results = x$weighted_results
            ))
        }
        x
    })
}

res_list <- list()

for(treat_model in c("a", "b", "c", "d")) {
    print(treat_model)
    res_temp <- list()
    for (mu_model in c("a", "b", "c", "d")) {
        file_name <- paste0("~/Downloads/sim_res/loop_", treat_model,
                            "_", mu_model, ".rds")
        rel_res <- stupid_name_fix(readRDS(file_name))
        res_temp[[mu_model]] <- reshape_parallel_sim(rel_res,
                                                     treat_model,
                                                     mu_model,
                                                     n_rows = 2000,
                                                     n_cols = 10,
                                                     p_cut = p_cut)
    }
    res_list[[treat_model]] <- dplyr::bind_rows(res_temp)
}

total_results <- dplyr::bind_rows(res_list)


temp_total <- lapply(c(500, 1000, 2000, 4000), function(x){
    temp_res <- total_results
    temp_res$n_rows <- x
    temp_res
}) %>%
    dplyr::bind_rows()
## temp_small <- temp_total %>% filter(mu_model != 'd')
temp_small <- temp_total
pdf("~/Downloads/temp.pdf", width = 12, height = 9)
plot_sims(all_results = temp_small, cut_y = 0.3)
dev.off()


##------------------------------------

boot_func <- function(vec, func, resamples = 100L) {
    rep_list <- lapply(1L:resamples, function(j){
        func(sample(vec, replace = TRUE))
    })
    unlist(rep_list)
}

rmse_func <- function(vec) {
    sqrt(mean((vec - 1)^2))
}

rmse_boot <- function(vec) {
    sd(boot_func(vec, rmse_func, 500L))
}

summary_res <- total_results %>%
    group_by(p_cut) %>%
    summarise(rmse_naive = rmse_func(naive_est),
              se_naive = rmse_boot(naive_est),
              rmse_propensity = rmse_func(propensity_est),
              se_propensity = rmse_boot(propensity_est),
              rmse_mahal = rmse_func(mahal_est),
              se_mahal = rmse_boot(mahal_est),
              rmse_weighted = rmse_func(weighted_est),
              se_weighted = rmse_boot(weighted_est))


summary_res <- total_results %>%
    group_by(p_cut, treat_model, mu_model) %>%
    summarise(rmse_naive = rmse_func(naive_est),
              se_naive = rmse_boot(naive_est),
              rmse_propensity = rmse_func(propensity_est),
              se_propensity = rmse_boot(propensity_est),
              rmse_mahal = rmse_func(mahal_est),
              se_mahal = rmse_boot(mahal_est),
              rmse_weighted = rmse_func(weighted_est),
              se_weighted = rmse_boot(weighted_est)) %>%
    ungroup()

temp <- summary_res %>%
    group_by(treat_model, mu_model) %>%
    summarise(se_naive_mean = mean(se_naive),
              se_propensity_mean = mean(se_propensity),
              se_mahal_mean = mean(se_mahal),
              se_weighted_mean = mean(se_weighted)) %>%
    ungroup()

temp %>%
    group_by(treat_model) %>%
    summarise(naive = mean(se_naive_mean),
              propensity = mean(se_propensity_mean),
              mahal = mean(se_mahal_mean),
              weighted = mean(se_weighted_mean))
temp %>%
    group_by(mu_model) %>%
    summarise(naive = mean(se_naive_mean),
              propensity = mean(se_propensity_mean),
              mahal = mean(se_mahal_mean),
              weighted = mean(se_weighted_mean))


pdf("~/Downloads/temp_rmse.pdf", height = 30, width = 30)
par(mfrow = c(4, 4))

for(t_model in c("a", "b", "c", "d")) {
    print(t_model)
    for (m_model in c("a", "b", "c", "d")) {
        temp_res <- summary_res %>%
            filter(treat_model == t_model,
                   mu_model == m_model)
        plot(p_cut, rep(0, length(p_cut)), type = 'l', ylim = c(0, 0.6),
             col = rgb(0, 0, 0, 0))
        text(x = 0.7, y = 0.4, labels = t_model)
        text(x = 0.7, y = 0.35, labels = m_model)
        points(temp_res$p_cut, temp_res$rmse_weighted, type = 'l', col = 1)
        points(temp_res$p_cut, temp_res$rmse_propensity, type = 'l', col = 2)
        points(temp_res$p_cut, temp_res$rmse_mahal, type = 'l', col = 3)
    }
}
dev.off()

pdf("~/Downloads/temp_se.pdf", height = 30, width = 30)
par(mfrow = c(4, 4))

for(t_model in c("a", "b", "c", "d")) {
    print(t_model)
    for (m_model in c("a", "b", "c", "d")) {
        temp_res <- summary_res %>%
            filter(treat_model == t_model,
                   mu_model == m_model)
        plot(p_cut, rep(0, length(p_cut)), type = 'l', ylim = c(0, 0.6),
             col = rgb(0, 0, 0, 0))
        text(x = 0.7, y = 0.4, labels = t_model)
        text(x = 0.7, y = 0.35, labels = m_model)
        points(temp_res$p_cut, temp_res$se_weighted, type = 'l', col = 1)
        points(temp_res$p_cut, temp_res$se_propensity, type = 'l', col = 2)
        points(temp_res$p_cut, temp_res$se_mahal, type = 'l', col = 3)
    }
}
dev.off()
