#' Plotting sim results
#'
#' long_results_list should have length four:
#' n = 500, 1000, 2000, 4000
#' each should then have length four:
#' p = 5, 10, 20, 40
#' each should then have length 16:
#' models (a, b, c, d) vs (a, b, c, d)
#'
#' plan: four plots, one for each N value
plot_sims <- function(all_results,
                      cut_y = 0.5) {
    plot(x = 0, y = 0,
         xlim = c(0, 3),
         ylim = c(0, 2.8),
         col = rgb(0, 0, 0, 0),
         xlab = "",
         ylab = "",
         axes = FALSE)

    par(xpd = TRUE)
    par(mar = c(1, 1, 1, 1))

    x_push <- 0.4
    y_push <- 0.6

    plot_results_block(all_results %>%
                       dplyr::filter(n_rows == 2000),
                       x_lim = c(0, 0.95),
                       y_lim = c(0, 0.95),
                       cut_y)
    plot_results_block(all_results %>%
                       dplyr::filter(n_rows == 4000),
                       x_lim = c(1.05, 2) + x_push,
                       y_lim = c(0, 0.95),
                       cut_y)
    plot_results_block(all_results %>%
                       dplyr::filter(n_rows == 500),
                       x_lim = c(0, 0.95),
                       y_lim = c(1.05, 2) + y_push,
                       cut_y)
    plot_results_block(all_results %>%
                       dplyr::filter(n_rows == 1000),
                       x_lim = c(1.05, 2) + x_push,
                       y_lim = c(1.05, 2) + y_push,
                       cut_y)

    ## add legends, overall title

    text(x = c(0.475, 1.525, 0.475, 1.525) + c(0, x_push, 0, x_push),
         y = c(0.95, 2, 2, 0.95) + c(0, y_push, y_push, 0) + 0.25,
         adj = c(0.5, 0.5),
         labels = paste("n = ", c(500, 1000, 2000, 4000)))

    y_seq <- seq(1.4, 2.3, length.out = 5)
    col_vec <- c(rgb(0, 0, 0, 0.5),
                 rgb(0, 0, 0, 1),
                 rgb(0.7, 0.3, 0.2, 1),
                 rgb(0.2, 0.8, 0.4),
                 rgb(0.3, 0.1, 0.9, 1))

    segments(x0 = rep(2 + x_push * 1.3, 5),
             x1 = rep(2 + x_push * 1.8, 5),
             y0 = y_seq,
             col = col_vec,
             lty = c(2, 2, 1, 1, 1))
    points(x = rep(2 + x_push * 1.55, 5),
           y = y_seq,
           pch = c(NA, 20, 20, 20, 20),,
           col = col_vec)
    text(x = rep(2 + x_push * 1.9, 5),
         y = y_seq,
         adj = c(0, 0.5),
         labels = c("naive CI bounds", "naive", "propensity",
                    "mahal", "Our Method"))

    rect(xleft = 2 + x_push * 1.535,
         xright = 2 + x_push * 1.565,
         ybottom = 0.8,
         ytop = 1.05,
         col = 1, border = NA)
    points(x = 2 + x_push * 1.55,
           y = 0.925,
           pch = 20,
           cex = 1.8)

    text(x = rep(2 + x_push * 1.9, 2),
         y = c(0.95, 0.88),
         adj = c(0, 0.5),
         cex = c(1, 0.8),
         labels = c("CI, ±2 se", "(bootstrapped)"))
}

plot_results_block <- function(n_results,
                               x_lim,
                               y_lim,
                               cut_y = 0.5,
                               rect_width = 0.006,
                               x_shift = 0.015) {
    rmse_boot <- function(vec) {
        sd(boot_func(vec, rmse_func, 500L))
    }

    summary_res <- n_results %>%
        group_by(p_cut) %>%
        summarise(rmse_naive = rmse_func(naive_est),
                  se_naive = rmse_boot(naive_est),
                  rmse_propensity = rmse_func(propensity_est),
                  se_propensity = rmse_boot(propensity_est),
                  rmse_mahal = rmse_func(mahal_est),
                  se_mahal = rmse_boot(mahal_est),
                  rmse_weighted = rmse_func(weighted_est),
                  se_weighted = rmse_boot(weighted_est)) %>%
        ungroup() %>%
        mutate(high_naive = pmin(cut_y * 1.6, rmse_naive + 2 * se_naive),
               low_naive = pmin(cut_y * 1.6, pmax(rmse_naive - 2 * se_naive, 0)),
               high_propensity = pmin(cut_y * 1.6, rmse_propensity + 2 * se_propensity),
               low_propensity = pmin(cut_y * 1.6, pmax(rmse_propensity - 2 * se_propensity, 0)),
               high_mahal = pmin(cut_y * 1.6, rmse_mahal + 2 * se_mahal),
               low_mahal = pmin(cut_y * 1.6, pmax(rmse_mahal - 2 * se_mahal, 0)),
               high_weighted = pmin(cut_y * 1.6, rmse_weighted + 2 * se_weighted),
               low_weighted = pmin(cut_y * 1.6, pmax(rmse_weighted - 2 * se_weighted, 0)))


    max_height <- min(max(max(summary_res[["high_propensity"]]),
                          max(summary_res[["high_mahal"]]),
                          max(summary_res[["high_weighted"]])) * 1.05,
                      cut_y)
    min_x <- min(summary_res[["p_cut"]])

    x_adj <- function(x) {
        ((x - min_x)/ (1 - min_x)) * (x_lim[2] - x_lim[1])+ x_lim[1]
    }
    y_adj <- function(y) {
        (y / max_height) * (y_lim[2] - y_lim[1])+ y_lim[1]
    }

    rect_width <- rect_width * (1 - min_x)
    x_shift <- x_shift * (1 - min_x)
    x_ax_adj <- 0.05 * (1 - min_x)
    y_ax_adj <- 0.02 * max_height

    ##------------------------------------

    ## axes
    segments(x0 = x_adj(min_x - x_ax_adj),
             x1 = x_adj(1),
             y0 = y_adj(0))
    x_tick_pre_seq <- unique(summary_res[["p_cut"]])
    x_tick_seq <- x_adj(x_tick_pre_seq)
    segments(x0 = x_tick_seq,
             y0 = y_adj(-y_ax_adj),
             y1 = y_adj(0))
    text(x = x_tick_seq,
         y = y_adj(-y_ax_adj * 2),
         labels = round(x_tick_pre_seq, 2),
         adj = c(0.5, 0.5),
         cex = 0.5)

    segments(x0 = x_adj(min_x - x_ax_adj),
             y0 = y_adj(0),
             y1 = y_adj(max_height))
    y_tick_pre_seq <- seq(from = 0,
                          to = max_height,
                          length.out = 5)
    y_tick_seq <- y_adj(y_tick_pre_seq)
    segments(y0 = y_tick_seq,
             x0 = x_adj(min_x - x_ax_adj * 2),
             x1 = x_adj(min_x - x_ax_adj))
    text(x = x_adj(min_x - x_ax_adj * 3),
         y = y_tick_seq,
         labels = round(y_tick_pre_seq, 2),
         adj = c(0.5, 0.5),
         cex = 0.5)

    par(xpd = FALSE)

    ##------------------------------------
    ## naive, may not even show!

    naive_mean <- summary_res[["rmse_naive"]][1]
    naive_high <- summary_res[["high_naive"]][1]
    naive_low <- summary_res[["low_naive"]][1]

    if (naive_mean / max_height < 1.1) {
        segments(x0 = x_lim[1],
                 x1 = x_lim[2],
                 y0 = y_adj(naive_mean),
                 lty = 2, col = 1, lwd = 1)
        segments(x0 = x_lim[1],
                 x1 = x_lim[2],
                 y0 = y_adj(naive_high),
                 lwd = 0.7,
                 lty = 2, col = rgb(0, 0, 0, 0.5))
        segments(x0 = x_lim[1],
                 x1 = x_lim[2],
                 y0 = y_adj(naive_low),
                 lwd = 0.7,
                 lty = 2, col = rgb(0, 0, 0, 0.5))
    }

    ##------------------------------------
    ## propensity

    points(x_adj(summary_res[["p_cut"]] - x_shift),
           y_adj(summary_res[["rmse_propensity"]]),
           type = "l", col = rgb(0.7, 0.3, 0.2, 1))
    points(x_adj(summary_res[["p_cut"]] - x_shift),
           y_adj(summary_res[["rmse_propensity"]]),
           pch = 20,
           col = rgb(0.7, 0.3, 0.2, 1))
    rect(xleft = x_adj(summary_res[["p_cut"]] - rect_width - x_shift),
         xright = x_adj(summary_res[["p_cut"]] + rect_width - x_shift),
         ybottom = y_adj(summary_res[["low_propensity"]]),
         ytop = y_adj(summary_res[["high_propensity"]]),
         col = rgb(0.7, 0.3, 0.2, 0.8),
         border = NA)

    ##------------------------------------
    ## mahal

    points(x_adj(summary_res[["p_cut"]] + x_shift),
           y_adj(summary_res[["rmse_mahal"]]),
           type = "l", col = rgb(0.2, 0.8, 0.4, 1))
    points(x_adj(summary_res[["p_cut"]] + x_shift),
           y_adj(summary_res[["rmse_mahal"]]),
           pch = 20,
           col = rgb(0.2, 0.8, 0.4, 1))
    rect(xleft = x_adj(summary_res[["p_cut"]] - rect_width + x_shift),
         xright = x_adj(summary_res[["p_cut"]] + rect_width + x_shift),
         ybottom = y_adj(summary_res[["low_mahal"]]),
         ytop = y_adj(summary_res[["high_mahal"]]),
         col = rgb(0.2, 0.8, 0.4, 0.8),
         border = NA)

    ##------------------------------------
    ## weighted

    points(x_adj(summary_res[["p_cut"]]),
           y_adj(summary_res[["rmse_weighted"]]),
           type = "l", col = rgb(0.3, 0.1, 0.9, 1))
    points(x_adj(summary_res[["p_cut"]]),
           y_adj(summary_res[["rmse_weighted"]]),
           pch = 20,
           col = rgb(0.3, 0.1, 0.9, 1))
    rect(xleft = x_adj(summary_res[["p_cut"]] - rect_width),
         xright = x_adj(summary_res[["p_cut"]] + rect_width),
         ybottom = y_adj(summary_res[["low_weighted"]]),
         ytop = y_adj(summary_res[["high_weighted"]]),
         col = rgb(0.3, 0.1, 0.9, 0.8),
         border = NA)

    par(xpd = TRUE)

    rect(xleft = x_adj(0.4),
         xright = x_adj(1.1),
         ybottom = y_adj(cut_y * 1.05),
         ytop = 50,
         border = NA, col = rgb(1, 1, 1, 1))
}

boot_func <- function(vec, func, resamples = 100L) {
    rep_list <- lapply(1L:resamples, function(j){
        func(sample(vec, replace = TRUE))
    })
    unlist(rep_list)
}

rmse_func <- function(vec) {
    sqrt(mean((vec - 1)^2))
}
