library(readr)
library(dplyr)
library(matchfinder)

sims <- read_csv("~/r-packages/data/matchfinder_sims_output.csv")

##------------------------------------
## while we shouldn't of course delete any rows,
## propensity has a really nasty value of -250 showing up, which I don't really
## care to count against it
## for reference, the largest weighted_est is like 30

plot_num <- function(num_delete) {
    bad_ids <- sims %>%
        group_by(id) %>%
        summarise(propensity_est = max(abs(propensity_est)),
                  mahal_est = max(abs(mahal_est)),
                  weighted_est = max(abs(weighted_est))) %>%
        filter(rank(-abs(propensity_est), ties.method = "first") <= num_delete |
               rank(-abs(mahal_est), ties.method = "first") <= num_delete |
               rank(-abs(weighted_est), ties.method = "first") <= num_delete) %>%
        pull(id) %>%
        unique()

    clean_sims <- sims %>%
        filter(!(id %in% bad_ids))

    full_plot_frame <- generate_p_cut_frame(
        clean_sims,
        p_cut_vals = seq(from = 0.4, to = 1, by = 0.1))

    plot_list <- generate_plot_list(full_plot_frame,
                                    rmse_func = rmse_from_one_func)

    file_name <- paste0("~/Downloads/sim_res/p", num_delete, ".pdf")

    pdf(file_name, height = 15, width = 15)
    plot_sims(plot_list)
    dev.off()
}

for (j in (0:10) * 20) {
    plot_num(j)
}

##------------------------------------

num_delete <- 10

bad_ids <- sims %>%
    filter(rank(-abs(propensity_est)) < num_delete |
           rank(-abs(weighted_est)) < num_delete |
           rank(-abs(mahal_est)) < num_delete) %>%
    pull(id) %>%
    unique()

clean_sims <- sims %>%
    filter(!(id %in% bad_ids))

full_plot_frame <- generate_p_cut_frame(
    clean_sims,
    p_cut_vals = seq(from = 0.4, to = 1, by = 0.1))

plot_list <- generate_plot_list(full_plot_frame,
                                rmse_func = rmse_from_one_func)

## just until we finish the sims!
plot_list[["500"]] <- plot_list[["1000"]]

pdf("~/Downloads/sim_res/p10.pdf", height = 20, width = 20)
plot_sims(plot_list)
dev.off()

##------------------------------------

num_delete <- 20

bad_ids <- sims %>%
    filter(rank(-abs(propensity_est)) < num_delete |
           rank(-abs(weighted_est)) < num_delete |
           rank(-abs(mahal_est)) < num_delete) %>%
    pull(id) %>%
    unique()

clean_sims <- sims %>%
    filter(!(id %in% bad_ids))

full_plot_frame <- generate_p_cut_frame(
    clean_sims,
    p_cut_vals = seq(from = 0.4, to = 1, by = 0.1))

plot_list <- generate_plot_list(full_plot_frame,
                                rmse_func = rmse_from_one_func)

## just until we finish the sims!
plot_list[["500"]] <- plot_list[["1000"]]

pdf("~/Downloads/sim_res/p20.pdf", height = 20, width = 20)
plot_sims(plot_list)
dev.off()

##------------------------------------

pdf("~/Downloads/test.pdf", width = 20, height = 20)
matchfinder::plot_sims(full_p_cut, cut_y = 0.5)
dev.off()

##------------------------------------
