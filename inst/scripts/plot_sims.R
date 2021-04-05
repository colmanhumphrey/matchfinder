library(readr)
library(dplyr)
library(matchfinder)

sims <- read_csv(paste0("~/r-packages/data/matchfinder/simulations/",
                        "matchfinder_sims_output.csv"))

full_plot_frame <- generate_p_cut_frame(
    clean_sims,
    p_cut_vals = seq(from = 0.4, to = 1, by = 0.1))

plot_list <- generate_plot_list(full_plot_frame,
                                rmse_func = rmse_from_one_func)

file_name <- paste0("~/Downloads/sim_res/p", num_delete, ".pdf")

pdf(file_name, height = 15, width = 15)
plot_sims(plot_list)
dev.off()

##------------------------------------
