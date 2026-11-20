# 11_summarise_sim_recovery.R
# Summarise whatever replicates 10_sim_recovery_cjs.R has finished so far.
# Safe to run while 10 is still going.
#
# Writes to the run folder:
#   summary_by_parameter.csv  coverage, bias, interval width, Rhat by model x parameter
#   sampler_summary.csv       time, divergences, treedepth by model
#   coverage_plot.png         95% interval coverage by parameter (target 0.95)
#   error_plot.png            estimate minus truth by parameter

library(tidyverse)
library(here)
source(here("scripts", "sims", "sim_functions.R"))

run <- "run1" # or "smoke"
out_dir <- here("results", "sims", "cjs_recovery", run)

res <- read_sim_results(out_dir)
stopifnot(nrow(res$estimates) > 0)

okabe_ito <- c(production = "#0072B2", noncentered = "#D55E00")

param_order <- res$estimates |> distinct(parameter) |> pull(parameter)

by_param <- res$estimates |>
  group_by(model, parameter) |>
  summarise(
    n_reps = n(),
    truth = first(truth),
    mean_estimate = mean(mean),
    bias = mean(error),
    rel_bias = if_else(abs(first(truth)) > 1e-8, bias / abs(first(truth)), NA_real_),
    coverage = mean(covered),
    mean_ci_width = mean(ci_width),
    median_rhat = median(Rhat, na.rm = TRUE),
    share_rhat_over_1.01 = mean(Rhat > 1.01, na.rm = TRUE),
    median_n_eff = median(n_eff, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(parameter = factor(parameter, levels = param_order)) |>
  arrange(parameter, model)

sampler <- res$sampler |>
  group_by(model) |>
  summarise(
    fits = n(),
    mean_minutes = mean(minutes),
    total_hours = sum(minutes) / 60,
    fits_with_divergences = sum(divergences > 0),
    mean_divergences = mean(divergences),
    fits_with_treedepth_hits = sum(treedepth_hits > 0),
    .groups = "drop"
  )

write_csv(by_param, file.path(out_dir, "summary_by_parameter.csv"))
write_csv(sampler, file.path(out_dir, "sampler_summary.csv"))
print(sampler)
print(by_param |> select(model, parameter, n_reps, truth, bias, coverage, mean_ci_width,
                         share_rhat_over_1.01), n = Inf)

n_reps_label <- paste0(n_distinct(res$estimates$rep), " replicates")

coverage_plot <- ggplot(by_param, aes(x = coverage, y = fct_rev(parameter),
                                      colour = model, shape = model)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = okabe_ito) +
  scale_shape_manual(values = c(production = 16, noncentered = 17)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Share of 95% intervals containing the true value", y = NULL,
       colour = NULL, shape = NULL,
       title = "Simulation recovery on the real design: coverage",
       subtitle = paste0("Dashed line = 0.95 target; ", n_reps_label)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

error_plot <- res$estimates |>
  mutate(parameter = factor(parameter, levels = param_order)) |>
  ggplot(aes(x = error, y = fct_rev(parameter), colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = okabe_ito) +
  scale_shape_manual(values = c(production = 16, noncentered = 17)) +
  labs(x = "Posterior mean minus truth", y = NULL, colour = NULL, shape = NULL,
       title = "Simulation recovery on the real design: error",
       subtitle = paste0("One point per replicate; ", n_reps_label)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "coverage_plot.png"), coverage_plot, width = 8, height = 9, dpi = 200)
ggsave(file.path(out_dir, "error_plot.png"), error_plot, width = 8, height = 9, dpi = 200)
message("Wrote summaries and plots to ", out_dir)
