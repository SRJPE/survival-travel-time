# ============================================================
# 52_summarise_lorenzen_sims.R
# Summarises results/reps.csv from 51_run_lorenzen_sims.R.
# Writes results/summary_by_scenario.csv and two figures.
# ============================================================

library(tidyverse)

res <- read.csv("results/reps.csv")

summ <- res |>
  group_by(scenario) |>
  summarise(
    reps = n(),
    recoveries_median = median(recoveries),
    c_true = first(c_true), c_est_median = median(c_med),
    c_cover90 = mean(c_true >= c_q05 & c_true <= c_q95),
    delta_true = first(delta_true), delta_est_median = median(delta_med),
    delta_cover90 = mean(delta_true >= delta_q05 & delta_true <= delta_q95),
    delta_flagged = mean(delta_q05 > 0 | delta_q95 < 0),       # 90% interval excludes 0
    delta_width = median(delta_q95 - delta_q05),
    cwt_surv_log_bias = median(surv_log_bias),                  # log(estimate / truth), groups averaged
    cwt_surv_cover90 = mean(surv_cover90),
    capture_log_bias = median(prec_log_bias),
    reps_with_divergences = sum(divergences > 0),
    minutes_median = median(minutes),
    .groups = "drop") |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write.csv(summ, "results/summary_by_scenario.csv", row.names = FALSE)
print(as.data.frame(summ))

lab <- c(S1_match = "S1 match (c = 1, delta = 0)", S2_cwt_lower = "S2 CWT hazard 1.5x (delta = 0.4)",
         S3_c_half = "S3 true c = 0.5", S4_per_km = "S4 mortality per km")
col <- c(S1_match = "#0072B2", S2_cwt_lower = "#D55E00", S3_c_half = "#009E73", S4_per_km = "#E69F00")

plot_par <- function(par, truth_col, ylab, file) {
  d <- res |> transmute(scenario, rep, med = .data[[paste0(par, "_med")]],
                        lo = .data[[paste0(par, "_q05")]], hi = .data[[paste0(par, "_q95")]],
                        truth = .data[[truth_col]])
  p <- ggplot(d, aes(rep, med, ymin = lo, ymax = hi, colour = scenario)) +
    geom_hline(aes(yintercept = truth), linetype = "dashed", colour = "grey30") +
    geom_pointrange(size = 0.25) +
    facet_wrap(~ scenario, labeller = as_labeller(lab), ncol = 2) +
    scale_colour_manual(values = col, guide = "none") +
    labs(x = "Replicate", y = ylab, title = paste("One-way Lorenzen model:", ylab),
         subtitle = "Posterior median and 90% interval; dashed = truth (S4 has no true delta: it shows the per-day misfit)") +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
  ggsave(file, p, width = 9, height = 6, dpi = 200, bg = "white")
  p
}
plot_par("delta", "delta_true", "CWT offset delta (log hazard)", "results/fig_delta.png")
plot_par("c", "c_true", "Lorenzen exponent c", "results/fig_c.png")
