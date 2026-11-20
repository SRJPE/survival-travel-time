# ============================================================
# 58_summarise_gamma_realsize.R (as 56, for 57)
# Summarises results/reps_gamma_realsize.csv from 55_run_lorenzen_gamma.R.
# Writes results/summary_gamma_realsize.csv and results/fig_gamma_realsize.png (Okabe-Ito).
# True delta = 0 and true c = 1 throughout.
# ============================================================
library(tidyverse)

res <- read.csv("results/reps_gamma_realsize.csv") |> distinct(truth, fit, rep, .keep_all = TRUE)
summ <- res |>
  group_by(truth, gamma_true, fit) |>
  summarise(reps = n(),
            gamma_est = median(gamma_med), gamma_cover90 = mean(gamma_true >= gamma_q05 & gamma_true <= gamma_q95),
            delta_est = median(delta_med), delta_flagged = mean(delta_q05 > 0 | delta_q95 < 0),
            c_est = median(c_med), c_cover90 = mean(1 >= c_q05 & 1 <= c_q95),
            surv_log_bias = median(surv_log_bias), surv_log_bias_small = median(surv_log_bias_small),
            surv_log_bias_large = median(surv_log_bias_large), surv_cover90 = mean(surv_cover90),
            capture_log_bias = median(prec_log_bias), reps_with_divergences = sum(divergences > 0),
            .groups = "drop") |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write.csv(summ, "results/summary_gamma_realsize.csv", row.names = FALSE)
print(as.data.frame(summ))

lab_t <- c(gamma_1 = "Truth: per day (gamma = 1)", `gamma_0.5` = "Truth: slow fish accumulate less (gamma = 0.5)", gamma_0 = "Truth: per trip (gamma = 0)")
lab_f <- c(per_day = "Fit: per day (gamma fixed at 1)", gamma = "Fit: gamma estimated")
long <- res |>
  select(truth, gamma_true, fit, rep, starts_with("gamma_"), starts_with("delta_"), starts_with("c_")) |>
  pivot_longer(c(gamma_med, delta_med, c_med), names_to = "par", values_to = "med") |>
  mutate(par = sub("_med", "", par),
         lo = case_when(par == "gamma" ~ gamma_q05, par == "delta" ~ delta_q05, TRUE ~ c_q05),
         hi = case_when(par == "gamma" ~ gamma_q95, par == "delta" ~ delta_q95, TRUE ~ c_q95),
         true = case_when(par == "gamma" ~ gamma_true, par == "delta" ~ 0, TRUE ~ 1),
         par = factor(par, levels = c("delta", "c", "gamma"), labels = c("CWT offset delta", "Size exponent c", "Exposure gamma")),
         cell = paste(lab_t[truth], lab_f[fit], sep = "\n")) |>
  filter(!(fit == "per_day" & par == "Exposure gamma"))
p <- ggplot(long, aes(rep, med, ymin = lo, ymax = hi, colour = fit)) +
  geom_hline(aes(yintercept = true), linetype = "dashed", colour = "grey30") +
  geom_pointrange(size = 0.2) +
  facet_grid(par ~ cell, scales = "free_y") +
  scale_colour_manual(values = c(per_day = "#D55E00", gamma = "#0072B2"), guide = "none") +
  labs(x = "Replicate", y = NULL, title = "Exposure simulations, realistic sizes and real number of CWT groups (360)",
       subtitle = "Posterior median and 90% interval; dashed = truth") +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave("results/fig_gamma_realsize.png", p, width = 15, height = 7, dpi = 200, bg = "white")
message("Wrote results/summary_gamma_realsize.csv and results/fig_gamma_realsize.png")
