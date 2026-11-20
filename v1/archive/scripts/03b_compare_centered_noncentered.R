# ============================================================
# 03b_compare_centered_noncentered.R
#
# Side by side: Flora's 25 Sep fit (centered, rstan) and the
# non-centered refit from 03_Run_CJS_Model_noncentered.R.
# Same data, likelihood and priors; only the parameterization
# (and the software that ran it) differs.
#
# Writes results/tables/compare_centered_noncentered.csv and
# results/figures/compare_centered_noncentered.png.
# ============================================================

library(rstan)
library(here)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

fits <- c(centered     = here("results", "fit_CovIndCont_MaxFlow_FL.Rdata"),
          noncentered  = here("results", "fit_CovIndCont_noncentered_MaxFlow_FL.Rdata"))

# population-level parameters (not the per-group or per-fish ones)
pars <- c("S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz",
          "S_bTrib", "T_bTrib", "S_bCovT", "TT_bCovT",
          "RE_sd", "RE_sdT", "TTRE_sd", "TTRE_sdT", "Pro_sd", "Pro_sdT",
          "muPb", "sdPb")

summarise_fit <- function(file, label) {
  e <- new.env()
  load(file, envir = e)
  f <- e$fit
  sp <- get_sampler_params(f, inc_warmup = FALSE)
  s <- summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |>
    as.data.frame() |>
    rownames_to_column("parameter") |>
    transmute(fit = label, parameter, mean, q2.5 = `2.5%`, median = `50%`, q97.5 = `97.5%`,
              rhat = Rhat, ess = n_eff)
  sampler <- tibble(fit = label,
                    draws = sum(map_int(sp, nrow)),
                    divergences = sum(map_dbl(sp, \(x) sum(x[, "divergent__"]))),
                    max_treedepth_seen = max(map_dbl(sp, \(x) max(x[, "treedepth__"]))),
                    hours_slowest_chain = max(rowSums(get_elapsed_time(f))) / 3600,
                    worst_rhat = max(s$rhat, na.rm = TRUE),
                    lowest_ess = min(s$ess, na.rm = TRUE))
  rm(e); gc()
  list(params = s, sampler = sampler)
}

res <- imap(fits, summarise_fit)

sampler_tbl <- map(res, "sampler") |> list_rbind()
print(as.data.frame(sampler_tbl), digits = 3)

compare_tbl <- map(res, "params") |> list_rbind() |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(compare_tbl, here("results", "tables", "compare_centered_noncentered.csv"), row.names = FALSE)

# wide view in the console: estimate, interval, Rhat, ESS for both fits
compare_tbl |>
  select(fit, parameter, median, q2.5, q97.5, rhat, ess) |>
  tidyr::pivot_wider(names_from = fit, values_from = c(median, q2.5, q97.5, rhat, ess)) |>
  as.data.frame() |>
  print()

# Figure: median and 95% interval, one row per parameter.
# Okabe-Ito blue / vermillion plus different shapes, so colour is not needed.
p <- ggplot(compare_tbl, aes(y = parameter, x = median, xmin = q2.5, xmax = q97.5,
                             colour = fit, shape = fit)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.3) +
  scale_colour_manual(values = c(centered = "#D55E00", noncentered = "#0072B2"),
                      labels = c("Flora's fit (centered, 25 Sep)", "Non-centered refit")) +
  scale_shape_manual(values = c(centered = 17, noncentered = 16),
                     labels = c("Flora's fit (centered, 25 Sep)", "Non-centered refit")) +
  scale_y_discrete(limits = rev) +
  labs(x = "Posterior median and 95% interval (model scale)", y = NULL,
       colour = NULL, shape = NULL,
       title = "CovIndCont (max flow + FL): centered vs non-centered random effects") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("results", "figures", "compare_centered_noncentered.png"), p,
       width = 8, height = 11, dpi = 200, bg = "white")
