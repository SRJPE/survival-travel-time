# ============================================================
# 03b_compare_three_fits.R
#
# Three fits side by side:
#   her_25Sep      : Flora's fit, centered, travel times MISALIGNED (pre-A1 data)
#   centered       : her same model on the current, aligned data (03c)
#   noncentered    : non-centered model on the current data (03 noncentered)
# her_25Sep vs centered     = effect of aligning the travel times
# centered vs noncentered   = effect of non-centering (should be none
#                             on estimates; only convergence)
#
# Writes results/tables/compare_three_fits.csv and
# results/figures/compare_three_fits.png.
# ============================================================

library(rstan)
library(here)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

fits <- c(her_25Sep   = here("results", "fit_CovIndCont_MaxFlow_FL.Rdata"),
          centered    = here("results", "fit_CovIndCont_centered_alignedTT_MaxFlow_FL.Rdata"),
          noncentered = here("results", "fit_CovIndCont_noncentered_MaxFlow_FL.Rdata"))
fit_labels <- c(her_25Sep   = "Flora's 25 Sep fit (misaligned travel times)",
                centered    = "Her model, aligned data",
                noncentered = "Non-centered, aligned data")

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
  mutate(across(where(is.numeric), \(x) round(x, 3)),
         fit = factor(fit, levels = names(fits)))
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(compare_tbl, here("results", "tables", "compare_three_fits.csv"), row.names = FALSE)

# wide view in the console: estimate, interval, Rhat, ESS for both fits
compare_tbl |>
  select(fit, parameter, median, q2.5, q97.5, rhat, ess) |>
  tidyr::pivot_wider(names_from = fit, values_from = c(median, q2.5, q97.5, rhat, ess)) |>
  as.data.frame() |>
  print()

# Figure: median and 95% interval, one row per parameter.
# Okabe-Ito colours plus different shapes, so colour is not needed.
p <- ggplot(compare_tbl, aes(y = parameter, x = median, xmin = q2.5, xmax = q97.5,
                             colour = fit, shape = fit)) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) +
  scale_colour_manual(values = c(her_25Sep = "#D55E00", centered = "#E69F00", noncentered = "#0072B2"),
                      labels = fit_labels) +
  scale_shape_manual(values = c(her_25Sep = 17, centered = 15, noncentered = 16),
                     labels = fit_labels) +
  scale_y_discrete(limits = rev) +
  labs(x = "Posterior median and 95% interval (model scale)", y = NULL,
       colour = NULL, shape = NULL,
       title = "CovIndCont (max flow + FL): three fits") +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("results", "figures", "compare_three_fits.png"), p,
       width = 8, height = 11, dpi = 200, bg = "white")
