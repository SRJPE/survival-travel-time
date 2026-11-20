# ============================================================
# scripts/refit/07_fit_cwt_gamma.R
# (was 07_Run_CWT_gamma.R, run standalone on the laptop 9 Oct 2026; now uses the
#  project folders instead of setwd(), and fits G1 only by default)
#
# One-way CWT model with an exposure exponent gamma (stan/CWT_oneway_gamma.stan):
#   cumulative hazard = M_r * (70 / L)^c * exp(delta + h) * T^gamma
# gamma = 1: mortality per day (as 06); gamma < 1: slow fish accumulate less;
# gamma = 0: per trip. Same telemetry anchor, CWT data and capture as 06.
#   G1: gamma + delta  (the handoff model: c feeds lorenzen_pars.csv.gz, step 13)
#   G0: gamma, no delta (comparison only; set runs_to_fit below)
# 9 Oct 2026 result (G1): gamma 0.017 (0.001 to 0.08), delta +0.13, c 0.46 (0.13 to 0.80).
#
# Inputs : results/fit_CWT_oneway_Hmed.Rdata (06: base_data, cwt_groups, H, starting values)
#          stan/CWT_oneway_gamma.stan
# Outputs: results/fit_CWT_gamma_G1.Rdata (and _G0 if fitted),
#          results/tables/cwt_gamma_parameters.csv, cwt_gamma_loo.csv, cwt_gamma_by_year.csv,
#          results/figures/fig_cwt_gamma_surv_vs_tt.png
# 4 chains x 1000/1000, about 20 to 40 min per model.
# ============================================================

library(tidyverse)
library(here)
library(rstan)
library(cmdstanr)
library(loo)

runs_to_fit <- c("G1")          # add "G0" for the comparison
nchains <- 4; nwarmup <- 1000; nsamp <- 1000; adapt_delta <- 0.95; seed <- 123
gamma_prior <- c(mu = 1, sd = 0.5)

ow_file <- here("results", "fit_CWT_oneway_Hmed.Rdata")
stopifnot("Run 06_fit_cwt_oneway.R first" = file.exists(ow_file))
walk(here("results", c("tables", "figures")), \(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))

ow <- new.env(); load(ow_file, envir = ow)   # fit, base_data, cwt_groups, H, ...
pmed <- \(f, p) apply(as.matrix(rstan::extract(f, p)[[1]]), 2, median)

# starting values: the one-way fit's posterior medians
start <- set_names(c("T_bReach", "TT_bCov", "T_bSz", "z_TT", "sd_cwtTTre", "Pro_sd",
                     "z_pcap", "c", "delta", "z_h", "sd_h"), identity) |>
  map(\(p) pmed(ow$fit, p))
make_inits <- function() {
  set.seed(seed)
  map(seq_len(nchains), \(ch) {
    x <- imap(start, \(v, nm) if (nm %in% c("sd_cwtTTre", "Pro_sd", "sd_h", "c")) v * exp(rnorm(length(v), 0, 0.05))
                              else v + rnorm(length(v), 0, 0.05))
    x$c <- pmin(pmax(x$c, 0.05), 2.9)
    c(x, list(gamma_raw = 1 + rnorm(1, 0, 0.05)))
  })
}

fix_csv <- function(f) {
  x <- readLines(f)
  com <- startsWith(x, "#")
  x[com] <- gsub("= true", "= 1", gsub("= false", "= 0", x[com]))
  out <- sub("\\.csv$", "_rstan.csv", f)
  writeLines(x, out)
  out
}
mod <- cmdstan_model(here("stan", "CWT_oneway_gamma.stan"))
runs <- list(G1 = list(UseOffset = 1L), G0 = list(UseOffset = 0L))[runs_to_fit]
fit_run <- function(sw, run) {
  dat <- modifyList(c(ow$base_data, list(H_tel = ow$H)),
                    c(sw, list(UseGamma = 1L, gamma_mu = gamma_prior[["mu"]], gamma_sd = gamma_prior[["sd"]])))
  cs_dir <- here("results", paste0("cmdstan_gamma_", run))
  dir.create(cs_dir, showWarnings = FALSE)
  fc <- mod$sample(data = dat, init = make_inits(), chains = nchains, parallel_chains = nchains,
                   iter_warmup = nwarmup, iter_sampling = nsamp, seed = seed,
                   adapt_delta = adapt_delta, refresh = 250, output_dir = cs_dir)
  message("Run ", run, ":"); print(fc$diagnostic_summary())
  fit <- rstan::read_stan_csv(map_chr(fc$output_files(), fix_csv))
  save(fit, file = here("results", paste0("fit_CWT_gamma_", run, ".Rdata")))
  fit
}
fits <- c(imap(runs, fit_run), list(oneway = ow$fit))

# ---- parameters -----------------------------------------------------------------------
pars <- c("gamma_raw", "delta", "c", "sd_h", "M_r", "T_bReach", "T_bSz", "Pro_sd")
par_tbl <- imap(fits, \(f, run) {
  have <- intersect(pars, f@model_pars)
  summary(f, pars = have, probs = c(0.025, 0.5, 0.975))$summary |>
    as.data.frame() |> rownames_to_column("parameter") |>
    transmute(run, parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`, rhat = Rhat, ess = n_eff)
}) |> list_rbind() |>
  mutate(parameter = recode(parameter, gamma_raw = "gamma"),
         across(where(is.numeric), \(x) signif(x, 3))) |>
  filter(!(run == "G0" & parameter == "delta"))   # delta is unused (prior only) in G0
write.csv(par_tbl, here("results", "tables", "cwt_gamma_parameters.csv"), row.names = FALSE)
par_tbl |> select(run, parameter, median, q2.5, q97.5, rhat) |>
  pivot_wider(names_from = run, values_from = c(median, q2.5, q97.5, rhat)) |>
  as.data.frame() |> print()

# ---- LOO (binomial recoveries per group) -----------------------------------------------
loos <- map(fits, \(f) {
  ll <- extract_log_lik(f, "log_lik", merge_chains = FALSE)
  loo(ll, r_eff = relative_eff(exp(ll)))
})
loo_tbl <- as.data.frame(loo_compare(loos))
if (!"model" %in% names(loo_tbl)) loo_tbl <- rownames_to_column(loo_tbl, "model")
loo_tbl <- loo_tbl |>
  mutate(looic = map_dbl(model, \(m) loos[[m]]$estimates["looic", "Estimate"]),
         pareto_k_over_0.7 = map_int(model, \(m) sum(pareto_k_values(loos[[m]]) > 0.7)))
write.csv(loo_tbl, here("results", "tables", "cwt_gamma_loo.csv"), row.names = FALSE)
print(loo_tbl)

# ---- survival by year and against travel time ------------------------------------------
grp <- imap(fits, \(f, run) tibble(run, grp = seq_len(ow$base_data$Ncwtgrp),
                                   surv = pmed(f, "surv_cwt"), tt = pmed(f, "pTT"))) |>
  list_rbind() |>
  left_join(ow$cwt_groups |> select(grp, year, avg_length), by = "grp")
by_year <- grp |> group_by(run, year) |> summarise(mean_surv = mean(surv), .groups = "drop") |>
  pivot_wider(names_from = run, values_from = mean_surv)
write.csv(by_year, here("results", "tables", "cwt_gamma_by_year.csv"), row.names = FALSE)
print(as.data.frame(by_year), digits = 3)

run_lab <- c(oneway = "One-way (gamma = 1)", G1 = "G1 gamma + delta", G0 = "G0 gamma only")
p <- ggplot(grp, aes(tt, surv, colour = run, shape = run)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(oneway = "#0072B2", G1 = "#D55E00", G0 = "#009E73"), labels = run_lab, name = NULL) +
  scale_shape_manual(values = c(oneway = 16, G1 = 17, G0 = 15), labels = run_lab, name = NULL) +
  labs(x = "Median travel time to KL (days, log scale)", y = "CWT survival to KL (log scale)",
       title = "CWT survival against travel time: per day (gamma = 1) vs estimated exposure",
       subtitle = "Each point = one CWT release group, posterior median") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "fig_cwt_gamma_surv_vs_tt.png"), p, width = 8, height = 5.5, dpi = 200, bg = "white")
message("Done. Next: export step 13 (13_export_lorenzen_pars.R).")
