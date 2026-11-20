# ============================================================
# 03h_diagnostics_production.R
# Convergence diagnostics for the production telemetry fit (03g), for the
# changelog. Loads only that fit. A few minutes; needs memory for one fit.
# Output: results/tables/diagnostics_03g.csv (key parameters) and a printed summary.
# ============================================================
library(rstan); library(tidyverse); library(here)

load(here("results", "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata"))   # object: fit

# sampler: divergences, tree depth, E-BFMI by chain
sp <- get_sampler_params(fit, inc_warmup = FALSE)
max_td <- 10   # cmdstanr default max_treedepth
samp <- tibble(chain = seq_along(sp),
               draws = map_int(sp, nrow),
               divergent = map_int(sp, \(x) sum(x[, "divergent__"])),
               at_max_treedepth = map_int(sp, \(x) sum(x[, "treedepth__"] >= max_td)),
               ebfmi = round(get_bfmi(fit), 2))
print(samp)

# all saved parameters: worst Rhat, lowest ESS
s_all <- summary(fit)$summary |> as.data.frame() |> rownames_to_column("parameter") |>
  filter(!str_starts(parameter, "lp__"), is.finite(Rhat))
cat("\nAll saved quantities:", nrow(s_all),
    "| Rhat > 1.01:", sum(s_all$Rhat > 1.01), "| ESS < 400:", sum(s_all$n_eff < 400), "\n")
cat("Worst Rhat:\n"); print(s_all |> slice_max(Rhat, n = 5) |> select(parameter, Rhat, n_eff))
cat("Lowest ESS:\n"); print(s_all |> slice_min(n_eff, n = 5) |> select(parameter, Rhat, n_eff))

# key parameters used in the forecasts
pars <- c("S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz", "S_bTrib", "T_bTrib",
          "S_bCovT", "TT_bCovT", "RE_sd", "RE_sdT", "TTRE_sd", "TTRE_sdT", "Pro_sd", "Pro_sdT")
key <- summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |> as.data.frame() |>
  rownames_to_column("parameter") |>
  transmute(parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`, rhat = Rhat, ess = n_eff) |>
  mutate(across(where(is.numeric), \(x) signif(x, 3)))
write.csv(key, here("results", "tables", "diagnostics_03g.csv"), row.names = FALSE)
cat("\nKey parameters: worst Rhat", max(key$rhat), "| lowest ESS", min(key$ess), "\n")
