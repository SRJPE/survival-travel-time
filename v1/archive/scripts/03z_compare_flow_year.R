# ============================================================
# 03z_compare_flow_year.R
#
# Compares the three telemetry models for the team's flow question:
#   M0 = flow, no year effect   (production fit 03g)
#   M1 = no flow, year effect   (03y_M1.R)
#   M2 = flow + year effect     (03y_M2.R)
# Questions: (1) does the flow slope survive once year-to-year survival
# differences are allowed (M2 vs M0)? (2) does flow add anything beyond
# a year effect (M2 vs M1)? (3) predictive fit (LOOIC, survival only).
# Run after 03g, M1 and M2 have all saved. Loads the large 03g fit, so
# don't run it while another fit is converting.
# Writes results/tables/flowyear_parameters.csv, flowyear_loo.csv,
# results/figures/flowyear_year_effects.png.
# ============================================================

library(rstan)
library(loo)
library(here)
library(tidyverse)

files <- c(M0 = here("results", "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata"),
           M1 = here("results", "fit_flowyear_M1.Rdata"),
           M2 = here("results", "fit_flowyear_M2.Rdata"))
labels <- c(M0 = "M0 flow (production)", M1 = "M1 no flow + year", M2 = "M2 flow + year")
pars <- c("S_bCov", "TT_bCov", "S_bCovT", "TT_bCovT", "S_bSz", "S_bReach", "RE_sd", "sd_yr")

one_model <- function(file, tag) {
  e <- new.env(); load(file, envir = e); f <- e$fit
  have <- intersect(pars, names(f@par_dims))
  par_tbl <- summary(f, pars = have, probs = c(0.025, 0.5, 0.975))$summary |>
    as.data.frame() |> rownames_to_column("parameter") |>
    transmute(model = tag, parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`,
              rhat = Rhat, ess = n_eff)
  ll <- loo::extract_log_lik(f, merge_chains = FALSE)
  lo <- loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
  yr <- if ("S_yr" %in% names(f@par_dims)) {
    summary(f, pars = "S_yr", probs = c(0.025, 0.5, 0.975))$summary |>
      as.data.frame() |>
      transmute(model = tag, year = e$Year_all, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`)
  } else NULL
  rm(e, f, ll); gc()
  list(par = par_tbl, loo = lo, yr = yr)
}
res <- imap(files, one_model)

par_tbl <- map(res, "par") |> list_rbind() |> mutate(across(where(is.numeric), \(x) round(x, 3)))
write.csv(par_tbl, here("results", "tables", "flowyear_parameters.csv"), row.names = FALSE)
par_tbl |>
  select(model, parameter, median, q2.5, q97.5, rhat) |>
  pivot_wider(names_from = model, values_from = c(median, q2.5, q97.5, rhat)) |>
  as.data.frame() |> print()

cmp <- loo::loo_compare(map(res, "loo"))
loo_tbl <- as.data.frame(cmp) |> rownames_to_column("model") |>
  transmute(model, elpd_diff, se_diff, looic, se_looic, p_loo) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))
write.csv(loo_tbl, here("results", "tables", "flowyear_loo.csv"), row.names = FALSE)
print(loo_tbl)

# year effects (logit scale, mainstem reaches 1-3), M1 vs M2
yr_tbl <- map(res, "yr") |> list_rbind()
p <- ggplot(yr_tbl, aes(factor(year), median, ymin = q2.5, ymax = q97.5, colour = model, shape = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c(M1 = "#E69F00", M2 = "#0072B2"), labels = labels[c("M1", "M2")]) +
  scale_shape_manual(values = c(M1 = 17, M2 = 16), labels = labels[c("M1", "M2")]) +
  labs(x = NULL, y = "Year effect on survival (logit)", colour = NULL, shape = NULL,
       title = "Year effects with and without flow",
       subtitle = "If flow explains the wet/dry contrast, M2's year effects shrink toward 0 relative to M1's") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here("results", "figures", "flowyear_year_effects.png"), p, width = 9, height = 5, dpi = 200, bg = "white")
p
