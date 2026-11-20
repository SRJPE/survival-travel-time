# ============================================================
# 06w_Run_CWT_Model_weekly_pcap.R
#
# CWT model with WEEKLY Knights Landing capture probability from
# BT-SPAS-X (scripts/CWT_weekly.stan), and, on exactly the same
# groups and recoveries, the current model with one pooled capture
# probability (scripts/CWT_nc.stan, as in 06). The difference
# between the two fits is the effect of weekly pCap alone.
#
# Release groups, covariates and travel times are built exactly as
# in 06 (same filters). Added:
#   1. Each group's weeks at KL from release to the end of the
#      BT-SPAS-X season (run week 31 = Julian week 22, early June),
#      as intervals of days since mid_release_date.
#   2. Weekly logit pCap prior (posterior mean and SD) for each
#      KL year x run week, from data/kl_weekly_pcap.csv (BT-SPAS-X
#      abundance fits, 25 KL years 1996-2025; the file used in the
#      simple simulations).
#   3. Groups dropped if their KL year has no BT-SPAS-X fit (2002,
#      2004-2006 have none) or they were released outside the season.
#      Recoveries after the season ends are left out of the counts
#      and travel times (both models, so the comparison is fair).
#
# Week convention: lubridate::week() (Julian weeks of 7 days from
# 1 Jan). Run weeks 1-9 = Julian weeks 45-53 of the previous
# calendar year, run weeks 10-31 = Julian weeks 1-22 (Korman et al.
# Table 1); KL year = the spring year.
#
# Outputs: results/fit_CWT_weekly.Rdata, results/fit_CWT_pooled_samegroups.Rdata,
# results/tables/cwt_weekly_*.csv, results/figures/cwt_weekly_*.png.
# Runs both models one after the other (4 chains each, same settings
# as 06). Fine alongside an overnight 03 run.
# ============================================================

rm(list = ls())

library(tidyverse)
library(here)
library(lubridate)
library(rstan)
library(cmdstanr)

source("scripts/02_GetData.R")

# Settings (as in 06) -----------------------------------------------------------------
use_size_effect <- 1
nchains <- 4
nwarmup <- 2000
nsamp <- 2000
adapt_delta <- 0.99
seed <- 123

release_distance_km <- c(Battle = 300.443064082, RBDD = 240.3631)
inv_logit <- function(x) exp(x) / (1 + exp(x))

kl_file <- here("data", "kl_weekly_pcap.csv")
if (!file.exists(kl_file)) stop("Copy kl_weekly_pcap.csv (from the simple sims inputs folder) into data/")
pcap_prior_file <- here("results", "pCap_one_site_skew_re_knights landing.rds")
if (!file.exists(pcap_prior_file)) stop("Pooled pCap prior not found: ", pcap_prior_file)

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

# Release groups: exactly as 06 ------------------------------------------------------
first_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

cwt_groups_06 <- drerelrec %>%
  filter(!is.na(relloc_area), relloc_area %in% names(release_distance_km),
         release_group_id != 873) %>%
  group_by(release_group_id) %>%
  summarise(
    release_location_name = first(release_location_name),
    relloc_area = first(relloc_area),
    run = first(run),
    year = first(year),
    month = first(month),
    mid_release_date = first(mid_release_date),
    Nrel = first(group_total_marked_N),
    avg_length = first(avg_length),
    Maxflowsac = first_finite(monthly_max_flow),
    cwt_recaptures = n_distinct(
      if_else(is.na(tag_code), NA_character_, paste(tag_code, date, forklength, sep = "_")),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ReachKM_cwt = unname(release_distance_km[relloc_area]),
         Rmult_cwt = ReachKM_cwt / 100) %>%
  filter(is.finite(Nrel), Nrel > 0, is.finite(avg_length),
         is.finite(Maxflowsac), is.finite(ReachKM_cwt)) %>%
  arrange(year, month, release_group_id)
message("CWT release groups in 06: ", nrow(cwt_groups_06))

# KL season weeks ----------------------------------------------------------------------
kl_weekly <- read.csv(kl_file) |>
  transmute(kl_year = as.integer(str_sub(fit, -4)), run_week = as.integer(week),
            lp_mean, lp_sd, p_median)

# date -> KL year and run week (NA outside the BT-SPAS-X season)
kl_week_of <- function(d) {
  jw <- as.integer(lubridate::week(d))
  yr <- as.integer(lubridate::year(d))
  tibble(kl_year = if_else(jw >= 45L, yr + 1L, yr),
         run_week = case_when(jw >= 45L ~ jw - 44L, jw <= 22L ~ jw + 9L, .default = NA_integer_))
}
season_end <- function(kl_year) make_date(kl_year, 1, 1) + days(7 * 22 - 1)   # last day of Julian week 22

grp_season <- cwt_groups_06 |>
  mutate(rel_date = as_date(mid_release_date), kl_week_of(rel_date),
         has_btspas = kl_year %in% kl_weekly$kl_year,
         keep = !is.na(run_week) & has_btspas)

# One row per distinct recovery (same definition as 06's count)
recs <- drerelrec |>
  filter(release_group_id %in% grp_season$release_group_id, !is.na(tag_code)) |>
  distinct(release_group_id, tag_code, date, forklength, .keep_all = TRUE) |>
  select(release_group_id, tag_code, date, forklength, cwt_tt) |>
  left_join(grp_season |> select(release_group_id, rel_date, rel_kl_year = kl_year),
            by = "release_group_id") |>
  mutate(date = as_date(date),
         in_season = is.na(date) | date <= season_end(rel_kl_year),
         kl_week_of(date)) |>
  left_join(kl_weekly |> select(kl_year, run_week, p_median), by = c("kl_year", "run_week"))

stopifnot(nrow(recs) == sum(cwt_groups_06$cwt_recaptures))   # same recoveries as 06's counts

# Check: is cwt_tt the recovery date minus mid_release_date?
recs |>
  mutate(tt_from_dates = as.numeric(date - rel_date), gap = tt_from_dates - cwt_tt) |>
  summarise(n = n(), n_gap_na = sum(is.na(gap)), gap_min = min(gap, na.rm = TRUE),
            gap_median = median(gap, na.rm = TRUE), gap_max = max(gap, na.rm = TRUE)) |>
  print()

# What the season rules remove, by year
filter_tbl <- grp_season |>
  left_join(recs |> group_by(release_group_id) |>
              summarise(rec_out_season = sum(!in_season),
                        rec_zero_pcap_week = sum(in_season & !is.na(p_median) & p_median < 1e-4)),
            by = "release_group_id") |>
  mutate(across(c(rec_out_season, rec_zero_pcap_week), \(x) replace_na(x, 0))) |>
  group_by(year) |>
  summarise(groups = n(), groups_kept = sum(keep),
            no_btspas_year = sum(!has_btspas), released_out_of_season = sum(is.na(run_week)),
            recoveries = sum(cwt_recaptures), rec_after_season = sum(rec_out_season[keep]),
            rec_in_unfished_week = sum(rec_zero_pcap_week[keep]), .groups = "drop")
print(as.data.frame(filter_tbl))
write.csv(filter_tbl, here("results", "tables", "cwt_weekly_filters_by_year.csv"), row.names = FALSE)

# Groups and counts used by BOTH models ------------------------------------------------
cwt_groups <- grp_season |>
  filter(keep) |>
  select(-cwt_recaptures) |>
  left_join(recs |> filter(in_season) |> count(release_group_id, name = "cwt_recaptures"),
            by = "release_group_id") |>
  mutate(cwt_recaptures = replace_na(cwt_recaptures, 0L)) |>
  arrange(year, month, release_group_id) |>
  mutate(grp = row_number())

size_mu <- mean(cwt_groups$avg_length); size_sd <- sd(cwt_groups$avg_length)
flow_mu <- mean(cwt_groups$Maxflowsac); flow_sd <- sd(cwt_groups$Maxflowsac)
cwt_groups <- cwt_groups |>
  mutate(Sz_cwt = (avg_length - size_mu) / size_sd,
         CovX_cwt = (Maxflowsac - flow_mu) / flow_sd)

# Each group's weeks at KL: release date to season end, as [lo, hi) in days since release.
# cwt_tt = (recovery date - mid_release_date) - 0.5 for every recovery (checked above),
# so a fish recovered on day `off` after release has cwt_tt = off - 0.5. Week edges are
# shifted by -0.5 day to match, otherwise recoveries on the first day of a week would be
# assigned to the week before.
make_segments <- function(grp, rel_date, kl_year) {
  tibble(day = seq(rel_date, season_end(kl_year), by = "day")) |>
    mutate(off = as.numeric(day - rel_date), kl_week_of(day)) |>
    group_by(kl_year, run_week) |>
    summarise(lo = pmax(min(off) - 0.5, 0), hi = max(off) + 0.5, .groups = "drop") |>
    arrange(lo) |>
    mutate(grp = grp)
}
segs <- cwt_groups |>
  select(grp, rel_date, kl_year) |>
  pmap(make_segments) |>
  list_rbind()

# weekly pCap parameters: every KL year x run week some group passes through
wk_tbl <- segs |>
  distinct(kl_year, run_week) |>
  arrange(kl_year, run_week) |>
  left_join(kl_weekly, by = c("kl_year", "run_week")) |>
  mutate(wk = row_number())
stopifnot(!anyNA(wk_tbl$lp_mean))
segs <- segs |>
  left_join(wk_tbl |> select(kl_year, run_week, wk), by = c("kl_year", "run_week")) |>
  arrange(grp, lo) |>
  mutate(seg = row_number())
seg_idx <- segs |> group_by(grp) |> summarise(seg_start = min(seg), seg_end = max(seg))
stopifnot(identical(seg_idx$grp, cwt_groups$grp))

# Travel times (as 06: finite, > 0), in season, with the week each fish was recovered
cwt_tt_obs <- recs |>
  filter(in_season, is.finite(cwt_tt), cwt_tt > 0) |>
  distinct(release_group_id, tag_code, date, forklength, cwt_tt) |>
  inner_join(cwt_groups |> select(release_group_id, grp), by = "release_group_id") |>
  inner_join(segs |> select(grp, lo, hi, wk), by = join_by(grp, cwt_tt >= lo, cwt_tt < hi)) |>
  arrange(grp, date, tag_code)
message("Travel-time observations: ", nrow(cwt_tt_obs))

cwt_data <- list(
  Ncwtgrp = nrow(cwt_groups), Nyrs = n_distinct(cwt_groups$year), Nobs = nrow(cwt_tt_obs),
  UseSizeEffect = use_size_effect,
  Rmult_cwt = cwt_groups$Rmult_cwt, Sz_cwt = cwt_groups$Sz_cwt,
  CovX_cwt = cwt_groups$CovX_cwt, ReachKM_cwt = cwt_groups$ReachKM_cwt,
  ObsTT = cwt_tt_obs$cwt_tt, TTind = cwt_tt_obs$grp,
  Nrel = as.integer(round(cwt_groups$Nrel)), cwt_recaptures = as.integer(cwt_groups$cwt_recaptures),
  # weekly pCap (CWT_weekly.stan)
  Nwk = nrow(wk_tbl), lp_mu = wk_tbl$lp_mean, lp_sd = wk_tbl$lp_sd,
  Nseg = nrow(segs), seg_lo = segs$lo, seg_hi = segs$hi, seg_wk = segs$wk,
  seg_start = seg_idx$seg_start, seg_end = seg_idx$seg_end, obs_wk = cwt_tt_obs$wk
)

# pooled pCap prior (CWT_nc.stan), as in 06
pcap_draws <- as.data.frame(readRDS(pcap_prior_file), pars = "logit_pCap") |> unlist()
pcap_draws <- pcap_draws[is.finite(pcap_draws)]
pooled_data <- c(cwt_data[c("Ncwtgrp", "Nyrs", "Nobs", "UseSizeEffect", "Rmult_cwt", "Sz_cwt",
                            "CovX_cwt", "ReachKM_cwt", "ObsTT", "TTind", "Nrel", "cwt_recaptures")],
                 list(lt_mu_pCap = mean(pcap_draws), lt_sd_pCap = sd(pcap_draws)))

save(cwt_data, pooled_data, cwt_groups, cwt_tt_obs, segs, wk_tbl,
     file = here("results", "cwt_weekly_stan_data.Rdata"))

# Fit both ------------------------------------------------------------------------------
fix_csv <- function(f) {
  x <- readLines(f)
  com <- startsWith(x, "#")
  x[com] <- gsub("= true", "= 1", gsub("= false", "= 0", x[com]))
  out <- sub("\\.csv$", "_rstan.csv", f)
  writeLines(x, out)
  out
}
fit_one <- function(stan_file, data, tag, init = NULL) {
  cs_dir <- here("results", paste0("cmdstan_cwt_", tag))
  dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
  mod <- cmdstan_model(here("scripts", stan_file))
  fit_cs <- mod$sample(data = data, chains = nchains, parallel_chains = nchains,
                       iter_warmup = nwarmup, iter_sampling = nsamp, seed = seed,
                       adapt_delta = adapt_delta, refresh = 500, output_dir = cs_dir,
                       init = init)
  print(fit_cs$diagnostic_summary())
  rstan::read_stan_csv(map_chr(fit_cs$output_files(), fix_csv))
}

# Pooled fit: reuse the saved one if it was fitted to exactly these data
pooled_file <- here("results", "fit_CWT_pooled_samegroups.Rdata")
if (file.exists(pooled_file)) {
  e <- new.env(); load(pooled_file, envir = e)
  stopifnot(identical(e$pooled_data, pooled_data))
  fit_cwt_pooled <- e$fit_cwt_pooled; rm(e)
  message("Pooled fit loaded from ", pooled_file)
} else {
  fit_cwt_pooled <- fit_one("CWT_nc.stan", pooled_data, "pooled_samegroups")
  save(fit_cwt_pooled, pooled_data, cwt_groups, cwt_tt_obs, file = pooled_file)
}

# Weekly fit starts near the pooled posterior (random starts can put travel
# times far in the tails, where the weekly-passage probabilities underflow).
# Small per-chain jitter so the chains still start apart; weekly pCap at the
# BT-SPAS-X means.
pmed <- function(p) apply(as.matrix(rstan::extract(fit_cwt_pooled, p)[[1]]), 2, median)
set.seed(seed)
init_weekly <- map(seq_len(nchains), \(ch) {
  jit <- \(x) x + rnorm(length(x), 0, 0.1)
  pos <- \(x) x * exp(rnorm(length(x), 0, 0.1))
  list(S_bReach = jit(pmed("S_bReach")), T_bReach = jit(pmed("T_bReach")),
       S_bCov = jit(pmed("S_bCov")), TT_bCov = jit(pmed("TT_bCov")),
       S_bSz = jit(pmed("S_bSz")), T_bSz = jit(pmed("T_bSz")),
       z_S = jit(pmed("z_S")), z_TT = jit(pmed("z_TT")),
       sd_cwtSre = pos(pmed("sd_cwtSre")), sd_cwtTTre = pos(pmed("sd_cwtTTre")),
       Pro_sd = pos(pmed("Pro_sd")), z_pcap = rep(0, cwt_data$Nwk))
})

fit_cwt_weekly <- fit_one("CWT_weekly.stan", cwt_data, "weekly", init = init_weekly)
save(fit_cwt_weekly, cwt_data, cwt_groups, cwt_tt_obs, segs, wk_tbl,
     file = here("results", "fit_CWT_weekly.Rdata"))

# Compare --------------------------------------------------------------------------------
pars <- c("S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz",
          "sd_cwtSre", "sd_cwtTTre", "Pro_sd")
fits <- list(pooled = fit_cwt_pooled, weekly = fit_cwt_weekly)

par_tbl <- imap(fits, \(f, nm)
  summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |>
    as.data.frame() |> rownames_to_column("parameter") |>
    transmute(model = nm, parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`,
              rhat = Rhat, ess = n_eff)) |>
  list_rbind() |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write.csv(par_tbl, here("results", "tables", "cwt_weekly_vs_pooled_parameters.csv"), row.names = FALSE)
par_tbl |>
  select(model, parameter, median, q2.5, q97.5, rhat, ess) |>
  pivot_wider(names_from = model, values_from = c(median, q2.5, q97.5, rhat, ess)) |>
  as.data.frame() |> print()

# P(recovered at KL | survived): pooled = one pCap; weekly = per group
p_rec <- bind_rows(
  tibble(model = "pooled", grp = cwt_groups$grp,
         median = median(inv_logit(rstan::extract(fit_cwt_pooled, "logit_pCap_Sim")[[1]]))),
  tibble(model = "weekly", grp = cwt_groups$grp,
         median = apply(rstan::extract(fit_cwt_weekly, "p_rec_given_surv")[[1]], 2, median))
)

# group survival to KL (posterior median), and its mean by year
surv_grp <- imap(fits, \(f, nm)
  tibble(model = nm, grp = cwt_groups$grp,
         surv = apply(rstan::extract(f, "surv_cwt")[[1]], 2, median))) |>
  list_rbind() |>
  left_join(p_rec |> rename(p_rec = median), by = c("model", "grp")) |>
  left_join(cwt_groups |> select(grp, release_group_id, year, relloc_area, cwt_recaptures, Nrel),
            by = "grp")
write.csv(surv_grp, here("results", "tables", "cwt_weekly_vs_pooled_groups.csv"), row.names = FALSE)

surv_year <- surv_grp |>
  group_by(model, year) |>
  summarise(groups = n(), zero_rec = sum(cwt_recaptures == 0), mean_surv = mean(surv),
            median_p_rec = median(p_rec), .groups = "drop")
write.csv(surv_year, here("results", "tables", "cwt_weekly_vs_pooled_by_year.csv"), row.names = FALSE)
surv_year |>
  pivot_wider(names_from = model, values_from = c(mean_surv, median_p_rec)) |>
  mutate(across(where(is.double), \(x) signif(x, 3))) |>
  as.data.frame() |> print()

# Figure: mean group survival to KL by year; Okabe-Ito colours plus shapes
p <- ggplot(surv_year, aes(year, mean_surv, colour = model, shape = model)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(pooled = "#D55E00", weekly = "#0072B2"),
                      labels = c(pooled = "One pooled pCap (06)", weekly = "Weekly BT-SPAS-X pCap")) +
  scale_shape_manual(values = c(pooled = 17, weekly = 16),
                     labels = c(pooled = "One pooled pCap (06)", weekly = "Weekly BT-SPAS-X pCap")) +
  labs(x = NULL, y = "Mean group survival to Knights Landing", colour = NULL, shape = NULL,
       title = "CWT survival by year: pooled vs weekly capture probability",
       subtitle = "Same groups and recoveries in both fits") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "cwt_weekly_vs_pooled_by_year.png"), p,
       width = 8, height = 5, dpi = 200, bg = "white")
p
