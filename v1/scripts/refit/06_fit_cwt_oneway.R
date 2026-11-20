# ---- srjpe_survival: refit step 06, one-way Lorenzen CWT model (median anchor)
# ---- made from scripts/06o_Run_CWT_oneway.R by setup/build_from_old_project.R (exact edits listed there)

# ============================================================
# 06o_Run_CWT_oneway.R
#
# ONE-WAY, telemetry-anchored CWT model (joint model plan, 8 Oct 2026).
# Telemetry is fitted alone (03g, the handoff stm) and is not changed.
# CWT survival comes from telemetry mortality per day through Lorenzen
# (scripts/CWT_oneway.stan); the CWT data estimate only the Lorenzen
# exponent c, a CWT offset, CWT group effects, CWT travel time and weekly
# KL capture probability (BT-SPAS-X priors).
#
# Telemetry anchor: H_tel = total hazard over mainstem telemetry fish
# (sum of -log SurvRelSac), with T_i = TT_RelSac and fork lengths from 02.
# Three runs: H_tel at its posterior median (main), 2.5% and 97.5%
# quantiles (sensitivity). Travel times T_i at their posterior medians.
#
# CWT data exactly as 06w (same groups, recoveries, weeks, priors).
# Each run: 4 chains x 1000/1000, adapt_delta 0.95, about 15 to 30 min.
# Needs: results/fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata (03g),
#        results/fit_CWT_weekly.Rdata (06w, for starting values and comparison),
#        data/kl_weekly_pcap.csv.
# Outputs: results/fit_CWT_oneway_<run>.Rdata, results/tables/cwt_oneway_*.csv,
#          results/figures/cwt_oneway_survival_by_year.png
# ============================================================

rm(list = ls())

library(tidyverse)
library(here)
library(lubridate)
library(rstan)
library(cmdstanr)

source(here("scripts", "refit", "02_get_data.R"))

nchains <- 4
nwarmup <- 1000
nsamp <- 1000
adapt_delta <- 0.95
seed <- 123
L_r <- 70                       # reference length (mm)
c_prior <- c(mu = 1, sd = 0.25) # Lorenzen exponent, about 0.5 to 1.5
use_offset <- 1L

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ---- Telemetry anchor (fixed) -------------------------------------------------------
e <- new.env(); load(here("results", "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata"), envir = e)
S_tel <- rstan::extract(e$fit, "SurvRelSac")[[1]]
T_tel_d <- rstan::extract(e$fit, "TT_RelSac")[[1]]
rm(e); gc()
stopifnot(ncol(S_tel) == length(FL))
ok <- is.finite(FL)
H_draws <- rowSums(-log(S_tel[, ok, drop = FALSE]))
T_tel <- apply(T_tel_d[, ok, drop = FALSE], 2, median)
L_tel <- FL[ok]
rm(S_tel, T_tel_d); gc()
H_runs <- c(Hmed = median(H_draws), Hlo = unname(quantile(H_draws, 0.025)),
            Hhi = unname(quantile(H_draws, 0.975)))
H_runs <- H_runs["Hmed"]   # handoff: median anchor only (2.5% / 97.5% anchors were sensitivity runs, 11 Oct)
print(round(H_runs, 1))
message("Telemetry fish: ", length(L_tel), "; length ", paste(range(L_tel), collapse = " to "), " mm")

# ---- CWT data: exactly as 06w ----------------------------------------------------------
release_distance_km <- c(Battle = 300.443064082, RBDD = 240.3631)
kl_file <- here("data", "kl_weekly_pcap.csv")
if (!file.exists(kl_file)) stop("Copy kl_weekly_pcap.csv into data/")

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
  mutate(ReachKM_cwt = unname(release_distance_km[relloc_area])) %>%
  filter(is.finite(Nrel), Nrel > 0, is.finite(avg_length),
         is.finite(Maxflowsac), is.finite(ReachKM_cwt)) %>%
  arrange(year, month, release_group_id)

kl_weekly <- read.csv(kl_file) |>
  transmute(kl_year = as.integer(str_sub(fit, -4)), run_week = as.integer(week), lp_mean, lp_sd)
kl_week_of <- function(d) {
  jw <- as.integer(lubridate::week(d))
  yr <- as.integer(lubridate::year(d))
  tibble(kl_year = if_else(jw >= 45L, yr + 1L, yr),
         run_week = case_when(jw >= 45L ~ jw - 44L, jw <= 22L ~ jw + 9L, .default = NA_integer_))
}
season_end <- function(kl_year) make_date(kl_year, 1, 1) + days(7 * 22 - 1)

grp_season <- cwt_groups_06 |>
  mutate(rel_date = as_date(mid_release_date), kl_week_of(rel_date),
         keep = !is.na(run_week) & kl_year %in% kl_weekly$kl_year)
recs <- drerelrec |>
  filter(release_group_id %in% grp_season$release_group_id, !is.na(tag_code)) |>
  distinct(release_group_id, tag_code, date, forklength, .keep_all = TRUE) |>
  select(release_group_id, tag_code, date, forklength, cwt_tt) |>
  left_join(grp_season |> select(release_group_id, rel_kl_year = kl_year), by = "release_group_id") |>
  mutate(date = as_date(date), in_season = is.na(date) | date <= season_end(rel_kl_year))
stopifnot(nrow(recs) == sum(cwt_groups_06$cwt_recaptures))

cwt_groups <- grp_season |>
  filter(keep) |>
  select(-cwt_recaptures) |>
  left_join(recs |> filter(in_season) |> count(release_group_id, name = "cwt_recaptures"),
            by = "release_group_id") |>
  mutate(cwt_recaptures = replace_na(cwt_recaptures, 0L)) |>
  arrange(year, month, release_group_id) |>
  mutate(grp = row_number(),
         Sz_cwt = (avg_length - mean(avg_length)) / sd(avg_length),        # as 06w (travel time only)
         CovX_cwt = (Maxflowsac - mean(Maxflowsac)) / sd(Maxflowsac))      # as 06w (travel time only)

make_segments <- function(grp, rel_date, kl_year) {
  tibble(day = seq(rel_date, season_end(kl_year), by = "day")) |>
    mutate(off = as.numeric(day - rel_date), kl_week_of(day)) |>
    group_by(kl_year, run_week) |>
    summarise(lo = pmax(min(off) - 0.5, 0), hi = max(off) + 0.5, .groups = "drop") |>
    arrange(lo) |>
    mutate(grp = grp)
}
segs <- cwt_groups |> select(grp, rel_date, kl_year) |> pmap(make_segments) |> list_rbind()
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

cwt_tt_obs <- recs |>
  filter(in_season, is.finite(cwt_tt), cwt_tt > 0) |>
  distinct(release_group_id, tag_code, date, forklength, cwt_tt) |>
  inner_join(cwt_groups |> select(release_group_id, grp), by = "release_group_id") |>
  inner_join(segs |> select(grp, lo, hi, wk), by = join_by(grp, cwt_tt >= lo, cwt_tt < hi)) |>
  arrange(grp, date, tag_code)
message("CWT groups: ", nrow(cwt_groups), "; recoveries: ", sum(cwt_groups$cwt_recaptures),
        "; travel times: ", nrow(cwt_tt_obs))   # 06w: 357, 3013, 3013

base_data <- list(
  Ncwtgrp = nrow(cwt_groups), Nyrs = n_distinct(cwt_groups$year), Nobs = nrow(cwt_tt_obs),
  UseSizeEffect = 1L,
  Sz_cwt = cwt_groups$Sz_cwt, CovX_cwt = cwt_groups$CovX_cwt, ReachKM_cwt = cwt_groups$ReachKM_cwt,
  ObsTT = cwt_tt_obs$cwt_tt, TTind = cwt_tt_obs$grp,
  Nrel = as.integer(round(cwt_groups$Nrel)), cwt_recaptures = as.integer(cwt_groups$cwt_recaptures),
  Nwk = nrow(wk_tbl), lp_mu = wk_tbl$lp_mean, lp_sd = wk_tbl$lp_sd,
  Nseg = nrow(segs), seg_lo = segs$lo, seg_hi = segs$hi, seg_wk = segs$wk,
  seg_start = seg_idx$seg_start, seg_end = seg_idx$seg_end, obs_wk = cwt_tt_obs$wk,
  Ntel = length(T_tel), T_tel = T_tel, log_Lratio_tel = log(L_r / L_tel),
  log_Lratio_cwt = log(L_r / cwt_groups$avg_length),
  c_mu = c_prior[["mu"]], c_sd = c_prior[["sd"]], UseOffset = use_offset)

# ---- Starting values: travel-time terms from the 06w weekly fit ----------------------
e <- new.env(); load(here("results", "fit_CWT_weekly.Rdata"), envir = e)
pmed <- \(p) apply(as.matrix(rstan::extract(e$fit_cwt_weekly, p)[[1]]), 2, median)
tt_start <- list(T_bReach = pmed("T_bReach"), TT_bCov = pmed("TT_bCov"), T_bSz = pmed("T_bSz"),
                 z_TT = pmed("z_TT"), sd_cwtTTre = pmed("sd_cwtTTre"), Pro_sd = pmed("Pro_sd"))
surv_06w <- apply(rstan::extract(e$fit_cwt_weekly, "surv_cwt")[[1]], 2, median)
rm(e); gc()
stopifnot(length(tt_start$z_TT) == base_data$Ncwtgrp, length(surv_06w) == base_data$Ncwtgrp)

set.seed(seed)
make_inits <- function() map(seq_len(nchains), \(ch) {
  jit <- \(x) x + rnorm(length(x), 0, 0.05)
  pos <- \(x) x * exp(rnorm(length(x), 0, 0.05))
  list(T_bReach = jit(tt_start$T_bReach), TT_bCov = jit(tt_start$TT_bCov), T_bSz = jit(tt_start$T_bSz),
       z_TT = jit(tt_start$z_TT), sd_cwtTTre = pos(tt_start$sd_cwtTTre), Pro_sd = pos(tt_start$Pro_sd),
       z_pcap = rep(0, base_data$Nwk), c = 1 + rnorm(1, 0, 0.05), delta = rnorm(1, 0, 0.05),
       z_h = rep(0, base_data$Ncwtgrp), sd_h = 0.2)
})

# ---- Fit the three runs --------------------------------------------------------------
fix_csv <- function(f) {
  x <- readLines(f)
  com <- startsWith(x, "#")
  x[com] <- gsub("= true", "= 1", gsub("= false", "= 0", x[com]))
  out <- sub("\\.csv$", "_rstan.csv", f)
  writeLines(x, out)
  out
}
mod <- cmdstan_model(here("stan", "CWT_oneway.stan"))
fit_run <- function(H, run) {
  cs_dir <- here("results", paste0("cmdstan_cwt_oneway_", run))
  dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
  fit_cs <- mod$sample(data = c(base_data, list(H_tel = H)), init = make_inits(),
                       chains = nchains, parallel_chains = nchains,
                       iter_warmup = nwarmup, iter_sampling = nsamp, seed = seed,
                       adapt_delta = adapt_delta, refresh = 250, output_dir = cs_dir)
  message("Run ", run, ":"); print(fit_cs$diagnostic_summary())
  fit <- rstan::read_stan_csv(map_chr(fit_cs$output_files(), fix_csv))
  save(fit, base_data, cwt_groups, segs, wk_tbl, H, file = here("results", paste0("fit_CWT_oneway_", run, ".Rdata")))
  fit
}
fits <- imap(H_runs, fit_run)

# ---- Summaries --------------------------------------------------------------------------
pars <- c("c", "delta", "sd_h", "M_r", "T_bReach", "T_bSz", "TT_bCov", "sd_cwtTTre", "Pro_sd")
par_tbl <- imap(fits, \(f, run)
  summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |>
    as.data.frame() |> rownames_to_column("parameter") |>
    transmute(run, parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`, rhat = Rhat, ess = n_eff)) |>
  list_rbind() |>
  mutate(across(where(is.numeric), \(x) signif(x, 3)))
write.csv(par_tbl, here("results", "tables", "cwt_oneway_parameters.csv"), row.names = FALSE)
par_tbl |> select(run, parameter, median, q2.5, q97.5, rhat) |>
  pivot_wider(names_from = run, values_from = c(median, q2.5, q97.5, rhat)) |>
  as.data.frame() |> print()

# BT-SPAS-X capture (prior only, arrival-weighted with the main run's travel times)
wk_tbl <- wk_tbl |> mutate(p_mean = map2_dbl(lp_mean, lp_sd, \(m, s) mean(plogis(rnorm(20000, m, s)))))
mu_main <- log(apply(rstan::extract(fits$Hmed, "pTT")[[1]], 2, median))
sd_main <- median(rstan::extract(fits$Hmed, "Pro_sd")[[1]])
p_btspas <- segs |>
  left_join(wk_tbl |> select(wk, p_mean), by = "wk") |>
  mutate(f = plnorm(hi, mu_main[grp], sd_main) - plnorm(lo, mu_main[grp], sd_main)) |>
  group_by(grp) |> summarise(p_btspas = sum(f * p_mean)) |> arrange(grp)

grp_tbl <- imap(fits, \(f, run) tibble(
  run, grp = cwt_groups$grp,
  surv = apply(rstan::extract(f, "surv_cwt")[[1]], 2, median),
  p_rec = apply(rstan::extract(f, "p_rec_given_surv")[[1]], 2, median),
  expected = colMeans(sweep(rstan::extract(f, "pKL")[[1]], 2, base_data$Nrel, `*`)))) |>
  list_rbind() |>
  left_join(cwt_groups |> select(grp, year, relloc_area, avg_length, Nrel, cwt_recaptures), by = "grp") |>
  left_join(p_btspas, by = "grp")
write.csv(grp_tbl, here("results", "tables", "cwt_oneway_groups.csv"), row.names = FALSE)

year_tbl <- grp_tbl |>
  group_by(run, year) |>
  summarise(groups = n(), observed = sum(cwt_recaptures), expected = sum(expected),
            mean_surv = mean(surv),
            capture_fitted = weighted.mean(p_rec, Nrel), capture_btspas = weighted.mean(p_btspas, Nrel),
            .groups = "drop") |>
  mutate(capture_shift = capture_fitted / capture_btspas)   # > 1: model raised capture above BT-SPAS-X
write.csv(year_tbl, here("results", "tables", "cwt_oneway_by_year.csv"), row.names = FALSE)
year_tbl |> filter(run == "Hmed") |> select(-run) |>
  mutate(across(where(is.double), \(x) signif(x, 3))) |> as.data.frame() |> print()

# ---- Figure: survival to KL by year, one-way (three anchors) vs 06w CWT-only ------------
fig_tbl <- bind_rows(
  imap(fits, \(f, run) {
    s <- rstan::extract(f, "surv_cwt")[[1]]
    Y <- model.matrix(~ factor(year) - 1, cwt_groups); Y <- sweep(Y, 2, colSums(Y), "/")
    m <- s %*% Y
    tibble(source = run, year = as.integer(sort(unique(cwt_groups$year))),
           median = apply(m, 2, median), q2.5 = apply(m, 2, quantile, 0.025),
           q97.5 = apply(m, 2, quantile, 0.975))
  }) |> list_rbind(),
  tibble(source = "cwt_only_06w", year = cwt_groups$year, s = surv_06w) |>
    group_by(source, year) |> summarise(median = mean(s), .groups = "drop") |>
    mutate(q2.5 = NA_real_, q97.5 = NA_real_))
src_lab <- c(Hmed = "One-way, telemetry anchor (median)", Hlo = "One-way, anchor 2.5%",
             Hhi = "One-way, anchor 97.5%", cwt_only_06w = "CWT only, weekly capture (06w)")
src_col <- c(Hmed = "#0072B2", Hlo = "#56B4E9", Hhi = "#009E73", cwt_only_06w = "#E69F00")
src_shp <- c(Hmed = 16, Hlo = 1, Hhi = 2, cwt_only_06w = 15)
p <- ggplot(fig_tbl, aes(year, median, ymin = q2.5, ymax = q97.5, colour = source, shape = source)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.3, na.rm = TRUE) +
  scale_colour_manual(values = src_col, labels = src_lab) +
  scale_shape_manual(values = src_shp, labels = src_lab) +
  labs(x = NULL, y = "Mean group survival to Knights Landing", colour = NULL, shape = NULL,
       title = "CWT survival anchored to telemetry (one-way Lorenzen model)",
       subtitle = "Telemetry mortality per day, scaled by length and travel time; telemetry not updated by CWT") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "cwt_oneway_survival_by_year.png"), p,
       width = 10, height = 5.5, dpi = 200, bg = "white")
p
