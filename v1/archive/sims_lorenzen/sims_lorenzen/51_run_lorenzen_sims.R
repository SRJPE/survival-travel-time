# ============================================================
# 51_run_lorenzen_sims.R
#
# Simulation check of the ONE-WAY, telemetry-anchored CWT model
# (CWT_oneway.stan, the same file used on real data in 06o).
# Noble's question: estimate survival from telemetry, predict survival of
# smaller CWT fish with Lorenzen, and see whether CWT recoveries match.
# Questions:
#   1. Does the one-way model recover the Lorenzen exponent c and the CWT
#      offset delta when the truth is known?
#   2. How large a CWT vs telemetry difference can it detect, given real
#      BT-SPAS-X uncertainty in Knights Landing capture?
#   3. What happens if mortality is really per km, not per day?
#
# Free standing: needs only this folder (two Stan files, inputs/ with the
# BT-SPAS-X weekly files) plus R with tidyverse and cmdstanr/CmdStan.
# Run from this folder: setwd() to it, then source("51_run_lorenzen_sims.R").
#
# Truth (each replicate):
#   - 10 years, each matched to a real KL year; true weekly capture is one
#     BT-SPAS-X posterior draw for that year; the model gets that week's
#     BT-SPAS-X posterior mean and SD of logit capture as a prior (as 06o).
#   - Telemetry: 2,500 fish, fork length about 100 mm (80 to 140), detection
#     at KL 0.9 (known).
#   - CWT: 40 groups (4 per year), 20,000 to 80,000 fish, mean length 45 to
#     85 mm, released at the start of run weeks 8 to 20.
#   - Travel time to KL (300 km), per fish lognormal: median
#     10 x (100 / L)^1.5 days (small fish slower), log SD 0.45, plus a CWT
#     group effect (SD 0.2).
#   - Survival per fish: exp(-h x T), h = M_r (L_r / L)^c (per day),
#     M_r = 0.13 per day at L_r = 70 mm; CWT groups also get exp(delta + e_g),
#     e_g ~ N(0, 0.2).
# Scenarios:
#   S1_match     c = 1,   delta = 0     (Lorenzen per day holds)
#   S2_cwt_lower c = 1,   delta = 0.4   (CWT hazard 1.5x: hatchery or tag effect)
#   S3_c_half    c = 0.5, delta = 0     (weaker size effect than the prior centre)
#   S4_per_km    mortality per km, not per day: S = exp(-k (L_r/L)^c D/100),
#                k set so a 100 mm fish has the same survival as in S1
# Fits per replicate: telemetry (sim_tel_lorenzen.stan) -> H and travel
# times -> one-way CWT model (CWT_oneway.stan), 4 chains x 500/500.
# Results are appended to results/reps.csv after every replicate, so an
# interruption loses at most one replicate; rerunning skips finished ones.
# Expect about 2 to 5 min per replicate; 4 scenarios x 20 replicates.
# ============================================================

library(tidyverse)
library(cmdstanr)

n_reps <- 20
scenarios <- tribble(
  ~scenario,      ~c_true, ~delta_true, ~mortality,
  "S1_match",      1.0,     0.0,        "per_day",
  "S2_cwt_lower",  1.0,     0.4,        "per_day",
  "S3_c_half",     0.5,     0.0,        "per_day",
  "S4_per_km",     1.0,     0.0,        "per_km")
truth <- list(M_r = 0.13, L_r = 70, D = 300, sigma_T = 0.45, sd_tt_grp = 0.2, sd_h = 0.2,
              p_tel = 0.9, n_tel = 2500, n_years = 10, grp_per_year = 4,
              rel_weeks = 8:20, last_week = 31)
med_T <- function(L) 10 * (100 / L)^1.5            # median days to KL (300 km)

dir.create("results", showWarnings = FALSE)
out_file <- "results/reps.csv"
done <- if (file.exists(out_file)) distinct(read.csv(out_file), scenario, rep) else tibble(scenario = character(), rep = integer())

kl_weekly <- read.csv("inputs/kl_weekly_pcap.csv")
kl_draws <- readRDS("inputs/kl_pcap_draws.rds")
mod_tel <- cmdstan_model("sim_tel_lorenzen.stan")
mod_cwt <- cmdstan_model("CWT_oneway.stan")

# ---- simulate one replicate ------------------------------------------------------------
simulate_rep <- function(sc, r) {
  set.seed(1000 * match(sc$scenario, scenarios$scenario) + r)
  hazard <- function(L, extra = 0) truth$M_r * (truth$L_r / L)^sc$c_true * exp(extra)
  k_km <- truth$M_r * med_T(100) / (truth$D / 100)   # per-km rate giving the per-day survival at 100 mm
  surv_fish <- function(L, T, extra = 0) {
    if (sc$mortality == "per_day") exp(-hazard(L, extra) * T)
    else exp(-k_km * (truth$L_r / L)^sc$c_true * exp(extra) * truth$D / 100) * T^0
  }

  # telemetry
  L_tel <- pmin(pmax(rnorm(truth$n_tel, 100, 12), 80), 140)
  T_tel <- rlnorm(truth$n_tel, log(med_T(L_tel)), truth$sigma_T)
  alive <- rbinom(truth$n_tel, 1, surv_fish(L_tel, T_tel))
  det <- rbinom(truth$n_tel, 1, alive * truth$p_tel)
  z_tel <- (L_tel - mean(L_tel)) / sd(L_tel)

  # CWT
  yrs <- sample(names(kl_draws), truth$n_years)
  grp <- tibble(year = rep(seq_len(truth$n_years), each = truth$grp_per_year)) |>
    mutate(kl_year = yrs[year],
           draw = rep(sample(200, truth$n_years, replace = TRUE), each = truth$grp_per_year),
           rel_week = sample(truth$rel_weeks, n(), replace = TRUE),
           Nrel = round(runif(n(), 20000, 80000)),
           L = runif(n(), 45, 85),
           tt_re = rnorm(n(), 0, truth$sd_tt_grp),
           h_re = sc$delta_true + rnorm(n(), 0, truth$sd_h),
           grp = row_number())
  sim_grp <- function(g) {
    T <- rlnorm(g$Nrel, log(med_T(g$L)) + g$tt_re, truth$sigma_T)
    s <- surv_fish(g$L, T, g$h_re)
    wk <- g$rel_week + floor(T / 7)
    p_week <- kl_draws[[g$kl_year]][g$draw, ]
    cap <- rbinom(g$Nrel, 1, s * ifelse(wk <= truth$last_week, p_week[pmin(wk, 31)], 0))
    list(true_surv = mean(s),
         true_prec = mean(ifelse(wk <= truth$last_week, p_week[pmin(wk, 31)], 0) * s) / mean(s),
         rec = tibble(grp = g$grp, T = T[cap == 1], run_week = wk[cap == 1]))
  }
  sims <- map(seq_len(nrow(grp)), \(i) sim_grp(grp[i, ]))
  grp <- grp |> mutate(true_surv = map_dbl(sims, "true_surv"), true_prec = map_dbl(sims, "true_prec"),
                       recoveries = map_int(sims, \(x) nrow(x$rec)))
  recs <- map(sims, "rec") |> list_rbind()

  # weeks: release week to the end of the season, shared by groups in the same year
  segs <- grp |>
    select(grp, kl_year, rel_week) |>
    mutate(run_week = map(rel_week, \(w) w:truth$last_week)) |>
    unnest(run_week) |>
    mutate(lo = (run_week - rel_week) * 7, hi = lo + 7)
  wk_tbl <- segs |> distinct(kl_year, run_week) |>
    left_join(kl_weekly |> select(fit, week, lp_mean, lp_sd), by = c("kl_year" = "fit", "run_week" = "week")) |>
    mutate(wk = row_number())
  segs <- segs |> left_join(wk_tbl |> select(kl_year, run_week, wk), by = c("kl_year", "run_week")) |>
    arrange(grp, lo) |> mutate(seg = row_number())
  seg_idx <- segs |> group_by(grp) |> summarise(seg_start = min(seg), seg_end = max(seg))
  recs <- recs |> left_join(grp |> select(grp, kl_year), by = "grp") |>
    left_join(wk_tbl |> select(kl_year, run_week, wk), by = c("kl_year", "run_week")) |>
    arrange(grp, T)

  list(grp = grp, recs = recs, segs = segs, wk_tbl = wk_tbl, seg_idx = seg_idx,
       tel = list(z = z_tel, L = L_tel, det = det, T_obs = T_tel[det == 1], idx = which(det == 1)))
}

# ---- fit one replicate -----------------------------------------------------------------
fit_rep <- function(sc, r) {
  sim <- simulate_rep(sc, r)
  t0 <- Sys.time()
  ft <- mod_tel$sample(data = list(N = length(sim$tel$z), z = sim$tel$z, det = sim$tel$det,
                                   p = truth$p_tel, Nobs = length(sim$tel$idx),
                                   T_obs = sim$tel$T_obs, idx = sim$tel$idx),
                       chains = 4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500,
                       seed = r, refresh = 0, show_messages = FALSE)
  dt <- ft$draws(c("H", "t0", "t1"), format = "draws_df")
  H_tel <- median(dt$H)
  T_tel_pred <- exp(median(dt$t0) + median(dt$t1) * sim$tel$z)
  g <- sim$grp
  cwt_data <- list(
    Ncwtgrp = nrow(g), Nyrs = truth$n_years, Nobs = nrow(sim$recs), UseSizeEffect = 1L,
    Sz_cwt = (g$L - mean(g$L)) / sd(g$L), CovX_cwt = rep(0, nrow(g)), ReachKM_cwt = rep(truth$D, nrow(g)),
    ObsTT = sim$recs$T, TTind = sim$recs$grp,
    Nrel = as.integer(g$Nrel), cwt_recaptures = as.integer(g$recoveries),
    Nwk = nrow(sim$wk_tbl), lp_mu = sim$wk_tbl$lp_mean, lp_sd = sim$wk_tbl$lp_sd,
    Nseg = nrow(sim$segs), seg_lo = sim$segs$lo, seg_hi = sim$segs$hi, seg_wk = sim$segs$wk,
    seg_start = sim$seg_idx$seg_start, seg_end = sim$seg_idx$seg_end, obs_wk = sim$recs$wk,
    H_tel = H_tel, Ntel = length(T_tel_pred), T_tel = T_tel_pred,
    log_Lratio_tel = log(truth$L_r / sim$tel$L), log_Lratio_cwt = log(truth$L_r / g$L),
    c_mu = 1, c_sd = 0.25, UseOffset = 1L)
  init_fun <- function() list(T_bReach = log(med_T(mean(g$L)) / 3) + rnorm(1, 0, 0.05), TT_bCov = 0, T_bSz = 0,
                              z_TT = rep(0, nrow(g)), sd_cwtTTre = 0.2, Pro_sd = 0.45,
                              z_pcap = rep(0, cwt_data$Nwk), c = 1, delta = 0,
                              z_h = rep(0, nrow(g)), sd_h = 0.2)
  fc <- mod_cwt$sample(data = cwt_data, init = init_fun, chains = 4, parallel_chains = 4,
                       iter_warmup = 500, iter_sampling = 500, adapt_delta = 0.95,
                       seed = r, refresh = 0, show_messages = FALSE)
  d <- unclass(fc$draws(c("c", "delta", "sd_h", "M_r", "surv_cwt", "p_rec_given_surv"), format = "draws_matrix"))
  d <- matrix(as.numeric(d), nrow = nrow(d), dimnames = list(NULL, colnames(d)))
  q <- \(x, p) unname(quantile(x, p))
  surv_med <- apply(d[, grep("^surv_cwt", colnames(d))], 2, median)
  surv_lo <- apply(d[, grep("^surv_cwt", colnames(d))], 2, quantile, 0.05)
  surv_hi <- apply(d[, grep("^surv_cwt", colnames(d))], 2, quantile, 0.95)
  prec_med <- apply(d[, grep("^p_rec_given_surv", colnames(d))], 2, median)
  diag <- fc$diagnostic_summary(quiet = TRUE)
  tibble(
    scenario = sc$scenario, rep = r, c_true = sc$c_true, delta_true = sc$delta_true,
    recoveries = sum(g$recoveries), tel_detected = sum(sim$tel$det),
    c_med = median(d[, "c"]), c_q05 = q(d[, "c"], 0.05), c_q95 = q(d[, "c"], 0.95),
    delta_med = median(d[, "delta"]), delta_q05 = q(d[, "delta"], 0.05), delta_q95 = q(d[, "delta"], 0.95),
    sd_h_med = median(d[, "sd_h"]), M_r_med = median(d[, "M_r"]), M_r_true = truth$M_r,
    surv_true_mean = mean(g$true_surv), surv_est_mean = mean(surv_med),
    surv_log_bias = mean(log(surv_med / g$true_surv)),
    surv_cover90 = mean(g$true_surv >= surv_lo & g$true_surv <= surv_hi),
    prec_log_bias = mean(log(prec_med / g$true_prec)),
    divergences = sum(diag$num_divergent), minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

# ---- run ------------------------------------------------------------------------------
todo <- crossing(scenarios, rep = seq_len(n_reps)) |> anti_join(done, by = c("scenario", "rep"))
message(nrow(todo), " replicates to run")
walk(seq_len(nrow(todo)), \(i) {
  sc <- todo[i, ]
  res <- tryCatch(fit_rep(sc, sc$rep), error = \(e) {
    message("Replicate failed: ", sc$scenario, " ", sc$rep, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    write.table(res, out_file, sep = ",", row.names = FALSE,
                col.names = !file.exists(out_file), append = file.exists(out_file))
    message(sprintf("%s rep %02d: c %.2f, delta %.2f (%.2f to %.2f), %d recoveries, %.1f min",
                    res$scenario, res$rep, res$c_med, res$delta_med, res$delta_q05, res$delta_q95,
                    res$recoveries, res$minutes))
  }
})
message("Done. Next: source('52_summarise_lorenzen_sims.R')")
