# ============================================================
# 57_run_lorenzen_gamma_realsize.R
# Same as 55, but with the REAL number of CWT groups: 18 years x 20 groups
# = 360 groups (real data: 357). Question: with as many groups as the real
# data, (1) is gamma = 0 recovered (the real fit gave 0.02, 0 to 0.08) and
# gamma = 1 not mistaken for 0, and (2) does fitting per day to a gamma = 0
# truth give delta and c like the real per-day fit (delta -0.47, c 0.28)?
# 2 truths (gamma 0 and 1) x 10 replicates x 2 fits. Each CWT fit has about
# 9 times more groups than in 55, so expect roughly 15 to 30 min per fit,
# about 10 to 20 hours in total. Restartable. Then source("58_summarise_gamma_realsize.R").
# ============================================================

library(tidyverse)
library(cmdstanr)

n_reps <- 10
truths <- tibble(truth = c("gamma_0", "gamma_1"), gamma_true = c(0, 1))
truth <- list(M_r = 0.13, L_r = 70, c = 1, D = 300, sigma_T = 0.45, sd_tt_grp = 0.2, sd_h = 0.2,
              p_tel = 0.9, n_tel = 2500, n_years = 18, grp_per_year = 20, small_share = 0.15,
              rel_weeks = 8:20, last_week = 31, L_typ = 88)
med_T <- function(L) 10 * (100 / L)^1.5            # median days to KL (300 km): 12 d at 88 mm, 33 d at 50 mm

dir.create("results", showWarnings = FALSE)
out_file <- "results/reps_gamma_realsize.csv"
done <- if (file.exists(out_file)) distinct(read.csv(out_file), truth, rep) else tibble(truth = character(), rep = integer())

kl_weekly <- read.csv("inputs/kl_weekly_pcap.csv")
kl_draws <- readRDS("inputs/kl_pcap_draws.rds")
mod_tel <- cmdstan_model("sim_tel_lorenzen.stan")
mod_cwt <- cmdstan_model("CWT_oneway_gamma.stan")

# ---- simulate one replicate ------------------------------------------------------------
simulate_rep <- function(g_true, seed) {
  set.seed(seed)
  # rate chosen so a typical telemetry fish has the same survival for any gamma
  M_g <- truth$M_r * med_T(truth$L_typ)^(1 - g_true)
  surv_fish <- function(L, T, extra = 0) exp(-M_g * (truth$L_r / L)^truth$c * exp(extra) * T^g_true)

  # telemetry: realistic sizes
  L_tel <- pmin(pmax(rnorm(truth$n_tel, 88, 7), 75), 115)
  T_tel <- rlnorm(truth$n_tel, log(med_T(L_tel)), truth$sigma_T)
  det <- rbinom(truth$n_tel, 1, rbinom(truth$n_tel, 1, surv_fish(L_tel, T_tel)) * truth$p_tel)
  z_tel <- (L_tel - mean(L_tel)) / sd(L_tel)

  # CWT: mostly 65 to 80 mm groups, about 15% small (45 to 60 mm) groups
  yrs <- sample(names(kl_draws), truth$n_years)
  ng <- truth$n_years * truth$grp_per_year
  grp <- tibble(year = rep(seq_len(truth$n_years), each = truth$grp_per_year)) |>
    mutate(kl_year = yrs[year],
           draw = rep(sample(200, truth$n_years, replace = TRUE), each = truth$grp_per_year),
           rel_week = sample(truth$rel_weeks, n(), replace = TRUE),
           Nrel = round(runif(n(), 7000, 27000)),   # scaled so total recoveries are near the real 3,013
           small = runif(n()) < truth$small_share,
           L = if_else(small, runif(n(), 45, 60), runif(n(), 65, 80)),
           tt_re = rnorm(n(), 0, truth$sd_tt_grp),
           h_re = rnorm(n(), 0, truth$sd_h),
           grp = row_number())
  sim_grp <- function(g) {
    T <- rlnorm(g$Nrel, log(med_T(g$L)) + g$tt_re, truth$sigma_T)
    s <- surv_fish(g$L, T, g$h_re)
    wk <- g$rel_week + floor(T / 7)
    p_week <- kl_draws[[g$kl_year]][g$draw, ]
    p_fish <- ifelse(wk <= truth$last_week, p_week[pmin(wk, 31)], 0)
    cap <- rbinom(g$Nrel, 1, s * p_fish)
    list(true_surv = mean(s), true_prec = mean(p_fish * s) / mean(s),
         rec = tibble(grp = g$grp, T = T[cap == 1], run_week = wk[cap == 1]))
  }
  sims <- map(seq_len(nrow(grp)), \(i) sim_grp(grp[i, ]))
  grp <- grp |> mutate(true_surv = map_dbl(sims, "true_surv"), true_prec = map_dbl(sims, "true_prec"),
                       recoveries = map_int(sims, \(x) nrow(x$rec)))
  recs <- map(sims, "rec") |> list_rbind()
  segs <- grp |> select(grp, kl_year, rel_week) |>
    mutate(run_week = map(rel_week, \(w) w:truth$last_week)) |> unnest(run_week) |>
    mutate(lo = (run_week - rel_week) * 7, hi = lo + 7)
  wk_tbl <- segs |> distinct(kl_year, run_week) |>
    left_join(kl_weekly |> select(fit, week, lp_mean, lp_sd), by = c("kl_year" = "fit", "run_week" = "week")) |>
    mutate(wk = row_number())
  segs <- segs |> left_join(wk_tbl |> select(kl_year, run_week, wk), by = c("kl_year", "run_week")) |>
    arrange(grp, lo) |> mutate(seg = row_number())
  seg_idx <- segs |> group_by(grp) |> summarise(seg_start = min(seg), seg_end = max(seg))
  recs <- recs |> left_join(grp |> select(grp, kl_year), by = "grp") |>
    left_join(wk_tbl |> select(kl_year, run_week, wk), by = c("kl_year", "run_week")) |> arrange(grp, T)
  list(grp = grp, recs = recs, segs = segs, wk_tbl = wk_tbl, seg_idx = seg_idx,
       tel = list(z = z_tel, L = L_tel, det = det, T_obs = T_tel[det == 1], idx = which(det == 1)))
}

# ---- fit one replicate (one telemetry fit, two CWT fits) --------------------------------
fit_rep <- function(tr, r) {
  seed <- 19000 + 1000 * match(tr$truth, truths$truth) + r
  sim <- simulate_rep(tr$gamma_true, seed)
  ft <- mod_tel$sample(data = list(N = length(sim$tel$z), z = sim$tel$z, det = sim$tel$det,
                                   p = truth$p_tel, Nobs = length(sim$tel$idx),
                                   T_obs = sim$tel$T_obs, idx = sim$tel$idx),
                       chains = 4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500,
                       seed = r, refresh = 0, show_messages = FALSE)
  dt <- ft$draws(c("H", "t0", "t1"), format = "draws_df")
  g <- sim$grp
  base <- list(
    Ncwtgrp = nrow(g), Nyrs = truth$n_years, Nobs = nrow(sim$recs), UseSizeEffect = 1L,
    Sz_cwt = (g$L - mean(g$L)) / sd(g$L), CovX_cwt = rep(0, nrow(g)), ReachKM_cwt = rep(truth$D, nrow(g)),
    ObsTT = sim$recs$T, TTind = sim$recs$grp,
    Nrel = as.integer(g$Nrel), cwt_recaptures = as.integer(g$recoveries),
    Nwk = nrow(sim$wk_tbl), lp_mu = sim$wk_tbl$lp_mean, lp_sd = sim$wk_tbl$lp_sd,
    Nseg = nrow(sim$segs), seg_lo = sim$segs$lo, seg_hi = sim$segs$hi, seg_wk = sim$segs$wk,
    seg_start = sim$seg_idx$seg_start, seg_end = sim$seg_idx$seg_end, obs_wk = sim$recs$wk,
    H_tel = median(dt$H), Ntel = length(sim$tel$z), T_tel = exp(median(dt$t0) + median(dt$t1) * sim$tel$z),
    log_Lratio_tel = log(truth$L_r / sim$tel$L), log_Lratio_cwt = log(truth$L_r / g$L),
    c_mu = 1, c_sd = 0.25, UseOffset = 1L, gamma_mu = 1, gamma_sd = 0.5)
  init_fun <- function() list(T_bReach = log(med_T(mean(g$L)) / 3) + rnorm(1, 0, 0.05), TT_bCov = 0, T_bSz = 0,
                              z_TT = rep(0, nrow(g)), sd_cwtTTre = 0.2, Pro_sd = 0.45,
                              z_pcap = rep(0, base$Nwk), c = 1, delta = 0,
                              z_h = rep(0, nrow(g)), sd_h = 0.2, gamma_raw = 1)
  fit_one <- function(use_gamma, label) {
    t0 <- Sys.time()
    fc <- mod_cwt$sample(data = c(base, list(UseGamma = as.integer(use_gamma))), init = init_fun,
                         chains = 4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500,
                         adapt_delta = 0.95, seed = r, refresh = 0, show_messages = FALSE)
    d <- unclass(fc$draws(c("c", "delta", "gamma_raw", "surv_cwt", "p_rec_given_surv"), format = "draws_matrix"))
    d <- matrix(as.numeric(d), nrow = nrow(d), dimnames = list(NULL, colnames(d)))
    q <- \(x, p) unname(quantile(x, p))
    sv <- d[, grep("^surv_cwt", colnames(d)), drop = FALSE]
    pr <- d[, grep("^p_rec_given_surv", colnames(d)), drop = FALSE]
    s_med <- apply(sv, 2, median); s_lo <- apply(sv, 2, quantile, 0.05); s_hi <- apply(sv, 2, quantile, 0.95)
    err <- log(s_med / g$true_surv)
    gam <- if (use_gamma) d[, "gamma_raw"] else rep(1, nrow(d))
    tibble(truth = tr$truth, gamma_true = tr$gamma_true, fit = label, rep = r,
           recoveries = sum(g$recoveries),
           gamma_med = median(gam), gamma_q05 = q(gam, 0.05), gamma_q95 = q(gam, 0.95),
           c_med = median(d[, "c"]), c_q05 = q(d[, "c"], 0.05), c_q95 = q(d[, "c"], 0.95),
           delta_med = median(d[, "delta"]), delta_q05 = q(d[, "delta"], 0.05), delta_q95 = q(d[, "delta"], 0.95),
           surv_log_bias = mean(err),
           surv_log_bias_small = mean(err[g$small]), surv_log_bias_large = mean(err[!g$small]),
           surv_cover90 = mean(g$true_surv >= s_lo & g$true_surv <= s_hi),
           prec_log_bias = mean(log(apply(pr, 2, median) / g$true_prec)),
           divergences = sum(fc$diagnostic_summary(quiet = TRUE)$num_divergent),
           minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  bind_rows(fit_one(FALSE, "per_day"), fit_one(TRUE, "gamma"))
}

# ---- run ------------------------------------------------------------------------------
todo <- crossing(truths, rep = seq_len(n_reps)) |> anti_join(done, by = c("truth", "rep"))
message(nrow(todo), " replicates to run (each = 1 telemetry fit + 2 CWT fits)")
walk(seq_len(nrow(todo)), \(i) {
  tr <- todo[i, ]
  res <- tryCatch(fit_rep(tr, tr$rep), error = \(e) {
    message("Replicate failed: ", tr$truth, " ", tr$rep, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    write.table(res, out_file, sep = ",", row.names = FALSE,
                col.names = !file.exists(out_file), append = file.exists(out_file))
    walk(seq_len(nrow(res)), \(k) message(sprintf("%s rep %02d fit %-7s: gamma %.2f, delta %.2f (%.2f to %.2f), c %.2f, surv bias %+.2f",
         res$truth[k], res$rep[k], res$fit[k], res$gamma_med[k], res$delta_med[k], res$delta_q05[k],
         res$delta_q95[k], res$c_med[k], res$surv_log_bias[k])))
  }
})
message("Done. Next: source('58_summarise_gamma_realsize.R')")
