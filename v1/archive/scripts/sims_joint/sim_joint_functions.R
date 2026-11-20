# Functions for the joint CJS + CWT simulation (draft, 30 Sep 2026)
#
# Question: can acoustic telemetry (CJS) and Knights Landing CWT recoveries be fit
# as ONE survival and travel-time process with two observation models, and does
# that fix the CJS vs CWT flow disagreement when KL trap efficiency depends on flow?
#
# Process (shared by all fish):
#   reach survival per 100 km: logit s_j = b_j + bQ*Q + bL*L + group RE (+ delta_AT for acoustic fish)
#   reach travel time: log days per 100 km tau_j = t_j + gQ*Q + gL*L + group RE
# Observation:
#   telemetry: detections at Woodson, Butte, Sacramento, Delta (CJS), travel times
#   CWT: daily recoveries at KL = survivors arriving that day x p_KL(day)
#   p_KL(day): logit = a0 + year RE + bp * daily KL flow; efficiency trials inform it
#
# Placeholder geography (river km from the Golden Gate). Replace with 02's
# ReachKM values before using results quantitatively.

library(dplyr)
library(tidyr)
library(purrr)

rkm <- c(Battle = 522, RBDD = 462, Woodson = 395, Butte = 318,
         KL = 222, Sacramento = 172, Delta = 71)

fl_mu <- 80   # fork-length standardization used for BOTH data sources
fl_sd <- 15

# Truth ------------------------------------------------------------------------------
default_truth <- function() {
  list(
    b = c(0.3, 1.0, 0.8, 1.0),        # logit survival per 100 km, reaches 1-4
    bQ = 0.5,                         # flow effect on survival (per SD of log flow)
    bL = 0.3,                         # size effect on survival (per 15 mm)
    delta_AT = -0.3,                  # acoustic tag effect on survival (logit, per 100 km)
    sd_S = 0.5,                       # release-group RE on survival
    t = c(1.2, 1.0, 1.0, 1.1),        # log days per 100 km, reaches 1-4
    gQ = -0.15, gL = -0.15,           # flow and size effects on travel time
    sd_T = 0.2,
    sigma_t = 0.35, sigma_c = 0.45,   # individual log travel-time SD, telemetry and CWT
    mu_p = c(2.0, 1.5, 2.0, 1.0),     # telemetry detection (logit), W, B, S, D
    sd_p = rep(0.5, 4),
    a0 = qlogis(0.003),               # KL trap efficiency at mean flow (logit)
    sd_a = 0.3,                       # year RE on KL efficiency
    bp = 0,                           # flow effect on KL efficiency (scenario)
    trial_sd = 0,                     # skew-normal trial-level noise scale (0 = none)
    trial_alpha = 0,                  # skew-normal shape
    cwt_noise = FALSE                 # TRUE: noise is real weekly efficiency variation
                                      # (CWT fish experience it); FALSE: trial artifact
  )
}

# Scenarios --------------------------------------------------------------------------
# S1: KL efficiency does not depend on flow
# S2: KL efficiency falls with flow (the leading explanation for the flat CWT flow slope)
# S3, S4: KL efficiency as in the (older) BT-SPAS-X KL fit, 2 Oct: intercept -7.35,
#   efficiency RISES with flow (b_flow +0.35), year SD 0.84, right-skewed trial noise
#   (scale 1.62, shape 4.5).
#   S3: that noise is a trial artifact; CWT fish experience efficiency at noise = 0
#   S4: that noise is real week-to-week efficiency variation shared by CWT fish and trials
scenarios <- list(
  S1_const_pKL = list(bp = 0),
  S2_flow_pKL  = list(bp = -0.6),
  S3_KLup_trialnoise = list(a0 = -7.35, bp = 0.35, sd_a = 0.84,
                            trial_sd = 1.62, trial_alpha = 4.5, cwt_noise = FALSE),
  S4_KLup_realnoise  = list(a0 = -7.35, bp = 0.35, sd_a = 0.84,
                            trial_sd = 1.62, trial_alpha = 4.5, cwt_noise = TRUE),
  # S3 with a smaller true tag effect (-0.1 per 100 km: tagged survival to KL about 11% lower,
  # vs about 30% at -0.3), closer to what the literature implies at under 5% tag burden
  S3_tag01 = list(a0 = -7.35, bp = 0.35, sd_a = 0.84, delta_AT = -0.1,
                  trial_sd = 1.62, trial_alpha = 4.5, cwt_noise = FALSE)
)

# skew-normal draws, location 0 (as the BT-SPAS-X process deviates)
rsn <- function(n, omega, alpha) {
  d <- alpha / sqrt(1 + alpha^2)
  omega * (d * abs(rnorm(n)) + sqrt(1 - d^2) * rnorm(n))
}

default_design <- function() {
  list(
    cwt_years = setdiff(2000:2024, c(2013, 2015, 2022)),  # 22 years, as in the filtered data
    tel_years = 2013:2024,
    n_tel_groups = 25,
    n_tel_fish = 150,            # per group (real data: about 360; reduced for speed)
    n_cwt_groups = 300,          # real data: 394
    cwt_N_median = 60000,        # released per CWT group (placeholder)
    battle_share_tel = 0.6,
    battle_share_cwt = 0.94,     # 372 of 394 in the real data
    cwt_size_confounded = TRUE,  # RBDD 43-60 mm, Battle 61-87 mm, as in the real data
    n_days = 260,                # daily flow series length (day 1 = 1 Dec)
    D = 60,                      # CWT recovery window (days after release)
    Q_window = 30,               # days after release averaged for the group flow covariate
    trials_per_year = 5,
    trial_m = 1500               # marked fish per efficiency trial
  )
}

# Flows ------------------------------------------------------------------------------
# log flow = base + year wetness + seasonal hump + AR(1) storms. Standardized on the
# log scale over all years and days (q). No loops: AR(1) via stats::filter.
sim_flows <- function(years, n_days) {
  wet <- rnorm(length(years))
  season <- 0.4 * sin(pi * seq_len(n_days) / n_days)
  map2(years, wet, \(y, w) {
    ar <- as.numeric(stats::filter(rnorm(n_days, 0, 0.15), 0.9, method = "recursive"))
    tibble(year = y, day = seq_len(n_days), logQ = log(10000) + 0.5 * w + season + ar)
  }) |>
    list_rbind() |>
    mutate(q = (logQ - mean(logQ)) / sd(logQ), cfs = exp(logQ))
}

# Release groups ---------------------------------------------------------------------
make_groups <- function(des, flows) {
  q_mat <- flows |>
    select(year, day, q) |>
    pivot_wider(names_from = day, values_from = q)
  yrs <- q_mat$year
  q_mat <- as.matrix(q_mat[, -1])
  rownames(q_mat) <- yrs

  group_Q <- function(year, rel_day) {
    mean(q_mat[as.character(year), rel_day:(rel_day + des$Q_window - 1)])
  }

  tel_years <- c(des$tel_years,
                 sample(des$tel_years, des$n_tel_groups - length(des$tel_years), TRUE))
  tel <- tibble(
    source = "tel",
    year = tel_years,
    site = if_else(runif(des$n_tel_groups) < des$battle_share_tel, "Battle", "RBDD"),
    rel_day = sample(20:120, des$n_tel_groups, TRUE),
    FL = pmin(pmax(rnorm(des$n_tel_groups, 90, 10), 70), 130),
    N = des$n_tel_fish
  )

  cwt_years <- c(des$cwt_years,
                 sample(des$cwt_years, des$n_cwt_groups - length(des$cwt_years), TRUE))
  cwt_site <- if_else(runif(des$n_cwt_groups) < des$battle_share_cwt, "Battle", "RBDD")
  cwt <- tibble(
    source = "cwt",
    year = cwt_years,
    site = cwt_site,
    rel_day = sample(20:150, des$n_cwt_groups, TRUE),
    FL = if (des$cwt_size_confounded) {
      if_else(site == "RBDD", runif(des$n_cwt_groups, 43, 60), runif(des$n_cwt_groups, 61, 87))
    } else {
      runif(des$n_cwt_groups, 43, 87)
    },
    N = round(rlnorm(des$n_cwt_groups, log(des$cwt_N_median), 0.5))
  )

  groups <- bind_rows(tel, cwt) |>
    mutate(
      Q = map2_dbl(year, rel_day, group_Q),
      L = (FL - fl_mu) / fl_sd,
      km1 = rkm[site] - rkm["Woodson"],
      km2 = rkm["Woodson"] - rkm["Butte"],
      km3 = if_else(source == "tel", rkm["Butte"] - rkm["Sacramento"], rkm["Butte"] - rkm["KL"]),
      km4 = if_else(source == "tel", rkm["Sacramento"] - rkm["Delta"], 0)
    )
  list(groups = groups, q_mat = q_mat)
}

# Process ----------------------------------------------------------------------------
add_truth_to_groups <- function(groups, tr) {
  n <- nrow(groups)
  groups |>
    mutate(
      zS = rnorm(n), zT = rnorm(n),
      # survival per 100 km and over each reach
      eta = pmap(list(Q, L, zS, source), \(Q, L, zS, src)
                 tr$b + tr$bQ * Q + tr$bL * L + tr$sd_S * zS + (src == "tel") * tr$delta_AT),
      km = pmap(list(km1, km2, km3, km4), c),
      phi = map2(eta, km, \(e, k) plogis(e)^(k / 100)),
      # cumulative median travel time to the end of each reach
      tau = pmap(list(Q, L, zT), \(Q, L, zT) tr$t + tr$gQ * Q + tr$gL * L + tr$sd_T * zT),
      cum_tt = map2(tau, km, \(ta, k) cumsum(exp(ta) * k / 100))
    )
}

sim_telemetry <- function(g, p_year) {
  # g: one telemetry group (row, as list); p_year: 4-vector of detection probs
  n <- g$N
  surv <- matrix(runif(n * 4), n, 4) < matrix(g$phi, n, 4, byrow = TRUE)
  alive <- t(apply(surv, 1, cumprod))
  det <- alive * (matrix(runif(n * 4), n, 4) < matrix(p_year, n, 4, byrow = TRUE))
  h <- as.integer(det %*% c(8, 4, 2, 1)) + 1L
  # travel times: individual multiplier on the cumulative median (marginally lognormal)
  eps <- rnorm(n, 0, g$sigma_t)
  logt <- outer(eps, log(g$cum_tt[1:3]), `+`)
  tt_det <- det[, 1:3] == 1
  list(
    n_hist = tabulate(h, 16),
    tt_n = colSums(tt_det),
    tt_s1 = colSums(logt * tt_det),
    tt_s2 = colSums(logt^2 * tt_det)
  )
}

sim_cwt <- function(g, q_row, a0, za_y, bp, D, dev_row = rep(0, length(q_row))) {
  # survivors to KL, their arrival days (multinomial over days 1..D and beyond),
  # then binomial capture with the daily KL efficiency
  S_KL <- prod(g$phi[1:3])
  n_surv <- rbinom(1, g$N, S_KL)
  mu <- log(g$cum_tt[3])
  F <- pnorm((log(seq_len(D)) - mu) / g$sigma_c)
  parr <- c(diff(c(0, F)), 1 - F[D])
  n_arr <- as.vector(rmultinom(1, n_surv, parr))[seq_len(D)]
  q_d <- q_row[g$rel_day + seq_len(D)]
  p_d <- plogis(a0 + za_y + bp * q_d + dev_row[g$rel_day + seq_len(D)])
  list(y = rbinom(D, n_arr, p_d), q = q_d)
}

# One replicate ----------------------------------------------------------------------
simulate_replicate <- function(des = default_design(), tr = default_truth()) {
  all_years <- sort(union(des$cwt_years, des$tel_years))
  flows <- sim_flows(all_years, des$n_days)
  gg <- make_groups(des, flows)
  groups <- add_truth_to_groups(gg$groups, tr) |>
    mutate(yr = match(year, all_years))
  q_mat <- gg$q_mat

  Y <- length(all_years)
  p_det <- plogis(matrix(tr$mu_p, Y, 4, byrow = TRUE) +
                    matrix(rnorm(Y * 4), Y, 4) * matrix(tr$sd_p, Y, 4, byrow = TRUE))
  za <- rnorm(Y, 0, tr$sd_a)
  # weekly KL efficiency deviates (year x week), expanded to days
  n_wk <- ceiling(des$n_days / 7)
  dev_wk <- matrix(rsn(Y * n_wk, tr$trial_sd, tr$trial_alpha), Y, n_wk)
  dev_day <- dev_wk[, ceiling(seq_len(des$n_days) / 7), drop = FALSE] * tr$cwt_noise

  tel <- groups |> filter(source == "tel")
  cwt <- groups |> filter(source == "cwt")

  tel_obs <- tel |>
    mutate(sigma_t = tr$sigma_t) |>
    pmap(list) |>
    map2(tel$yr, \(g, y) sim_telemetry(g, p_det[y, ]))

  cwt_obs <- cwt |>
    mutate(sigma_c = tr$sigma_c) |>
    pmap(list) |>
    map(\(g) sim_cwt(g, q_mat[as.character(g$year), ], tr$a0, za[g$yr], tr$bp, des$D,
                     dev_day[g$yr, ]))

  # efficiency trials in every year
  trials <- tibble(year = rep(all_years, each = des$trials_per_year)) |>
    mutate(yr = match(year, all_years),
           day = sample(20:180, n(), TRUE),
           q = map2_dbl(year, day, \(y, d) q_mat[as.character(y), d]),
           m = des$trial_m,
           # trial noise: the week's real deviate (S4) or an independent artifact (S3)
           dev = if (tr$cwt_noise) {
             map2_dbl(yr, day, \(y, d) dev_wk[y, ceiling(d / 7)])
           } else {
             rsn(n(), tr$trial_sd, tr$trial_alpha)
           },
           r = rbinom(n(), m, plogis(tr$a0 + za[yr] + tr$bp * q + dev)))

  list(des = des, truth = tr, all_years = all_years, flows = flows,
       tel = tel, cwt = cwt, tel_obs = tel_obs, cwt_obs = cwt_obs, trials = trials)
}

# True estimands ---------------------------------------------------------------------
estimand_settings <- function() {
  list(km_est = unname(c(rkm["Battle"] - rkm["Woodson"], rkm["Woodson"] - rkm["Butte"],
                         rkm["Butte"] - rkm["KL"])),
       L_est = (75 - fl_mu) / fl_sd)
}

true_values <- function(tr) {
  es <- estimand_settings()
  tibble(
    parameter = c("bQ", "bL", "delta_AT", "gQ", "gL", "sd_S", "a0", "bp", "sdk",
                  "S_KL", "S_KL_tagged", "TT_KL", "p_KL_mean"),
    truth = c(tr$bQ, tr$bL, tr$delta_AT, tr$gQ, tr$gL, tr$sd_S, tr$a0, tr$bp, tr$trial_sd,
              prod(plogis(tr$b[1:3] + tr$bL * es$L_est)^(es$km_est / 100)),
              prod(plogis(tr$b[1:3] + tr$bL * es$L_est + tr$delta_AT)^(es$km_est / 100)),
              sum(exp(tr$t[1:3] + tr$gL * es$L_est) * es$km_est / 100),
              plogis(tr$a0))
  )
}

# Models -----------------------------------------------------------------------------
model_flags <- list(
  A_cwt_only_constp = list(use_cjs = 0, use_cwt = 1, use_trials = 0, p_flow = 0, p_year = 0, est_delta = 0, pool_reaches = 1, trial_re = 0, cwt_noise_mean = 0),
  B_joint_constp    = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 0, p_year = 1, est_delta = 1, pool_reaches = 0, trial_re = 0, cwt_noise_mean = 0),
  C_joint_flowp     = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 1, p_year = 1, est_delta = 1, pool_reaches = 0, trial_re = 0, cwt_noise_mean = 0),
  D_cjs_only        = list(use_cjs = 1, use_cwt = 0, use_trials = 0, p_flow = 0, p_year = 0, est_delta = 0, pool_reaches = 0, trial_re = 0, cwt_noise_mean = 0),
  # E, F: joint, flow pCap, skew-normal trial noise (BT-SPAS-X KL structure)
  E_joint_trialnoise_artifact = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 1, p_year = 1, est_delta = 1, pool_reaches = 0, trial_re = 1, cwt_noise_mean = 0),
  F_joint_trialnoise_real     = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 1, p_year = 1, est_delta = 1, pool_reaches = 0, trial_re = 1, cwt_noise_mean = 1),
  # G: as E and F but tag effect fixed at 0 (exchangeability, as in Perry et al. 2010, 2018;
  #    Hance et al. 2021): telemetry anchors the survival level for all fish. The sim truth
  #    has delta_AT = -0.3, so G also shows the bias from assuming 0 under a modest real effect.
  G_tag0_noise_artifact = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 1, p_year = 1, est_delta = 0, pool_reaches = 0, trial_re = 1, cwt_noise_mean = 0),
  G_tag0_noise_real     = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 1, p_year = 1, est_delta = 0, pool_reaches = 0, trial_re = 1, cwt_noise_mean = 1),
  # E2 (3 Oct): as E with a literature-informed prior on the tag effect, delta_AT ~ normal(0, 0.15)
  #    ("soft exchangeability": 95% prior mass within +-0.3 per 100 km). Other models use delta_sd = 1.
  E2_tagprior_artifact  = list(use_cjs = 1, use_cwt = 1, use_trials = 1, p_flow = 1, p_year = 1, est_delta = 1, pool_reaches = 0, trial_re = 1, cwt_noise_mean = 0, delta_sd = 0.15)
)

# Models fitted in each scenario
scenario_models <- list(
  S1_const_pKL = c("A_cwt_only_constp", "B_joint_constp", "C_joint_flowp", "D_cjs_only"),
  S2_flow_pKL  = c("A_cwt_only_constp", "B_joint_constp", "C_joint_flowp", "D_cjs_only"),
  S3_KLup_trialnoise = c("A_cwt_only_constp", "C_joint_flowp",
                         "E_joint_trialnoise_artifact", "F_joint_trialnoise_real",
                         "G_tag0_noise_artifact", "G_tag0_noise_real",
                         "E2_tagprior_artifact"),
  S4_KLup_realnoise  = c("A_cwt_only_constp", "C_joint_flowp",
                         "E_joint_trialnoise_artifact", "F_joint_trialnoise_real",
                         "G_tag0_noise_artifact", "G_tag0_noise_real"),
  S3_tag01 = c("E_joint_trialnoise_artifact", "E2_tagprior_artifact", "G_tag0_noise_artifact")
)

make_stan_data <- function(rep, flags, prior_a0_sd = 1.39, bin_days = 5) {
  flags <- modifyList(list(delta_sd = 1), flags)   # tag-effect prior SD; 1 unless the model sets it
  tel <- rep$tel
  cwt <- rep$cwt
  stack <- \(x, nm) do.call(rbind, map(x, nm))
  # daily CWT recoveries -> time bins (sum of recoveries, mean of daily flow)
  D_days <- rep$des$D
  bin <- ceiling(seq_len(D_days) / bin_days)
  D <- max(bin)
  bin_mat <- outer(bin, seq_len(D), `==`) * 1          # days x bins
  y_bin <- stack(rep$cwt_obs, "y") %*% bin_mat
  q_bin <- unname(stack(rep$cwt_obs, "q")) %*% sweep(bin_mat, 2, colSums(bin_mat), `/`)
  c(
    flags,
    list(
      # prior centred on the mean trial-level logit efficiency (as 06's prior is built from
      # BT-SPAS-X trial estimates): truth a0 plus the mean of the skewed trial noise
      prior_a0_mu = with(rep$truth, a0 + trial_sd * trial_alpha / sqrt(1 + trial_alpha^2) * sqrt(2 / pi)),
      prior_a0_sd = prior_a0_sd,
      Y = length(rep$all_years),
      G_t = nrow(tel), yr_t = tel$yr, Q_t = tel$Q, L_t = tel$L,
      km_t = as.matrix(tel[, c("km1", "km2", "km3", "km4")]),
      n_hist = stack(rep$tel_obs, "n_hist"),
      tt_n = matrix(as.integer(stack(rep$tel_obs, "tt_n")), nrow(tel)),
      tt_s1 = stack(rep$tel_obs, "tt_s1"),
      tt_s2 = stack(rep$tel_obs, "tt_s2"),
      G_c = nrow(cwt), D = D, edge = as.numeric(tapply(seq_len(D_days), bin, max)),
      yr_c = cwt$yr, Q_c = cwt$Q, L_c = cwt$L,
      km_c = as.matrix(cwt[, c("km1", "km2", "km3")]),
      N_c = as.numeric(cwt$N),
      y_c = matrix(as.integer(round(y_bin)), nrow(cwt)),
      q_c = q_bin,
      K = nrow(rep$trials), yr_k = rep$trials$yr, q_k = rep$trials$q,
      m_k = as.integer(rep$trials$m), r_k = as.integer(rep$trials$r)
    ),
    estimand_settings()
  )
}

# Fitting ----------------------------------------------------------------------------
# backend "cmdstanr" (default; works on the laptop) or "rstan" (used for testing here)
fit_model <- function(stan_data, model, backend = c("cmdstanr", "rstan"),
                      chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                      adapt_delta = 0.9, seed = 1) {
  backend <- match.arg(backend)
  pars <- c("bQ", "bL", "delta_AT", "gQ", "gL", "sd_S", "a0", "bp", "sdk",
            "S_KL", "S_KL_tagged", "TT_KL", "p_KL_mean")
  t0 <- Sys.time()
  if (backend == "cmdstanr") {
    fit <- model$sample(data = stan_data, chains = chains, parallel_chains = chains,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta, seed = seed, refresh = 0,
                        show_messages = FALSE)
    draws <- posterior::as_draws_df(fit$draws(pars))
    sd_ <- fit$sampler_diagnostics(format = "df")
    n_div <- sum(sd_$divergent__)
    n_td <- sum(sd_$treedepth__ >= 10)
  } else {
    fit <- rstan::sampling(model, data = stan_data, chains = chains, cores = chains,
                           warmup = iter_warmup, iter = iter_warmup + iter_sampling,
                           control = list(adapt_delta = adapt_delta), seed = seed,
                           refresh = 0, pars = pars)
    draws <- posterior::as_draws_df(rstan::extract(fit, pars = pars, permuted = FALSE))
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    n_div <- sum(map_dbl(sp, \(x) sum(x[, "divergent__"])))
    n_td <- sum(map_dbl(sp, \(x) sum(x[, "treedepth__"] >= 10)))
  }
  summ <- posterior::summarise_draws(
    posterior::subset_draws(draws, variable = pars),
    mean, median, sd,
    q2.5 = \(x) unname(quantile(x, 0.025)), q97.5 = \(x) unname(quantile(x, 0.975)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
  ) |>
    rename(parameter = variable)
  list(summary = summ, n_div = n_div, n_treedepth = n_td,
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
