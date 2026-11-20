# sim_functions.R
options(cmdstanr_warn_inits = FALSE)
# Simulation-recovery functions for the production CJS model (CovIndCont.stan)
# and its non-centered revision (CovIndCont_nc.stan). Sourced by
# 10_sim_recovery_cjs.R and 11_summarise_sim_recovery.R.
#
# The simulator mirrors CovIndCont.stan exactly: mainstem fish with 4 reaches
# and 5 detection locations, Butte/Feather fish with 2 reaches and 3 locations,
# year-by-reach detection, release-group-by-reach random effects on survival and
# travel time, and lognormal cumulative travel times observed at each detection.

library(purrr)
library(dplyr)
library(tibble)

# ---------------------------------------------------------------------------
# Truth: posterior means from the corrected 03 fit (CovIndCont run 15,
# flow + fork length, 25 Sep 2026). sdPb is set to 1 on purpose: the fitted
# sdPb (2.4 to 3.7, flat prior) implies years with detection near 0 or 1,
# which is not a useful test design.
# ---------------------------------------------------------------------------
truth_03 <- list(
  S_bReach = c(0.44, 0.04, 0.18, -0.08),
  T_bReach = c(0.71, 0.55, 0.53, 0.93),
  S_bCov = 0.70, TT_bCov = -0.15, S_bSz = 0.10, T_bSz = -0.025,
  S_bTrib = c(-1.42, -2.11), T_bTrib = c(1.39, 1.44),
  S_bCovT = c(2.11, 2.34), TT_bCovT = c(0.80, -0.06),
  RE_sd = c(2.02, 1.01, 0.93, 1.26), RE_sdT = c(0.80, 1.10),
  TTRE_sd = c(0.51, 0.30, 0.31, 0.09), TTRE_sdT = c(0.46, 0.57),
  Pro_sd = 0.63, Pro_sdT = 0.51,
  muPb = c(3.70, 2.52, 4.76, 1.23), sdPb = c(1, 1, 1, 1)
)

tracked_pars <- c(
  "S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz",
  "S_bTrib", "T_bTrib", "S_bCovT", "TT_bCovT",
  "RE_sd", "RE_sdT", "TTRE_sd", "TTRE_sdT", "Pro_sd", "Pro_sdT",
  "muPb", "sdPb"
)

# ---------------------------------------------------------------------------
# Design: the real fish, release groups, years, sizes, flows and distances,
# taken from the objects 02_GetData.R builds (same settings as 03 run 15).
# ---------------------------------------------------------------------------
make_design <- function() {
  musz <- mean(c(FL, FL_T), na.rm = TRUE)
  sdsz <- sd(c(FL, FL_T), na.rm = TRUE)
  list(
    main = list(
      N = Nind, rg = as.integer(rgind), yr = as.integer(yrind),
      Sz = (FL - musz) / sdsz,
      CovX = unname(as.matrix(cbind(MaxflowSac.z, MaxflowSac.z, MaxflowSac.z, MaxflowDelta.z))),
      Rmult = unname(as.matrix(Rmult)),
      KM = unname(as.matrix(ReachKM_ind))
    ),
    trib = list(
      N = NindT, rg = as.integer(rgindT), yr = as.integer(yrindT), trib = as.integer(trib_ind),
      Sz = (FL_T - musz) / sdsz,
      CovX = unname(as.matrix(cbind(MaxflowT.z, MaxflowDeltaT.z))),
      Rmult = as.numeric(RmultT),
      KM = unname(as.matrix(ReachKMT_ind))
    ),
    Nyrs = as.integer(Nyrs), Nrg = as.integer(Nrg), NrgT = as.integer(NrgT)
  )
}

# Keep a fraction of fish within every release group (at least min_n per group,
# or all of them if the group is smaller), so every group, year and flow value
# in the real design stays in the simulation.
thin_part <- function(part, frac, min_n) {
  keep <- split(seq_len(part$N), part$rg) |>
    map(\(ix) if (length(ix) <= min_n) ix else sample(ix, max(min_n, round(frac * length(ix))))) |>
    unlist(use.names = FALSE) |>
    sort()
  c(list(N = length(keep)),
    map(part[setdiff(names(part), "N")],
        \(x) if (is.matrix(x)) x[keep, , drop = FALSE] else x[keep]))
}

thin_design <- function(design, frac, min_n = 20) {
  design$main <- thin_part(design$main, frac, min_n)
  design$trib <- thin_part(design$trib, frac, min_n)
  design
}

# ---------------------------------------------------------------------------
# Random effects for one replicate
# ---------------------------------------------------------------------------
draw_effects <- function(tr, design) {
  re_mat <- function(n, sds) matrix(rnorm(n * length(sds), 0, rep(sds, each = n)), n, length(sds))
  c(tr, list(
    P_b = matrix(rnorm(design$Nyrs * 4, rep(tr$muPb, each = design$Nyrs),
                       rep(tr$sdPb, each = design$Nyrs)), design$Nyrs, 4),
    S_RE = re_mat(design$Nrg, tr$RE_sd),
    TT_RE = re_mat(design$Nrg, tr$TTRE_sd),
    S_REt = re_mat(design$NrgT, tr$RE_sdT),
    TT_RET = re_mat(design$NrgT, tr$TTRE_sdT)
  ))
}

# ---------------------------------------------------------------------------
# Simulate capture histories and travel times
# ---------------------------------------------------------------------------
sim_mainstem <- function(tr, d) {
  n <- d$N
  surv <- plogis(sweep(tr$S_RE[d$rg, , drop = FALSE], 2, tr$S_bReach, "+") +
                   tr$S_bCov * d$CovX + tr$S_bSz * d$Sz)^d$Rmult
  pcap <- plogis(tr$P_b[d$yr, , drop = FALSE]) # detection at locations 2 to 5
  alive <- accumulate(1:4, \(a, j) a * rbinom(n, 1, surv[, j]), .init = rep(1L, n)) |>
    do.call(what = cbind)
  det <- alive[, 2:5] * matrix(rbinom(n * 4, 1, pcap), n, 4)
  tt_reach <- exp(sweep(tr$TT_RE[d$rg, , drop = FALSE], 2, tr$T_bReach, "+") +
                    tr$T_bSz * d$Sz + tr$TT_bCov * d$CovX) * d$KM / 100
  cum <- t(apply(tt_reach, 1, cumsum))
  obs <- matrix(rlnorm(n * 4, log(cum), tr$Pro_sd), n, 4)
  obs[det == 0] <- NA
  CH <- cbind(1L, det)
  storage.mode(CH) <- "integer"
  list(CH = CH, TT = obs, true_surv_rel_sac = rowMeans(surv[, 1:3]))
}

sim_trib <- function(tr, d, rmult_delta) {
  n <- d$N
  s1 <- plogis(tr$S_bTrib[d$trib] + tr$S_bCovT[d$trib] * d$CovX[, 1] + tr$S_bSz * d$Sz +
                 tr$S_REt[cbind(d$rg, 1)])^d$Rmult
  s2 <- plogis(tr$S_bReach[4] + tr$S_bCov * d$CovX[, 2] + tr$S_bSz * d$Sz +
                 tr$S_REt[cbind(d$rg, 2)])^rmult_delta
  a2 <- rbinom(n, 1, s1)
  a3 <- a2 * rbinom(n, 1, s2)
  det <- cbind(a2 * rbinom(n, 1, plogis(tr$P_b[cbind(d$yr, 3)])),
               a3 * rbinom(n, 1, plogis(tr$P_b[cbind(d$yr, 4)])))
  t1 <- exp(tr$T_bTrib[d$trib] + tr$T_bSz * d$Sz + tr$TT_bCovT[d$trib] * d$CovX[, 1] +
              tr$TT_RET[cbind(d$rg, 1)]) * d$KM[, 1] / 100
  t2 <- t1 + exp(tr$T_bReach[4] + tr$T_bSz * d$Sz + tr$TT_bCov * d$CovX[, 2] +
                   tr$TT_RET[cbind(d$rg, 2)]) * d$KM[, 2] / 100
  obs <- cbind(rlnorm(n, log(t1), tr$Pro_sdT), rlnorm(n, log(t2), tr$Pro_sdT))
  obs[det == 0] <- NA
  CH <- cbind(1L, det)
  storage.mode(CH) <- "integer"
  list(CH = CH, TT = obs)
}

# Same row-wise indexing as 02_GetData.R: ObsTT holds the observed times in
# row order; TTind[i, j] points to fish i's reach-j observation, 0 if none.
tt_index <- function(TT) {
  v <- as.vector(t(TT))
  idx <- integer(length(v))
  idx[!is.na(v)] <- seq_len(sum(!is.na(v)))
  list(Nobs = sum(!is.na(v)), Obs = v[!is.na(v)],
       ind = matrix(idx, nrow(TT), byrow = TRUE))
}

last_cap <- function(CH) max.col(CH * col(CH), ties.method = "last")

make_stan_data <- function(design, simM, simT) {
  m <- design$main
  t <- design$trib
  ttM <- tt_index(simM$TT)
  ttT <- tt_index(simT$TT)
  list(
    Nind = m$N, Nreaches = 4L, Ndetlocs = 5L, Nyrs = design$Nyrs, Nrg = design$Nrg,
    UseSizeEffect = 1L, Ntribs = 2L, NindT = t$N, NreachesT = 2L, NrgT = design$NrgT,
    RmultT = t$Rmult, Rmult = m$Rmult, CH = simM$CH,
    yrind = m$yr, rgind = m$rg, rch_covind = 1:4,
    firstCap = rep(1L, m$N), lastCap = last_cap(simM$CH),
    trib_ind = t$trib, firstCapT = rep(1L, t$N), lastCapT = last_cap(simT$CH),
    yrindT = t$yr, rgindT = t$rg, CHT = simT$CH,
    Sz = m$Sz, SzT = t$Sz, CovX = m$CovX, CovXT = t$CovX,
    Nobs = ttM$Nobs, ObsTT = ttM$Obs, TTind = ttM$ind,
    NobsT = ttT$Nobs, ObsTTT = ttT$Obs, TTindT = ttT$ind,
    ReachKM_ind = m$KM, ReachKMT_ind = t$KM
  )
}

truth_vector <- function(tr, pars = tracked_pars) {
  map(pars, \(p) {
    v <- tr[[p]]
    set_names(v, if (length(v) == 1) p else paste0(p, "[", seq_along(v), "]"))
  }) |>
    unlist()
}

# ---------------------------------------------------------------------------
# One replicate: thin, simulate, fit, summarise
# ---------------------------------------------------------------------------
inits_for <- function(model_name, design, n_chains) {
  base <- list(muPb = rep(0, 4), sdPb = rep(1, 4), S_bReach = c(1.5, 0.5, 0.5, 0.5))
  one <- if (model_name == "production") {
    c(base, list(P_b = matrix(2.2, design$Nyrs, 4)))
  } else {
    base
  }
  rep(list(one), n_chains)
}

run_replicate <- function(rep_id, model_name, stan_model_obj, design_full, cfg, out_dir) {
  out_file <- file.path(out_dir, model_name, sprintf("rep_%03d.rds", rep_id))
  if (file.exists(out_file)) {
    message("skip ", model_name, " rep ", rep_id, " (already done)")
    return(invisible(readRDS(out_file)))
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

  # same simulated data set for both models within a replicate
  set.seed(cfg$base_seed + rep_id)
  design <- thin_design(design_full, cfg$thin_frac, cfg$min_per_group)
  tr <- draw_effects(cfg$truth, design)
  simM <- sim_mainstem(tr, design$main)
  simT <- sim_trib(tr, design$trib, design$main$Rmult[1, 4])
  sdata <- make_stan_data(design, simM, simT)

  t0 <- Sys.time()
  fitted <- if (cfg$backend == "cmdstanr") {
    fit_cmdstanr(stan_model_obj, sdata, cfg, model_name, design, rep_id)
  } else {
    fit_rstan(stan_model_obj, sdata, cfg, model_name, design, rep_id)
  }
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  truth <- truth_vector(tr)
  est <- fitted$summary |>
    mutate(model = model_name, rep = rep_id,
           truth = unname(truth[parameter]), .before = 1) |>
    mutate(covered = truth >= q2.5 & truth <= q97.5, error = mean - truth,
           ci_width = q97.5 - q2.5)

  sampler <- tibble(
    model = model_name, rep = rep_id, minutes = minutes, backend = cfg$backend,
    n_main = design$main$N, n_trib = design$trib$N,
    divergences = fitted$divergences,
    treedepth_hits = fitted$treedepth_hits,
    post_warmup_draws = fitted$draws
  )

  res <- list(estimates = est, sampler = sampler, cfg = cfg[setdiff(names(cfg), "truth")])
  saveRDS(res, out_file)
  message(sprintf("%s rep %d: %.1f min, %d divergences, max Rhat %.3f",
                  model_name, rep_id, minutes, sampler$divergences, max(est$Rhat, na.rm = TRUE)))
  invisible(res)
}


# ---------------------------------------------------------------------------
# Backends. Both return the same summary columns:
# parameter, mean, q2.5, q50, q97.5, n_eff (bulk ESS for cmdstanr), Rhat
# ---------------------------------------------------------------------------
fit_rstan <- function(model, sdata, cfg, model_name, design, rep_id) {
  fit <- rstan::sampling(
    model, data = sdata,
    chains = cfg$chains, cores = cfg$cores, iter = cfg$iter, warmup = cfg$iter / 2,
    init = inits_for(model_name, design, cfg$chains),
    pars = tracked_pars, include = TRUE,
    seed = cfg$base_seed + rep_id, refresh = max(1, cfg$iter / 5)
  )
  s <- rstan::summary(fit, pars = tracked_pars, probs = c(0.025, 0.5, 0.975))$summary
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  list(
    summary = tibble(parameter = rownames(s), mean = s[, "mean"], q2.5 = s[, "2.5%"],
                     q50 = s[, "50%"], q97.5 = s[, "97.5%"], n_eff = s[, "n_eff"], Rhat = s[, "Rhat"]),
    divergences = sum(map_dbl(sp, \(x) sum(x[, "divergent__"]))),
    treedepth_hits = sum(map_dbl(sp, \(x) sum(x[, "treedepth__"] >= 10))),
    draws = sum(map_int(sp, nrow))
  )
}

fit_cmdstanr <- function(model, sdata, cfg, model_name, design, rep_id) {
  tmp <- file.path(tempdir(), paste0("simfit_", model_name, "_", rep_id))
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE) # CSVs are big; keep only the summary
  fit <- model$sample(
    data = sdata, seed = cfg$base_seed + rep_id,
    chains = cfg$chains, parallel_chains = cfg$cores,
    iter_warmup = cfg$iter / 2, iter_sampling = cfg$iter / 2,
    init = inits_for(model_name, design, cfg$chains),
    refresh = max(1, cfg$iter / 5), output_dir = tmp, show_messages = FALSE
  )
  s <- fit$summary(variables = tracked_pars, "mean",
                   ~ quantile(.x, probs = c(0.025, 0.5, 0.975)), "rhat", "ess_bulk")
  d <- fit$diagnostic_summary(quiet = TRUE)
  list(
    summary = tibble(parameter = s$variable, mean = s$mean, q2.5 = s$`2.5%`, q50 = s$`50%`,
                     q97.5 = s$`97.5%`, n_eff = s$ess_bulk, Rhat = s$rhat),
    divergences = sum(d$num_divergent),
    treedepth_hits = sum(d$num_max_treedepth),
    draws = cfg$chains * cfg$iter / 2
  )
}

compile_models <- function(model_paths, backend) {
  if (backend == "cmdstanr") {
    map(model_paths, \(f) cmdstanr::cmdstan_model(f))
  } else {
    map(model_paths, \(f) rstan::stan_model(file = f))
  }
}

# ---------------------------------------------------------------------------
# Read everything written so far
# ---------------------------------------------------------------------------
read_sim_results <- function(out_dir) {
  files <- list.files(out_dir, pattern = "^rep_\\d+\\.rds$", recursive = TRUE, full.names = TRUE)
  res <- map(files, readRDS)
  list(estimates = map(res, "estimates") |> bind_rows(),
       sampler = map(res, "sampler") |> bind_rows())
}
