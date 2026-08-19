# Run a combined Sacramento mainstem CJS + CWT survival and travel-time model.
# CJS uses Sacramento River fish only with the original multiple mainstem
# reaches. Butte and Feather tributary CJS likelihoods are excluded. CWT release
# groups are modeled to Knights Landing. Fork-length and flow effects are
# partially shared, with a shared mean plus CJS- and CWT-specific deviations.

rm(list = ls())

library(tidyverse)
library(here)
library(rstan)

source("scripts/02_GetData.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Settings ---------------------------------------------------------------------------
use_size_effect <- 1
nchains <- 3
niter <- 1000
nwarmup <- 500
adapt_delta <- 0.995
max_treedepth <- 15
seed <- 1234

release_distance_km <- c(
  Battle = 300.443064082,
  RBDD = 240.3631
)

inv_logit <- function(x) exp(x) / (1 + exp(x))

default_pcap_prior <- list(
  lt_mu_pCap = qlogis(0.5),
  lt_sd_pCap = 1.5,
  pcap_draws = tibble()
)

extract_knights_pcap_estimate <- function(pcap_file) {
  pcap_kl <- readRDS(pcap_file)

  as.data.frame(pcap_kl, pars = "logit_pCap") %>%
    pivot_longer(
      cols = matches("^logit_pCap\\["),
      names_to = "trial",
      values_to = "logit_estimate"
    ) %>%
    mutate(
      trial = parse_number(trial),
      estimate = inv_logit(logit_estimate)
    )
}

extract_knights_pcap_prior <- function(pcap_file) {
  pcap_draws <- extract_knights_pcap_estimate(pcap_file)
  logit_estimates <- pcap_draws$logit_estimate[is.finite(pcap_draws$logit_estimate)]

  if (length(logit_estimates) <= 1) {
    return(default_pcap_prior)
  }

  pcap_prior <- list(
    lt_mu_pCap = mean(logit_estimates, na.rm = TRUE),
    lt_sd_pCap = sd(logit_estimates, na.rm = TRUE),
    pcap_draws = pcap_draws
  )

  if (!is.finite(pcap_prior$lt_mu_pCap) ||
      !is.finite(pcap_prior$lt_sd_pCap) ||
      pcap_prior$lt_sd_pCap <= 0) {
    return(default_pcap_prior)
  }

  pcap_prior
}

first_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  x[1]
}

check_stan_convergence <- function(fit, pars, rhat_threshold = 1.01, min_ess = 400,
                                   max_treedepth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))
  max_treedepth_hits <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "treedepth__"] >= max_treedepth),
    numeric(1)
  ))
  bfmi <- vapply(
    sampler_params,
    function(x) {
      energy <- x[, "energy__"]
      mean(diff(energy)^2, na.rm = TRUE) / var(energy, na.rm = TRUE)
    },
    numeric(1)
  )
  low_bfmi_chains <- sum(bfmi < 0.3, na.rm = TRUE)

  parameter_summary <- summary(fit, pars = pars)$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    transmute(
      parameter,
      n_eff,
      Rhat,
      low_ess = is.na(n_eff) | n_eff < min_ess,
      high_rhat = is.na(Rhat) | Rhat > rhat_threshold
    )

  list(
    sampler = tibble(
      diagnostic = c(
        "divergent_transitions",
        "max_treedepth_hits",
        "low_bfmi_chains",
        paste0("bfmi_chain_", seq_along(bfmi))
      ),
      value = c(divergences, max_treedepth_hits, low_bfmi_chains, bfmi),
      failed = c(divergences > 0, max_treedepth_hits > 0, low_bfmi_chains > 0, bfmi < 0.3)
    ),
    parameters = parameter_summary,
    converged = divergences == 0 &&
      max_treedepth_hits == 0 &&
      low_bfmi_chains == 0 &&
      !any(parameter_summary$high_rhat, na.rm = TRUE) &&
      !any(parameter_summary$low_ess, na.rm = TRUE)
  )
}

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

# Sacramento mainstem CJS covariates -------------------------------------------------
CovX <- data.frame(cbind(MaxflowSac.z, MaxflowSac.z, MaxflowSac.z, MaxflowDelta.z))

# CWT data ---------------------------------------------------------------------------
cwt_groups <- drerelrec %>%
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
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ReachKM_cwt = unname(release_distance_km[relloc_area]),
    Rmult_cwt = ReachKM_cwt / 100
  ) %>%
  filter(
    is.finite(Nrel), Nrel > 0,
    is.finite(avg_length),
    is.finite(Maxflowsac),
    is.finite(ReachKM_cwt)
  ) %>%
  arrange(year, month, release_group_id)

cwt_tt_obs <- drerelrec %>%
  filter(
    release_group_id %in% cwt_groups$release_group_id,
    release_group_id != 873,
    is.finite(cwt_tt), cwt_tt > 0
  ) %>%
  distinct(release_group_id, tag_code, date, forklength, cwt_tt) %>%
  mutate(TTind = match(release_group_id, cwt_groups$release_group_id)) %>%
  filter(!is.na(TTind)) %>%
  arrange(TTind, date, tag_code)

# Standardize fork length on the combined CJS + CWT size range so the shared
# coefficient has the same interpretation in both likelihoods.
musz <- mean(c(FL, cwt_groups$avg_length), na.rm = TRUE)
sdsz <- sd(c(FL, cwt_groups$avg_length), na.rm = TRUE)
Sz <- (FL - musz) / sdsz

shared_flow_mu <- mean(MaxflowSac, na.rm = TRUE)
shared_flow_sd <- sd(MaxflowSac, na.rm = TRUE)

cwt_groups <- cwt_groups %>%
  mutate(
    Sz_cwt = (avg_length - musz) / sdsz,
    CovX_cwt = (Maxflowsac - shared_flow_mu) / shared_flow_sd
  )

pcap_prior_file <- here("results", "pCap_one_site_skew_re_knights landing.rds")
pcap_prior <- if (file.exists(pcap_prior_file)) {
  extract_knights_pcap_prior(pcap_prior_file)
} else {
  default_pcap_prior
}

pred_fork_length <- seq(
  min(c(FL, cwt_groups$avg_length), na.rm = TRUE),
  max(c(FL, cwt_groups$avg_length), na.rm = TRUE),
  length.out = 100
)

pred_flow <- seq(
  min(c(MaxflowSac, cwt_groups$Maxflowsac), na.rm = TRUE),
  max(c(MaxflowSac, cwt_groups$Maxflowsac), na.rm = TRUE),
  length.out = 100
)

pred_reach_km_combined <- mean(c(rowSums(as.matrix(ReachKM_ind), na.rm = TRUE),
                                cwt_groups$ReachKM_cwt), na.rm = TRUE)

combined_data <- list(
  Nind = Nind,
  Nreaches = Nreaches,
  Ndetlocs = Ndetlocs,
  Nyrs = Nyrs,
  Nrg = Nrg,
  UseSizeEffect = use_size_effect,
  Rmult = Rmult,
  CH = CH,
  yrind = yrind,
  rgind = rgind,
  rch_covind = rch_covind,
  firstCap = firstCap,
  lastCap = lastCap,
  Sz = Sz,
  CovX = CovX,
  Nobs = Nobs,
  ObsTT = ObsTT,
  TTind = TTind,
  ReachKM_ind = ReachKM_ind,
  Ncwtgrp = nrow(cwt_groups),
  Nobs_cwt = nrow(cwt_tt_obs),
  Rmult_cwt = cwt_groups$Rmult_cwt,
  Sz_cwt = cwt_groups$Sz_cwt,
  CovX_cwt = cwt_groups$CovX_cwt,
  ReachKM_cwt = cwt_groups$ReachKM_cwt,
  ObsTT_cwt = cwt_tt_obs$cwt_tt,
  TTind_cwt = cwt_tt_obs$TTind,
  Nrel_cwt = as.integer(round(cwt_groups$Nrel)),
  cwt_recaptures = as.integer(cwt_groups$cwt_recaptures),
  lt_mu_pCap_cwt = pcap_prior$lt_mu_pCap,
  lt_sd_pCap_cwt = pcap_prior$lt_sd_pCap,
  Npred_size = length(pred_fork_length),
  pred_size_z = (pred_fork_length - musz) / sdsz,
  pred_cov_mean = 0,
  Npred_flow = length(pred_flow),
  pred_cov_combined = (pred_flow - shared_flow_mu) / shared_flow_sd,
  pred_size_mean_z = 0,
  pred_reach_km_combined = pred_reach_km_combined
)

saveRDS(combined_data, here("results", "combined_cjs_cwt_stan_data.rds"))
write.csv(cwt_groups, here("results", "tables", "combined_cwt_release_groups.csv"), row.names = FALSE)
write.csv(cwt_tt_obs, here("results", "tables", "combined_cwt_travel_time_obs.csv"), row.names = FALSE)

inits1 <- list(
  P_b = matrix(data = 2.2, nrow = Nyrs, ncol = Nreaches),
  muPb = rep(0, Nreaches),
  sdPb = rep(1.0, Nreaches),
  S_bCJS = rep(0.5, Nreaches),
  T_bCJS = rep(0, Nreaches),
  S_bCWT = 0.5,
  T_bCWT = 0
)
inits <- replicate(nchains, inits1, simplify = FALSE)

fit_combined <- stan(
  file = here("scripts", "Combined_CJS_CWT_FL_SharedMainstem.stan"),
  data = combined_data,
  init = inits,
  chains = nchains,
  iter = niter,
  warmup = nwarmup,
  seed = seed,
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
)

save(
  fit_combined, combined_data, cwt_groups, cwt_tt_obs, pred_fork_length, pred_flow,
  file = here("results", "fit_combined_CJS_CWT_FL.Rdata")
)

diagnostic_pars <- c(
  "S_bCJS", "S_bCWT",
  "S_bCov_shared", "S_bCov_cjs", "S_bCov_cwt",
  "S_bSz_shared", "S_bSz_cjs", "S_bSz_cwt",
  "T_bCJS", "T_bCWT",
  "TT_bCov_cjs", "TT_bCov_cwt",
  "T_bSz_shared", "T_bSz_cjs", "T_bSz_cwt",
  "sigma_S_bCov_model", "sigma_S_bSz_model", "sigma_T_bSz_model",
  "RE_sd", "sd_cwtSre",
  "TTRE_sd", "sd_cwtTTre", "Pro_sd", "Pro_sd_cwt"
)

diagnostic_summary <- summary(fit_combined, pars = diagnostic_pars)$summary %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  transmute(
    parameter,
    mean,
    sd,
    q2.5 = `2.5%`,
    q50 = `50%`,
    q97.5 = `97.5%`,
    n_eff,
    Rhat,
    low_ess = n_eff < 400,
    high_rhat = is.na(Rhat) | Rhat > 1.01
  ) %>%
  arrange(desc(high_rhat), desc(low_ess), Rhat, n_eff)

write.csv(
  diagnostic_summary,
  here("results", "tables", "combined_cjs_cwt_convergence_diagnostics.csv"),
  row.names = FALSE
)

write.csv(
  diagnostic_summary %>% filter(low_ess | high_rhat),
  here("results", "tables", "combined_cjs_cwt_convergence_diagnostics_flagged.csv"),
  row.names = FALSE
)

sampler_params <- get_sampler_params(fit_combined, inc_warmup = FALSE)
n_divergent <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
n_transitions <- sum(vapply(sampler_params, nrow, integer(1)))
n_treedepth <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)))

sampler_diagnostics <- tibble(
  diagnostic = c(
    "post_warmup_transitions",
    "divergent_transitions",
    "divergent_transition_rate",
    "max_treedepth_hits",
    "max_treedepth_hit_rate"
  ),
  value = c(
    n_transitions,
    n_divergent,
    n_divergent / n_transitions,
    n_treedepth,
    n_treedepth / n_transitions
  )
)

sampler_parameter_summary <- bind_rows(
  lapply(seq_along(sampler_params), function(chain_id) {
    as.data.frame(sampler_params[[chain_id]]) %>%
      summarise(
        across(
          everything(),
          list(mean = mean, min = min, max = max),
          .names = "{.col}_{.fn}"
        )
      ) %>%
      mutate(chain = chain_id, .before = 1)
  })
)

write.csv(
  sampler_diagnostics,
  here("results", "tables", "combined_cjs_cwt_sampler_diagnostics.csv"),
  row.names = FALSE
)

write.csv(
  sampler_parameter_summary,
  here("results", "tables", "combined_cjs_cwt_sampler_parameter_summary.csv"),
  row.names = FALSE
)

convergence_check <- check_stan_convergence(
  fit_combined,
  diagnostic_pars,
  max_treedepth = max_treedepth
)

write.csv(
  convergence_check$sampler,
  here("results", "tables", "combined_cjs_cwt_convergence_check_sampler.csv"),
  row.names = FALSE
)

write.csv(
  convergence_check$parameters,
  here("results", "tables", "combined_cjs_cwt_convergence_check_parameters.csv"),
  row.names = FALSE
)

if (!isTRUE(convergence_check$converged)) {
  warning(
    "Combined CJS-CWT model convergence diagnostics were flagged. ",
    "Check results/tables/combined_cjs_cwt_convergence_diagnostics_flagged.csv, ",
    "combined_cjs_cwt_convergence_check_sampler.csv, and ",
    "results/figures/combined_cjs_cwt_trace_diagnostics.png."
  )
}

trace_pars <- c(
  "S_bCov_shared", "S_bCov_cjs", "S_bCov_cwt",
  "S_bSz_shared", "S_bSz_cjs", "S_bSz_cwt",
  "TT_bCov_cjs", "TT_bCov_cwt",
  "T_bSz_shared", "T_bSz_cjs", "T_bSz_cwt",
  "sigma_S_bCov_model",
  "sigma_S_bSz_model", "sigma_T_bSz_model",
  "S_bCWT", "T_bCWT", "sd_cwtSre", "sd_cwtTTre",
  "Pro_sd", "Pro_sd_cwt"
)

png(
  filename = here("results", "figures", "combined_cjs_cwt_trace_diagnostics.png"),
  width = 1400,
  height = 1000,
  res = 130
)
print(traceplot(fit_combined, pars = trace_pars))
dev.off()

size_prediction_pars <- c(
  "pred_combined_survival_by_size",
  "pred_combined_travel_time_by_size"
)

flow_prediction_pars <- c(
  "pred_combined_survival_by_flow",
  "pred_combined_travel_time_by_flow"
)

size_prediction_summary <- bind_rows(lapply(size_prediction_pars, function(par_name) {
  summary(fit_combined, pars = par_name)$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    mutate(
      prediction = par_name,
      fork_length = pred_fork_length,
      prediction_model = "Hierarchical shared mean relationship"
    ) %>%
    select(prediction, fork_length, prediction_model, mean, sd, `2.5%`, `50%`, `97.5%`, n_eff, Rhat)
}))

write.csv(
  size_prediction_summary,
  here("results", "tables", "combined_cjs_cwt_predictions_by_fork_length.csv"),
  row.names = FALSE
)

flow_prediction_summary <- bind_rows(lapply(flow_prediction_pars, function(par_name) {
  summary(fit_combined, pars = par_name)$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    mutate(
      prediction = par_name,
      flow = pred_flow,
      prediction_model = "Hierarchical shared mean relationship"
    ) %>%
    select(prediction, flow, prediction_model, mean, sd, `2.5%`, `50%`, `97.5%`, n_eff, Rhat)
}))

write.csv(
  flow_prediction_summary,
  here("results", "tables", "combined_cjs_cwt_predictions_by_flow.csv"),
  row.names = FALSE
)

size_prediction_plot <- size_prediction_summary %>%
  mutate(
    outcome = case_when(
      grepl("survival", prediction) ~ "Survival",
      TRUE ~ "Travel time"
    )
  ) %>%
  ggplot(aes(x = fork_length, y = mean, ymin = `2.5%`, ymax = `97.5%`,
             colour = prediction_model, fill = prediction_model)) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(
    x = "Fork length (mm)",
    y = "Posterior prediction",
    colour = NULL,
    fill = NULL
  )

ggsave(
  here("results", "figures", "combined_cjs_cwt_predictions_by_fork_length.png"),
  plot = size_prediction_plot,
  width = 9,
  height = 8,
  dpi = 350
)

flow_prediction_plot <- flow_prediction_summary %>%
  mutate(
    outcome = case_when(
      grepl("survival", prediction) ~ "Survival",
      TRUE ~ "Travel time"
    )
  ) %>%
  ggplot(aes(x = flow, y = mean, ymin = `2.5%`, ymax = `97.5%`,
             colour = prediction_model, fill = prediction_model)) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(
    x = "Flow",
    y = "Posterior prediction",
    colour = NULL,
    fill = NULL
  )

ggsave(
  here("results", "figures", "combined_cjs_cwt_predictions_by_flow.png"),
  plot = flow_prediction_plot,
  width = 9,
  height = 8,
  dpi = 350
)
