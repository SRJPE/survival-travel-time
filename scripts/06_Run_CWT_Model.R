# Run the CWT survival and travel-time model.
# prepares the CWT.stan data list, fits the model, and saves posterior summaries.

rm(list = ls())

library(tidyverse)
library(here)
library(lubridate)
library(readxl)
library(rstan)

source("scripts/02_GetData.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Settings ---------------------------------------------------------------------------
use_size_effect <- 1
nchains <- 4
niter <- 1000
nwarmup <- 500
adapt_delta <- 0.99
seed <- 123

# Distance from each CWT release area to Knights Landing, in river km.
# These values match the RST-to-Knights distances used in GetData_flora.R:
# Battle = Battle Creek, RBDD = Red Bluff Diversion Dam area.
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

# Prior for Knights Landing capture probability on the logit scale.
# Use the posterior logit pCap draws from the Knights Landing pCap model.
pcap_prior_file <- here("results", "pCap_one_site_skew_re_knights landing.rds")
pcap_prior <- if (file.exists(pcap_prior_file)) {
  extract_knights_pcap_prior(pcap_prior_file)
} else {
  default_pcap_prior
}

lt_mu_pCap <- pcap_prior$lt_mu_pCap
lt_sd_pCap <- pcap_prior$lt_sd_pCap

pcap_prior_summary <- tibble(
  parameter = c("lt_mu_pCap", "lt_sd_pCap"),
  value = c(lt_mu_pCap, lt_sd_pCap)
)

if (nrow(pcap_prior$pcap_draws) > 0) {
  write.csv(
    pcap_prior$pcap_draws,
    here("results", "tables", "knights_pcap_prior_draws.csv"),
    row.names = FALSE
  )
}

write.csv(
  pcap_prior_summary,
  here("results", "tables", "knights_pcap_prior_summary.csv"),
  row.names = FALSE
)

# Prepare CWT data -------------------------------------------------------------------
# Collapse to one row per release group for the survival likelihood, and keep only
# the Battle and RBDD release areas
first_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  x[1]
}

cwt_groups <- drerelrec %>%
  filter(!is.na(relloc_area), relloc_area %in% names(release_distance_km),
         release_group_id != 873) %>% # Avg_length is inconsistent so remove for now
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


# One travel-time observation per unique recaptured fish/date/length record.
cwt_tt_obs <- drerelrec %>%
  filter(
    release_group_id %in% cwt_groups$release_group_id,
    release_group_id != 873, # Avg_length is inconsistent so remove for now
    is.finite(cwt_tt), cwt_tt > 0
  ) %>%
  distinct(release_group_id, tag_code, date, forklength, cwt_tt) %>%
  mutate(TTind = match(release_group_id, cwt_groups$release_group_id)) %>%
  filter(!is.na(TTind)) %>%
  arrange(TTind, date, tag_code)

# Standardize covariates using the CWT release groups in this model.
size_mu <- mean(cwt_groups$avg_length, na.rm = TRUE)
size_sd <- sd(cwt_groups$avg_length, na.rm = TRUE)
flow_mu <- mean(cwt_groups$Maxflowsac, na.rm = TRUE)
flow_sd <- sd(cwt_groups$Maxflowsac, na.rm = TRUE)

cwt_groups <- cwt_groups %>%
  mutate(
    Sz_cwt = (avg_length - size_mu) / size_sd,
    CovX_cwt = (Maxflowsac - flow_mu) / flow_sd
  )

cwt_data <- list(
  Ncwtgrp = nrow(cwt_groups),
  Nyrs = n_distinct(cwt_groups$year),
  Nobs = nrow(cwt_tt_obs),
  UseSizeEffect = use_size_effect,
  Rmult_cwt = cwt_groups$Rmult_cwt,
  Sz_cwt = cwt_groups$Sz_cwt,
  CovX_cwt = cwt_groups$CovX_cwt,
  ReachKM_cwt = cwt_groups$ReachKM_cwt,
  ObsTT = cwt_tt_obs$cwt_tt,
  TTind = cwt_tt_obs$TTind,
  Nrel = as.integer(round(cwt_groups$Nrel)),
  cwt_recaptures = as.integer(cwt_groups$cwt_recaptures),
  lt_mu_pCap = lt_mu_pCap,
  lt_sd_pCap = lt_sd_pCap
)

saveRDS(cwt_data, here("results", "cwt_stan_data.rds"))
write.csv(cwt_groups, here("results/tables", "cwt_release_groups.csv"), row.names = FALSE)
write.csv(cwt_tt_obs, here("results/tables", "cwt_travel_time_obs.csv"), row.names = FALSE)

# Fit model --------------------------------------------------------------------------
fit_cwt <- stan(
  file = here("scripts", "CWT.stan"),
  data = cwt_data,
  chains = nchains,
  iter = niter,
  warmup = nwarmup,
  seed = seed,
  control = list(adapt_delta = adapt_delta)
)

save(fit_cwt, cwt_data, cwt_groups, cwt_tt_obs,
     file = here("results", "fit_CWT.Rdata"))

# Model diagnostics ------------------------------------------------------------------
diagnostic_pars <- c(
  "S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz",
  "sd_cwtSre", "sd_cwtTTre", "logit_pCap_Sim", "Pro_sd"
)

diagnostic_summary <- summary(fit_cwt, pars = diagnostic_pars)$summary %>%
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
  arrange(desc(high_rhat), Rhat, n_eff)

sampler_params <- get_sampler_params(fit_cwt, inc_warmup = FALSE)
n_divergent <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
n_transitions <- sum(vapply(sampler_params, nrow, integer(1)))
n_treedepth <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 10), numeric(1)))

sampler_diagnostics <- tibble(
  diagnostic = c("post_warmup_transitions", "divergent_transitions",
                 "divergent_transition_rate", "max_treedepth_hits",
                 "max_treedepth_hit_rate"),
  value = c(n_transitions, n_divergent, n_divergent / n_transitions,
            n_treedepth, n_treedepth / n_transitions)
)

sampler_param_summary <- bind_rows(
  lapply(seq_along(sampler_params), function(chain_id) {
    as.data.frame(sampler_params[[chain_id]]) %>%
      summarise(across(everything(),
                       list(mean = mean, min = min, max = max),
                       .names = "{.col}_{.fn}")) %>%
      mutate(chain = chain_id, .before = 1)
  })
)

write.csv(diagnostic_summary,
          here("results/tables", "cwt_convergence_diagnostics.csv"),
          row.names = FALSE)
write.csv(sampler_diagnostics,
          here("results/tables", "cwt_sampler_diagnostics.csv"),
          row.names = FALSE)
write.csv(sampler_param_summary,
          here("results/tables", "cwt_sampler_parameter_summary.csv"),
          row.names = FALSE)

png(filename = here("results/figures", "cwt_trace_diagnostics.png"),
    width = 1200, height = 900, res = 130)
print(traceplot(fit_cwt, pars = diagnostic_pars))
dev.off()

png(filename = here("results/figures", "cwt_pairs_diagnostics.png"),
    width = 1200, height = 900, res = 130)
pairs(fit_cwt, pars = c("S_bReach", "T_bReach", "S_bCov", "TT_bCov",
                        "S_bSz", "T_bSz", "sd_cwtSre", "sd_cwtTTre",
                        "Pro_sd"))
dev.off()

# Posterior summaries ----------------------------------------------------------------
pars_to_save <- c(
  "S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz",
  "sd_cwtSre", "sd_cwtTTre", "logit_pCap_Sim", "Pro_sd",
  "surv_cwt", "pKL", "pTT", "log_lik"
)

cwt_summary <- summary(fit_cwt, pars = pars_to_save)$summary
write.csv(cwt_summary, here("results/tables", "cwt_model_summary.csv"))

surv_summary <- summary(fit_cwt, pars = "surv_cwt")$summary %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(group_index = row_number()) %>%
  bind_cols(cwt_groups %>% select(release_group_id, year, month, release_location_name, relloc_area)) %>%
  select(group_index, release_group_id, year, month, release_location_name, relloc_area,
         mean, sd, `2.5%`, `50%`, `97.5%`, n_eff, Rhat)

tt_summary <- summary(fit_cwt, pars = "pTT")$summary %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  mutate(group_index = row_number()) %>%
  bind_cols(cwt_groups %>% select(release_group_id, year, month, release_location_name, relloc_area)) %>%
  select(group_index, release_group_id, year, month, release_location_name, relloc_area,
         mean, sd, `2.5%`, `50%`, `97.5%`, n_eff, Rhat)

write.csv(surv_summary, here("results/tables", "cwt_survival_summary.csv"), row.names = FALSE)
write.csv(tt_summary, here("results/tables", "cwt_travel_time_summary.csv"), row.names = FALSE)

