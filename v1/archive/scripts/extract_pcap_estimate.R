library(tidyverse)
library(rstan)

inv_logit <- function(x) exp(x) / (1 + exp(x))

extract_knights_pcap_estimate <- function(
  pcap_file = "results/pCap_one_site_skew_re_knights landing.rds"
) {
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

summarise_knights_pcap <- function(
  pcap_file = "results/pCap_one_site_skew_re_knights landing.rds"
) {
  extract_knights_pcap_estimate(pcap_file) %>%
    group_by(trial) %>%
    summarise(
      lcl_10 = quantile(estimate, 0.1, na.rm = TRUE),
      median = median(estimate, na.rm = TRUE),
      mean = mean(estimate, na.rm = TRUE),
      ucl_90 = quantile(estimate, 0.9, na.rm = TRUE),
      .groups = "drop"
    )
}

extract_knights_pcap_prior <- function(
  pcap_file = "results/pCap_one_site_skew_re_knights landing.rds"
) {
  pcap_draws <- extract_knights_pcap_estimate(pcap_file)
  logit_estimates <- pcap_draws$logit_estimate[is.finite(pcap_draws$logit_estimate)]
  
  if (length(logit_estimates) <= 1) {
    return(list(
      lt_mu_pCap = qlogis(0.5),
      lt_sd_pCap = 1.5,
      pcap_by_trial = summarise_knights_pcap(pcap_file)
    ))
  }
  
  list(
    lt_mu_pCap = mean(logit_estimates, na.rm = TRUE),
    lt_sd_pCap = sd(logit_estimates, na.rm = TRUE),
    pcap_by_trial = summarise_knights_pcap(pcap_file)
  )
}

if (sys.nframe() == 0) {
  pcap_prior <- extract_knights_pcap_prior()
  pcap_by_trial <- pcap_prior$pcap_by_trial
  
  print(pcap_by_trial)
  cat("lt_mu_pCap:", pcap_prior$lt_mu_pCap, "\n")
  cat("lt_sd_pCap:", pcap_prior$lt_sd_pCap, "\n")
}
