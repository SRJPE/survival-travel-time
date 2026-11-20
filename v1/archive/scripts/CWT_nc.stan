// CWT survival and travel time model: NON-CENTERED version of CWT.stan
// (2026-09). Identical likelihood and priors to CWT.stan; only the release-group
// random effects are reparameterized:
//   S_cwtre  = sd_cwtSre  * z_S,   z_S  ~ std_normal()
//   TT_cwtre = sd_cwtTTre * z_TT,  z_TT ~ std_normal()
// With a centered parameterization, sd_cwtSre, S_bReach and logit_pCap_Sim
// sat on a correlated ridge (Rhat ~1.07, ESS ~55 after filtering the
// no-reading years). Non-centering usually removes the funnel between a
// group-level SD and its effects.
data {
  int Ncwtgrp;
  int Nyrs;
  int Nobs;
  int UseSizeEffect;
  array[Ncwtgrp] real Rmult_cwt;
  array[Ncwtgrp] real Sz_cwt;
  array[Ncwtgrp] real CovX_cwt;
  array[Ncwtgrp] real ReachKM_cwt;

  vector[Nobs] ObsTT;
  array[Nobs] int TTind;

  array[Ncwtgrp] int Nrel;
  array[Ncwtgrp] int cwt_recaptures;

  real<lower=-10, upper=10> lt_mu_pCap;
  real<lower=-10, upper=10> lt_sd_pCap;
}

parameters {
  real<lower=-5, upper=5> S_bReach;
  real<lower=-5, upper=5> T_bReach;
  real<lower=-10, upper=10> S_bCov;
  real<lower=-10, upper=10> TT_bCov;
  real<lower=-10, upper=10> S_bSz;
  real<lower=-10, upper=10> T_bSz;
  vector[Ncwtgrp] z_S;
  real<lower=0.001> sd_cwtSre;
  vector[Ncwtgrp] z_TT;
  real<lower=0.001> sd_cwtTTre;
  real logit_pCap_Sim;
  real<lower=0.001, upper=5> Pro_sd;
}

transformed parameters {
  vector[Ncwtgrp] S_cwtre = sd_cwtSre * z_S;
  vector[Ncwtgrp] TT_cwtre = sd_cwtTTre * z_TT;
  array[Ncwtgrp] real surv_cwt;
  array[Ncwtgrp] real pKL;
  array[Ncwtgrp] real pTT;

  for (i in 1:Ncwtgrp) {
    surv_cwt[i] = inv_logit(S_bReach + UseSizeEffect * S_bSz * Sz_cwt[i] +
                            S_bCov * CovX_cwt[i] + S_cwtre[i])^Rmult_cwt[i];
    pKL[i] = surv_cwt[i] * inv_logit(logit_pCap_Sim);
    pTT[i] = exp(T_bReach + UseSizeEffect * T_bSz * Sz_cwt[i] +
                 TT_bCov * CovX_cwt[i] + TT_cwtre[i]) * ReachKM_cwt[i] / 100;
  }
}

model {
  S_bReach ~ normal(0, 1.5);
  T_bReach ~ normal(0, 1.5);
  S_bCov ~ normal(0, 1);
  TT_bCov ~ normal(0, 1);
  S_bSz ~ normal(0, 1);
  T_bSz ~ normal(0, 1);
  sd_cwtSre ~ normal(0, 1);
  sd_cwtTTre ~ normal(0, 1);
  Pro_sd ~ normal(0, 1);

  z_S ~ std_normal();
  z_TT ~ std_normal();

  logit_pCap_Sim ~ normal(lt_mu_pCap, lt_sd_pCap); // from BT-SPAS-X, read in as data

  cwt_recaptures ~ binomial(Nrel, pKL);

  for (i in 1:Nobs) {
    ObsTT[i] ~ lognormal(log(pTT[TTind[i]]), Pro_sd);
  }
}

generated quantities {
  vector[Ncwtgrp] log_lik;
  for (i in 1:Ncwtgrp) {
    log_lik[i] = binomial_lpmf(cwt_recaptures[i] | Nrel[i], pKL[i]);
  }
}
