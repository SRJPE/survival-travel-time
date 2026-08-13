// Combined Sacramento mainstem CJS + CWT survival and travel-time model.
// CJS uses Sacramento River mainstem fish only with multiple mainstem reaches.
// Butte and Feather tributary CJS likelihoods are excluded. CWT release groups
// are modeled to Knights Landing. Fork-length and flow effects are partially
// shared: each effect has a shared mean plus CJS- and CWT-specific deviations.
data {
  int Nind;
  int Nreaches;
  int Ndetlocs;
  int Nyrs;
  int Nrg;
  int UseSizeEffect;
  array[Nind, Nreaches] real Rmult;
  array[Nind, Ndetlocs] int CH;
  array[Nind] int yrind;
  array[Nind] int rgind;
  array[Nreaches] int rch_covind;
  array[Nind] int firstCap;
  array[Nind] int lastCap;
  array[Nind] real Sz;
  array[Nind, Nreaches] real CovX;

  int Nobs;
  vector[Nobs] ObsTT;
  array[Nind, Nreaches] int TTind;
  array[Nind, Nreaches] real ReachKM_ind;

  int Ncwtgrp;
  int Nobs_cwt;
  array[Ncwtgrp] real Rmult_cwt;
  array[Ncwtgrp] real Sz_cwt;
  array[Ncwtgrp] real CovX_cwt;
  array[Ncwtgrp] real ReachKM_cwt;
  vector[Nobs_cwt] ObsTT_cwt;
  array[Nobs_cwt] int TTind_cwt;
  array[Ncwtgrp] int Nrel_cwt;
  array[Ncwtgrp] int cwt_recaptures;
  real<lower=-10, upper=10> lt_mu_pCap_cwt;
  real<lower=0> lt_sd_pCap_cwt;

  int Npred_size;
  array[Npred_size] real pred_size_z;
  real pred_cov_mean;
  int Npred_flow;
  array[Npred_flow] real pred_cov_combined;
  real pred_size_mean_z;
  real pred_reach_km_combined;
}

parameters {
  array[Nyrs] vector<lower=-10, upper=10>[Nreaches] P_b;
  vector<lower=-10, upper=10>[Nreaches] muPb;
  vector<lower=0.001>[Nreaches] sdPb;

  vector<lower=-5, upper=5>[Nreaches] S_bCJS;
  vector<lower=-5, upper=5>[Nreaches] T_bCJS;
  real<lower=-5, upper=5> S_bCWT;
  real<lower=-5, upper=5> T_bCWT;

  real<lower=-10, upper=10> S_bCov_shared;
  real<lower=-10, upper=10> TT_bCov_shared;
  real<lower=-10, upper=10> S_bSz_shared;
  real<lower=-10, upper=10> T_bSz_shared;
  vector[2] z_S_bCov_model;
  vector[2] z_TT_bCov_model;
  vector[2] z_S_bSz_model;
  vector[2] z_T_bSz_model;
  real<lower=0.001> sigma_S_bCov_model;
  real<lower=0.001> sigma_TT_bCov_model;
  real<lower=0.001> sigma_S_bSz_model;
  real<lower=0.001> sigma_T_bSz_model;

  array[Nrg, Nreaches] real<lower=-10, upper=10> S_RE;
  vector<lower=0.001>[Nreaches] RE_sd;
  array[Nrg, Nreaches] real<lower=-10, upper=10> TT_RE;
  vector<lower=0.001>[Nreaches] TTRE_sd;

  vector[Ncwtgrp] S_cwtre;
  real<lower=0.001> sd_cwtSre;
  vector[Ncwtgrp] TT_cwtre;
  real<lower=0.001> sd_cwtTTre;
  real logit_pCap_cwt;

  real<lower=0.001, upper=5> Pro_sd;
  real<lower=0.001, upper=5> Pro_sd_cwt;
}

transformed parameters {
  real S_bCov_cjs = S_bCov_shared + sigma_S_bCov_model * z_S_bCov_model[1];
  real S_bCov_cwt = S_bCov_shared + sigma_S_bCov_model * z_S_bCov_model[2];
  real TT_bCov_cjs = TT_bCov_shared + sigma_TT_bCov_model * z_TT_bCov_model[1];
  real TT_bCov_cwt = TT_bCov_shared + sigma_TT_bCov_model * z_TT_bCov_model[2];
  real S_bSz_cjs = S_bSz_shared + sigma_S_bSz_model * z_S_bSz_model[1];
  real S_bSz_cwt = S_bSz_shared + sigma_S_bSz_model * z_S_bSz_model[2];
  real T_bSz_cjs = T_bSz_shared + sigma_T_bSz_model * z_T_bSz_model[1];
  real T_bSz_cwt = T_bSz_shared + sigma_T_bSz_model * z_T_bSz_model[2];

  array[Nind] vector[Nreaches] surv;
  array[Nind] vector[Ndetlocs] Pcap;
  array[Nind] vector[Ndetlocs] chi;
  array[Nind] vector[Nreaches] pTT;
  vector[Nobs] lg_pTT;
  array[Ncwtgrp] real surv_cwt;
  array[Ncwtgrp] real pKL_cwt;
  array[Ncwtgrp] real pTT_cwt;

  for (i in 1:Nind) {
    real TT;
    TT = exp(T_bCJS[1] + UseSizeEffect * T_bSz_cjs * Sz[i] +
             TT_bCov_cjs * CovX[i, 1] + TT_RE[rgind[i], 1]);
    pTT[i, 1] = TT * ReachKM_ind[i, 1] / 100;
    if (TTind[i, 1] > 0) lg_pTT[TTind[i, 1]] = log(pTT[i, 1]);

    for (j in 2:Nreaches) {
      TT = exp(T_bCJS[j] + UseSizeEffect * T_bSz_cjs * Sz[i] +
               TT_bCov_cjs * CovX[i, j] + TT_RE[rgind[i], j]);
      pTT[i, j] = pTT[i, j - 1] + TT * ReachKM_ind[i, j] / 100;
      if (TTind[i, j] > 0) lg_pTT[TTind[i, j]] = log(pTT[i, j]);
    }
  }

  for (i in 1:Nind) {
    for (j in 1:Ndetlocs) {
      if (j < Ndetlocs) {
        surv[i, j] = inv_logit(
          S_bCJS[rch_covind[j]] + S_bCov_cjs * CovX[i, j] +
          UseSizeEffect * S_bSz_cjs * Sz[i] + S_RE[rgind[i], j]
        ) ^ Rmult[i, j];
      }
      if (j > 1) Pcap[i, j] = inv_logit(P_b[yrind[i], j - 1]);
    }
    chi[i, Ndetlocs] = 1.0;
    for (j in 1:Nreaches) {
      int r_curr = Ndetlocs - j;
      int r_next = r_curr + 1;
      chi[i, r_curr] = (1 - surv[i, r_curr]) +
        surv[i, r_curr] * (1 - Pcap[i, r_next]) * chi[i, r_next];
    }
  }

  for (i in 1:Ncwtgrp) {
    surv_cwt[i] = inv_logit(
      S_bCWT + S_bCov_cwt * CovX_cwt[i] +
      UseSizeEffect * S_bSz_cwt * Sz_cwt[i] + S_cwtre[i]
    ) ^ Rmult_cwt[i];
    pKL_cwt[i] = surv_cwt[i] * inv_logit(logit_pCap_cwt);
    pTT_cwt[i] = exp(
      T_bCWT + TT_bCov_cwt * CovX_cwt[i] +
      UseSizeEffect * T_bSz_cwt * Sz_cwt[i] + TT_cwtre[i]
    ) * ReachKM_cwt[i] / 100;
  }
}

model {
  for (j in 1:Nreaches) {
    P_b[1:Nyrs, j] ~ normal(muPb[j], sdPb[j]);
    S_RE[, j] ~ normal(0, RE_sd[j]);
    TT_RE[, j] ~ normal(0, TTRE_sd[j]);
  }
  muPb ~ normal(0, 1.5);
  sdPb ~ normal(0, 1);

  S_bCJS ~ normal(0, 1.5);
  T_bCJS ~ normal(0, 1.5);
  S_bCWT ~ normal(0, 1.5);
  T_bCWT ~ normal(0, 1.5);
  S_bCov_shared ~ normal(0, 1);
  TT_bCov_shared ~ normal(0, 1);
  S_bSz_shared ~ normal(0, 1);
  T_bSz_shared ~ normal(0, 1);
  z_S_bCov_model ~ normal(0, 1);
  z_TT_bCov_model ~ normal(0, 1);
  z_S_bSz_model ~ normal(0, 1);
  z_T_bSz_model ~ normal(0, 1);
  sigma_S_bCov_model ~ normal(0, 0.5);
  sigma_TT_bCov_model ~ normal(0, 0.5);
  sigma_S_bSz_model ~ normal(0, 0.5);
  sigma_T_bSz_model ~ normal(0, 0.5);

  RE_sd ~ normal(0, 1);
  TTRE_sd ~ normal(0, 1);
  S_cwtre ~ normal(0, sd_cwtSre);
  TT_cwtre ~ normal(0, sd_cwtTTre);
  sd_cwtSre ~ normal(0, 1);
  sd_cwtTTre ~ normal(0, 1);
  Pro_sd ~ normal(0, 1);
  Pro_sd_cwt ~ normal(0, 1);

  logit_pCap_cwt ~ normal(lt_mu_pCap_cwt, lt_sd_pCap_cwt);

  ObsTT ~ lognormal(lg_pTT, Pro_sd);
  for (i in 1:Nobs_cwt) {
    ObsTT_cwt[i] ~ lognormal(log(pTT_cwt[TTind_cwt[i]]), Pro_sd_cwt);
  }

  for (i in 1:Nind) {
    for (j in (firstCap[i] + 1):lastCap[i]) {
      1 ~ bernoulli(surv[i, j - 1]);
      CH[i, j] ~ bernoulli(Pcap[i, j]);
    }
    1 ~ bernoulli(chi[i, lastCap[i]]);
  }

  for (i in 1:Ncwtgrp) {
    cwt_recaptures[i] ~ binomial(Nrel_cwt[i], pKL_cwt[i]);
  }
}

generated quantities {
  vector[Nind + Ncwtgrp] log_lik;
  vector[Npred_size] pred_combined_survival_by_size;
  vector[Npred_size] pred_combined_travel_time_by_size;
  vector[Npred_flow] pred_combined_survival_by_flow;
  vector[Npred_flow] pred_combined_travel_time_by_flow;

  for (k in 1:Npred_size) {
    pred_combined_survival_by_size[k] = inv_logit(
      mean(S_bCJS) * 0.5 + 0.5 * S_bCWT +
      S_bCov_shared * pred_cov_mean +
      UseSizeEffect * S_bSz_shared * pred_size_z[k]
    ) ^ (pred_reach_km_combined / 100);
    pred_combined_travel_time_by_size[k] = exp(
      mean(T_bCJS) * 0.5 + 0.5 * T_bCWT +
      TT_bCov_shared * pred_cov_mean +
      UseSizeEffect * T_bSz_shared * pred_size_z[k]
    ) * pred_reach_km_combined / 100;
  }

  for (k in 1:Npred_flow) {
    pred_combined_survival_by_flow[k] = inv_logit(
      mean(S_bCJS) * 0.5 + 0.5 * S_bCWT +
      S_bCov_shared * pred_cov_combined[k] +
      UseSizeEffect * S_bSz_shared * pred_size_mean_z
    ) ^ (pred_reach_km_combined / 100);
    pred_combined_travel_time_by_flow[k] = exp(
      mean(T_bCJS) * 0.5 + 0.5 * T_bCWT +
      TT_bCov_shared * pred_cov_combined[k] +
      UseSizeEffect * T_bSz_shared * pred_size_mean_z
    ) * pred_reach_km_combined / 100;
  }

  for (i in 1:Nind) {
    log_lik[i] = 0;
    for (j in (firstCap[i] + 1):lastCap[i]) {
      log_lik[i] += bernoulli_lpmf(1 | surv[i, j - 1]);
      log_lik[i] += bernoulli_lpmf(CH[i, j] | Pcap[i, j]);
    }
    log_lik[i] += bernoulli_lpmf(1 | chi[i, lastCap[i]]);
  }

  for (i in 1:Ncwtgrp) {
    log_lik[Nind + i] = binomial_lpmf(cwt_recaptures[i] | Nrel_cwt[i], pKL_cwt[i]);
  }
}
