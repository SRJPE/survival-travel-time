// CWT_oneway_gamma.stan (9 Oct 2026): CWT_oneway.stan with ONE change, an
// exposure exponent gamma on time at large:
//   cumulative hazard_g = M_r * (L_r / L_g)^c * exp(delta + h_g) * T_g^gamma
//   M_r(c, gamma) = H_tel / sum_i T_i^gamma (L_r / L_i)^c   (telemetry, same form)
// gamma = 1 is the original per-day model. gamma < 1: slow fish accumulate
// less mortality than days at large imply (rearing, holding). gamma = 0:
// mortality per trip, i.e. per distance for a fixed reach (like sim S4).
// UseGamma = 0 fixes gamma at 1 (reproduces CWT_oneway.stan). In the code the
// exponent is "gam" (gamma is a Stan function name); gamma_raw is the parameter.
// CWT_oneway.stan (bp, 8 Oct 2026): ONE-WAY, telemetry-anchored CWT model.
//
// CWT_weekly.stan (06w) with ONE change: CWT survival is no longer a fitted
// rate. It comes from telemetry mortality per day through Lorenzen's
// size-mortality relationship:
//   hazard_g = M_r * (L_r / L_g)^c * exp(delta + h_g)        (per day)
//   surv_g   = exp(-hazard_g * T_g)
// L_g = group mean fork length, L_r = 70 mm, T_g = the group's median travel
// time to KL (estimated here, as in 06w), delta = optional CWT offset
// (hatchery or tag differences), h_g = CWT group effect.
// M_r comes from the telemetry production fit and is FIXED (one-way: CWT
// cannot change telemetry):
//   M_r(c) = H_tel / sum_i T_i (L_r / L_i)^c
// H_tel = total telemetry hazard (sum of -log survival release to Sacramento),
// T_i, L_i = telemetry fish travel time and fork length (posterior medians).
// M_r depends on c, so it is recomputed as c is sampled.
// Unchanged from 06w: travel time model, weekly KL capture probability with
// BT-SPAS-X priors, the conditioning term for recovered travel times, log_lik.
// Simplifications: survival uses the median travel time (not averaged over
// the travel-time distribution); no flow term on CWT survival.

functions {
  // log P(t_lo < T < t_hi) for lognormal T with log-scale location mu and
  // scale sigma. Uses only std_normal_lcdf, which is accurate far into the
  // lower tail (std_normal_lccdf and lognormal_lccdf give -inf and an
  // infinite gradient beyond about 8 SDs). Left of the median: CDF
  // difference; right: upper-tail difference via lccdf(z) = lcdf(-z). Beyond
  // 20 SDs (mass < 1e-88) the density at the interval midpoint times its
  // width is used, so the log probability and its gradient stay finite.
  real log_lognormal_interval_prob(real t_lo, real t_hi, real mu, real sigma) {
    real z_hi = (log(t_hi) - mu) / sigma;
    real z_lo = t_lo > 0 ? (log(t_lo) - mu) / sigma : negative_infinity();
    if (z_hi < -20 || z_lo > 20) {
      return lognormal_lpdf(0.5 * (t_lo + t_hi) | mu, sigma) + log(t_hi - t_lo);
    }
    if (t_lo <= 0) {
      return std_normal_lcdf(z_hi);
    }
    if (z_hi <= 0) {
      return log_diff_exp(std_normal_lcdf(z_hi), std_normal_lcdf(z_lo));
    }
    return log_diff_exp(std_normal_lcdf(-z_lo), std_normal_lcdf(-z_hi));
  }
}

data {
  int Ncwtgrp;
  int Nyrs;
  int Nobs;
  int UseSizeEffect;
  array[Ncwtgrp] real Sz_cwt;
  array[Ncwtgrp] real CovX_cwt;
  array[Ncwtgrp] real ReachKM_cwt;

  vector[Nobs] ObsTT;
  array[Nobs] int TTind;

  array[Ncwtgrp] int Nrel;
  array[Ncwtgrp] int cwt_recaptures;

  int Nwk;
  vector[Nwk] lp_mu;
  vector<lower=0>[Nwk] lp_sd;
  int Nseg;
  array[Nseg] real seg_lo;
  array[Nseg] real seg_hi;
  array[Nseg] int seg_wk;
  array[Ncwtgrp] int seg_start;
  array[Ncwtgrp] int seg_end;
  array[Nobs] int obs_wk;

  // ONE-WAY: telemetry anchor (fixed) and Lorenzen terms
  real<lower=0> H_tel;              // total telemetry hazard, sum of -log(survival)
  int Ntel;
  vector<lower=0>[Ntel] T_tel;      // telemetry travel time release to Sacramento (days)
  vector[Ntel] log_Lratio_tel;      // log(L_r / L_i), telemetry fish
  vector[Ncwtgrp] log_Lratio_cwt;   // log(L_r / L_g), CWT groups
  real c_mu;                        // prior on the Lorenzen exponent c
  real<lower=0> c_sd;
  int<lower=0, upper=1> UseOffset;  // 1 = estimate a CWT offset on the hazard
  int<lower=0, upper=1> UseGamma;   // 1 = estimate the exposure exponent gamma
  real gamma_mu;                    // prior on gamma
  real<lower=0> gamma_sd;
}

parameters {
  real<lower=-5, upper=5> T_bReach;
  real<lower=-10, upper=10> TT_bCov;
  real<lower=-10, upper=10> T_bSz;
  vector[Ncwtgrp] z_TT;
  real<lower=0.001> sd_cwtTTre;
  vector[Nwk] z_pcap;
  real<lower=0.001, upper=5> Pro_sd;
  real<lower=0, upper=3> c;         // ONE-WAY: Lorenzen exponent
  real delta;                       // ONE-WAY: CWT offset on log hazard
  vector[Ncwtgrp] z_h;              // ONE-WAY: CWT group effects on log hazard
  real<lower=0.001> sd_h;
  real<lower=0, upper=2> gamma_raw; // exposure exponent (used when UseGamma = 1)
}

transformed parameters {
  real gam = UseGamma ? gamma_raw : 1.0;
  real log_Mr = log(H_tel) - log_sum_exp(gam * log(T_tel) + c * log_Lratio_tel);   // per day^gamma at L_r
  vector[Ncwtgrp] TT_cwtre = sd_cwtTTre * z_TT;
  vector[Nwk] logit_pcap = lp_mu + lp_sd .* z_pcap;
  vector[Ncwtgrp] mu_TT;
  vector[Ncwtgrp] log_hazard;
  vector[Ncwtgrp] surv_cwt;
  vector[Ncwtgrp] log_prec;
  vector[Ncwtgrp] pKL;

  for (i in 1:Ncwtgrp) {
    mu_TT[i] = T_bReach + UseSizeEffect * T_bSz * Sz_cwt[i] +
               TT_bCov * CovX_cwt[i] + TT_cwtre[i] + log(ReachKM_cwt[i] / 100);
    log_hazard[i] = log_Mr + c * log_Lratio_cwt[i] + UseOffset * delta + sd_h * z_h[i];
    surv_cwt[i] = exp(-exp(log_hazard[i] + gam * mu_TT[i]));
    {
      int n = seg_end[i] - seg_start[i] + 1;
      vector[n] lw;
      for (k in 1:n) {
        int s = seg_start[i] + k - 1;
        lw[k] = log_lognormal_interval_prob(seg_lo[s], seg_hi[s], mu_TT[i], Pro_sd) +
                log_inv_logit(logit_pcap[seg_wk[s]]);
      }
      log_prec[i] = log_sum_exp(lw);
    }
    pKL[i] = surv_cwt[i] * exp(log_prec[i]);
  }
}

model {
  T_bReach ~ normal(0, 1.5);
  TT_bCov ~ normal(0, 1);
  T_bSz ~ normal(0, 1);
  sd_cwtTTre ~ normal(0, 1);
  Pro_sd ~ normal(0, 1);
  z_TT ~ std_normal();
  z_pcap ~ std_normal();
  c ~ normal(c_mu, c_sd);
  delta ~ normal(0, 0.5);
  gamma_raw ~ normal(gamma_mu, gamma_sd);
  sd_h ~ normal(0, 0.5);
  z_h ~ std_normal();

  cwt_recaptures ~ binomial(Nrel, pKL);
  ObsTT ~ lognormal(mu_TT[TTind], Pro_sd);
  target += sum(log_inv_logit(logit_pcap[obs_wk]) - log_prec[TTind]);
}

generated quantities {
  real M_r = exp(log_Mr);
  vector[Ncwtgrp] pTT = exp(mu_TT);
  vector[Ncwtgrp] p_rec_given_surv = exp(log_prec);
  vector[Ncwtgrp] log_lik;
  for (i in 1:Ncwtgrp) {
    log_lik[i] = binomial_lpmf(cwt_recaptures[i] | Nrel[i], pKL[i]);
  }
}
