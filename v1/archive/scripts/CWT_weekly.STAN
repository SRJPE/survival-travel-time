// CWT survival and travel-time model with WEEKLY Knights Landing
// capture probability from BT-SPAS-X (Oct 2026, bp).
//
// This is CWT_nc.stan with ONE change: the single capture probability
// (logit_pCap_Sim, one value for every group and week, prior from the
// pooled BT-SPAS-X trials) is replaced by a capture probability for
// each KL year x run week, with a prior from that week's BT-SPAS-X
// posterior (mean and SD of logit pCap). Weeks are shared by every
// group passing KL that week.
//
// A surviving fish of group i passes KL in week w with probability
// f_iw, from the group's own lognormal travel-time distribution
// (median pTT[i], SD Pro_sd, as in CWT_nc.stan). So
//
//   P(recovered | survived) = sum_w f_iw * pcap[w]
//   pKL[i] = surv_cwt[i] * P(recovered | survived)
//
// Travel time: recovered fish are over-represented in weeks when the
// trap catches well, so each observed travel time gets the term
//   log pcap[week of recovery] - log P(recovered | survived)
// on top of CWT_nc.stan's lognormal. Together with the binomial this
// is the full likelihood of the recoveries and their dates (bp's
// simple joint model used the same idea with weekly counts).
//
// Weeks after the end of the BT-SPAS-X season (run week 31, Julian
// week 22) are not modelled; recoveries in them are left out in 06w.
// Unchanged from CWT_nc.stan: survival and travel-time structure,
// non-centered group effects, priors, log_lik (binomial per group).

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
  array[Ncwtgrp] real Rmult_cwt;
  array[Ncwtgrp] real Sz_cwt;
  array[Ncwtgrp] real CovX_cwt;
  array[Ncwtgrp] real ReachKM_cwt;

  vector[Nobs] ObsTT;
  array[Nobs] int TTind;

  array[Ncwtgrp] int Nrel;
  array[Ncwtgrp] int cwt_recaptures;

  // NEW: weekly KL capture probability (BT-SPAS-X), one per KL year x run week used
  int Nwk;
  vector[Nwk] lp_mu;            // posterior mean of logit pCap
  vector<lower=0>[Nwk] lp_sd;   // posterior SD of logit pCap

  // NEW: each group's weeks at KL, as intervals of days since release
  int Nseg;
  array[Nseg] real seg_lo;
  array[Nseg] real seg_hi;
  array[Nseg] int seg_wk;       // index into lp_mu / lp_sd
  array[Ncwtgrp] int seg_start; // first and last segment of each group
  array[Ncwtgrp] int seg_end;
  array[Nobs] int obs_wk;       // week (index into lp_mu) of each recovery
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
  vector[Nwk] z_pcap;           // NEW (replaces logit_pCap_Sim)
  real<lower=0.001, upper=5> Pro_sd;
}

transformed parameters {
  vector[Ncwtgrp] S_cwtre = sd_cwtSre * z_S;
  vector[Ncwtgrp] TT_cwtre = sd_cwtTTre * z_TT;
  vector[Nwk] logit_pcap = lp_mu + lp_sd .* z_pcap;
  array[Ncwtgrp] real surv_cwt;
  array[Ncwtgrp] real pTT;
  vector[Ncwtgrp] mu_TT;        // log of median travel time (log pTT), kept on the log scale
  vector[Ncwtgrp] log_prec;     // log P(recovered at KL | survived)
  array[Ncwtgrp] real pKL;

  for (i in 1:Ncwtgrp) {
    surv_cwt[i] = inv_logit(S_bReach + UseSizeEffect * S_bSz * Sz_cwt[i] +
                            S_bCov * CovX_cwt[i] + S_cwtre[i])^Rmult_cwt[i];
    mu_TT[i] = T_bReach + UseSizeEffect * T_bSz * Sz_cwt[i] +
               TT_bCov * CovX_cwt[i] + TT_cwtre[i] + log(ReachKM_cwt[i] / 100);
    pTT[i] = exp(mu_TT[i]);   // same as CWT_nc.stan
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
  z_pcap ~ std_normal();        // logit_pcap ~ normal(lp_mu, lp_sd), from BT-SPAS-X

  cwt_recaptures ~ binomial(Nrel, pKL);

  ObsTT ~ lognormal(mu_TT[TTind], Pro_sd);   // as CWT_nc.stan, log(pTT) = mu_TT
  // NEW: travel times of recovered fish, given recovery
  target += sum(log_inv_logit(logit_pcap[obs_wk]) - log_prec[TTind]);
}

generated quantities {
  vector[Ncwtgrp] log_lik;
  vector[Ncwtgrp] p_rec_given_surv = exp(log_prec);
  for (i in 1:Ncwtgrp) {
    log_lik[i] = binomial_lpmf(cwt_recaptures[i] | Nrel[i], pKL[i]);
  }
}
