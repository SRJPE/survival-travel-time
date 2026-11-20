// Joint CJS (acoustic telemetry) + CWT (Knights Landing recoveries) model
// Draft for simulation testing, 30 Sep 2026.
//
// One survival and travel-time process for all fish, two observation models:
//   Telemetry: CJS over Woodson (W), Butte (B), Sacramento (S), Delta (D)
//   CWT:       daily recoveries in the Knights Landing (KL) rotary screw trap
// KL lies inside reach 3 (Butte to Sacramento).
//
// Survival per 100 km in reach j (logit scale):
//   eta_j = b[j] + bQ*Q + bL*L + group RE + delta_AT (acoustic fish only)
//   S over a reach of km_j = inv_logit(eta_j)^(km_j/100)
// Travel time: log days per 100 km tau_j = t[j] + gQ*Q + gL*L + group RE
//   cumulative median to a station = sum_j exp(tau_j) * km_j/100,
//   individual times lognormal around it (Flora's structure), on the log scale.
// KL trap: logit p_KL(day) = a0 + year RE + bp * daily KL flow.
// CWT recoveries by time bin: y[g,d] ~ Poisson(N_g * S_KL,g * P(arrive in bin d) * p_KL(d)),
// the time-stratified idea of Hance et al. (2021), at group level. Bins (default 5 days)
// instead of days keep the likelihood fast; p_KL uses the bin-mean flow.
//
// Switches select the four models compared in the simulation:
//   A CWT only, constant pCap with prior (like 06)
//   B Joint, constant-in-flow pCap, efficiency trials
//   C Joint, flow-dependent pCap, efficiency trials
//   D CJS only (acoustic offset not estimable, fixed at 0)

data {
  // model switches
  int<lower=0, upper=1> use_cjs;
  int<lower=0, upper=1> use_cwt;
  int<lower=0, upper=1> use_trials;
  int<lower=0, upper=1> p_flow;     // KL efficiency depends on daily flow
  int<lower=0, upper=1> p_year;     // year random effect on KL efficiency
  int<lower=0, upper=1> est_delta;  // estimate acoustic-tag survival offset
  int<lower=0, upper=1> pool_reaches; // one survival and travel-time rate for reaches 1-3
                                      // (CWT-only fit, like 06: reaches not separable)
  real prior_a0_mu;                 // prior on logit KL efficiency (BT-SPAS-X style)
  real<lower=0> prior_a0_sd;

  int<lower=1> Y;                   // years (all data sources)

  // telemetry release groups
  int<lower=0> G_t;
  array[G_t] int<lower=1, upper=Y> yr_t;
  vector[G_t] Q_t;                  // standardized flow covariate (group level)
  vector[G_t] L_t;                  // standardized fork length (group mean)
  matrix[G_t, 4] km_t;              // km in reaches 1-4 (release-W, W-B, B-S, S-D)
  array[G_t, 16] int<lower=0> n_hist;   // counts of detection histories over W,B,S,D
  array[G_t, 3] int<lower=0> tt_n;      // travel-time sufficient statistics, W,B,S
  matrix[G_t, 3] tt_s1;                 // sum of log travel times
  matrix[G_t, 3] tt_s2;                 // sum of squared log travel times

  // CWT release groups
  int<lower=0> G_c;
  int<lower=1> D;                   // recovery-time bins after release (e.g. 12 bins of 5 days)
  vector<lower=0>[D] edge;          // upper edge of each bin (days after release)
  array[G_c] int<lower=1, upper=Y> yr_c;
  vector[G_c] Q_c;
  vector[G_c] L_c;
  matrix[G_c, 3] km_c;              // km in reaches 1, 2 and the part of 3 above KL
  vector[G_c] N_c;                  // number released
  array[G_c, D] int<lower=0> y_c;   // KL recoveries per bin
  matrix[G_c, D] q_c;               // standardized KL flow, mean over each bin

  // KL efficiency trials (independent marked releases)
  int<lower=0> K;
  array[K] int<lower=1, upper=Y> yr_k;
  vector[K] q_k;
  array[K] int<lower=0> m_k;
  array[K] int<lower=0> r_k;

  // estimand: untagged fish from the Battle site, mean flow
  vector[3] km_est;
  real L_est;
}

transformed data {
  // detection history h = 1..16 as bits over stations W, B, S, D
  array[16, 4] int H;
  array[16] int last_det;           // last station detected (0 = none)
  array[4] int pw = {8, 4, 2, 1};
  for (h in 1:16) {
    int v = h - 1;
    last_det[h] = 0;
    for (k in 1:4) {
      H[h, k] = (v %/% pw[k]) % 2;
      if (H[h, k] == 1) last_det[h] = k;
    }
  }
  vector[D] log_day = log(edge);
}

parameters {
  // survival
  vector[4] b;
  real bQ;
  real bL;
  real delta_raw;
  real<lower=0> sd_S;
  vector[G_t] zS_t;
  vector[G_c] zS_c;
  // travel time
  vector[4] t;
  real gQ;
  real gL;
  real<lower=0> sd_T;
  vector[G_t] zT_t;
  vector[G_c] zT_c;
  real<lower=0> sigma_t;            // individual log travel-time SD, telemetry
  real<lower=0> sigma_c;            // individual log travel-time SD, CWT
  // telemetry detection by station and year
  vector[4] mu_p;
  vector<lower=0>[4] sd_p;
  matrix[Y, 4] zp;
  // KL trap efficiency
  real a0;
  real<lower=0> sd_a;
  vector[Y] za;
  real bp_raw;
}

transformed parameters {
  real delta_AT = est_delta * delta_raw;
  real bp = p_flow * bp_raw;
  vector[3] b3 = pool_reaches == 1 ? rep_vector(b[1], 3) : b[1:3];
  vector[3] t3 = pool_reaches == 1 ? rep_vector(t[1], 3) : t[1:3];
}

model {
  // priors
  b ~ normal(0, 1.5);
  bQ ~ normal(0, 1);
  bL ~ normal(0, 1);
  delta_raw ~ normal(0, 1);
  sd_S ~ normal(0, 1);
  zS_t ~ std_normal();
  zS_c ~ std_normal();
  t ~ normal(1, 1);
  gQ ~ normal(0, 1);
  gL ~ normal(0, 1);
  sd_T ~ normal(0, 0.5);
  zT_t ~ std_normal();
  zT_c ~ std_normal();
  sigma_t ~ normal(0, 1);
  sigma_c ~ normal(0, 1);
  mu_p ~ normal(1, 1.5);
  sd_p ~ normal(0, 1);
  to_vector(zp) ~ std_normal();
  a0 ~ normal(prior_a0_mu, prior_a0_sd);
  sd_a ~ normal(0, 0.5);
  za ~ std_normal();
  bp_raw ~ normal(0, 1);

  // telemetry: CJS on grouped detection histories + travel times
  if (use_cjs == 1) {
    for (g in 1:G_t) {
      vector[4] eta = b + bQ * Q_t[g] + bL * L_t[g] + sd_S * zS_t[g] + delta_AT;
      vector[4] logphi;
      vector[4] p;
      vector[5] chi;                // chi[k+1]: P(never seen after occasion k)
      vector[4] tau = t + gQ * Q_t[g] + gL * L_t[g] + sd_T * zT_t[g];
      for (j in 1:4) {
        logphi[j] = km_t[g, j] / 100 * log_inv_logit(eta[j]);
        p[j] = inv_logit(mu_p[j] + sd_p[j] * zp[yr_t[g], j]);
      }
      chi[5] = 1;
      for (r in 1:4) {
        int k = 5 - r;              // occasion k = 4..1 (station index), k = 0 release
        chi[k] = (1 - exp(logphi[k])) + exp(logphi[k]) * (1 - p[k]) * chi[k + 1];
      }
      for (h in 1:16) {
        if (n_hist[g, h] > 0) {
          real ll = 0;
          for (k in 1:last_det[h]) {
            ll += logphi[k];
            ll += H[h, k] == 1 ? log(p[k]) : log1m(p[k]);
          }
          ll += log(chi[last_det[h] + 1]);
          target += n_hist[g, h] * ll;
        }
      }
      // cumulative median travel time to W, B, S
      {
        vector[4] lseg = tau + log(km_t[g]' / 100);
        for (k in 1:3) {
          if (tt_n[g, k] > 0) {
            real mu = log_sum_exp(lseg[1:k]);
            target += -tt_n[g, k] * log(sigma_t)
                      - (tt_s2[g, k] - 2 * mu * tt_s1[g, k] + tt_n[g, k] * square(mu))
                        / (2 * square(sigma_t));
          }
        }
      }
    }
  }

  // CWT: daily recoveries at KL
  if (use_cwt == 1) {
    for (g in 1:G_c) {
      vector[3] eta = b3 + bQ * Q_c[g] + bL * L_c[g] + sd_S * zS_c[g];
      real logS = dot_product(km_c[g]' / 100, log_inv_logit(eta));
      vector[3] tau = t3 + gQ * Q_c[g] + gL * L_c[g] + sd_T * zT_c[g];
      real mu = log_sum_exp(tau + log(km_c[g]' / 100));
      vector[D] F = Phi((log_day - mu) / sigma_c);
      vector[D] parr = F - append_row(0, F[1:(D - 1)]);
      vector[D] logit_p = a0 + p_year * sd_a * za[yr_c[g]] + bp * q_c[g]';
      vector[D] log_lambda = log(N_c[g]) + logS + log(parr + 1e-12)
                             + log_inv_logit(logit_p);
      y_c[g] ~ poisson_log(log_lambda);
    }
  }

  // KL efficiency trials
  if (use_trials == 1 && K > 0) {
    r_k ~ binomial_logit(m_k, a0 + p_year * sd_a * za[yr_k] + bp * q_k);
  }
}

generated quantities {
  // survival from the Battle site to KL, untagged fish, typical group, mean flow
  real S_KL = exp(dot_product(km_est / 100, log_inv_logit(b3 + bL * L_est)));
  // same for an acoustic-tagged fish (what a CJS-only fit estimates)
  real S_KL_tagged = exp(dot_product(km_est / 100,
                                     log_inv_logit(b3 + bL * L_est + delta_AT)));
  // median travel time to KL (days)
  real TT_KL = exp(log_sum_exp(t3 + gL * L_est + log(km_est / 100)));
  real p_KL_mean = inv_logit(a0);
}
