// sim_tel_lorenzen.stan: telemetry fit for the Lorenzen one-way simulation.
// Stands in for the production telemetry model: survival logit-linear in fork
// length (as Flora's), lognormal travel time with a length effect, detection
// at the end of the reach known (p). The simulated truth is Lorenzen per-day
// mortality, so this model is deliberately not the truth, as in real data.
data {
  int<lower=1> N;
  vector[N] z;                       // standardized fork length
  array[N] int<lower=0, upper=1> det;
  real<lower=0, upper=1> p;
  int<lower=0> Nobs;
  vector<lower=0>[Nobs] T_obs;       // travel times of detected fish (days)
  array[Nobs] int<lower=1, upper=N> idx;
}
parameters {
  real a;
  real b;
  real t0;
  real t1;
  real<lower=0> sig;
}
model {
  a ~ normal(0, 2);
  b ~ normal(0, 1);
  t0 ~ normal(log(10), 1.5);
  t1 ~ normal(0, 1);
  sig ~ normal(0, 1);
  det ~ bernoulli(inv_logit(a + b * z) * p);
  T_obs ~ lognormal(t0 + t1 * z[idx], sig);
}
generated quantities {
  real H = -sum(log_inv_logit(a + b * z));   // total telemetry hazard, as 06o
}
