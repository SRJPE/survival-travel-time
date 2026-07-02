// CWT survival and travel time model
data {
  int Ncwtgrp;
  int Nyrs;
  int Nobs;// # number of individuals*observed detection locations, which allows travel time from release to detection location
  int UseSizeEffect; // turn on and off the fish size effect
  array[Ncwtgrp] real Rmult_cwt; 
  array[Ncwtgrp] real Sz_cwt;
  array[Ncwtgrp] real CovX_cwt; 
  array[Ncwtgrp] real ReachKM_cwt;

  vector[Nobs] ObsTT;
  array[Nobs] int TTind; // Index which specifies the element of ObsTT given the 'i' individual in the jth reach.
 
  array[Ncwtgrp] int Nrel;
  array[Ncwtgrp] int cwt_recaptures;
  
  real <lower=-10,upper=10> lt_mu_pCap;
  real <lower=-10,upper=10> lt_sd_pCap;
 }

parameters {
  real<lower=-5, upper=5> S_bReach;
  real<lower=-5, upper=5> T_bReach;
  real<lower=-10, upper=10> S_bCov;
  real<lower=-10, upper=10> TT_bCov;
  real<lower=-10, upper=10> S_bSz;
  real<lower=-10, upper=10> T_bSz;
  vector[Ncwtgrp] S_cwtre;
  real<lower=0.001> sd_cwtSre;
  vector[Ncwtgrp] TT_cwtre;
  real<lower=0.001> sd_cwtTTre;
  real logit_pCap_Sim;
  real<lower=0.001, upper=5> Pro_sd;
}

transformed parameters{
  array[Ncwtgrp] real surv_cwt;
  array[Ncwtgrp] real pKL;
  array[Ncwtgrp] real pTT;
  real TT;
  
  for(i in 1:Ncwtgrp){
    surv_cwt[i]=inv_logit(S_bReach +  UseSizeEffect*S_bSz*Sz_cwt[i] +
                S_bCov*CovX_cwt[i] +
                 S_cwtre[i])^Rmult_cwt[i];
    pKL[i]= surv_cwt[i] * inv_logit(logit_pCap_Sim);
    
    TT  = exp(T_bReach + UseSizeEffect*T_bSz * Sz_cwt[i] + 
              TT_bCov * CovX_cwt[i] + 
              TT_cwtre[i]);
    pTT[i] = TT * ReachKM_cwt[i] / 100;

    }
}

model {
  
  // Priors
  S_bReach ~ normal(0, 1.5);
  T_bReach ~ normal(0, 1.5);
  S_bCov ~ normal(0, 1);
  TT_bCov ~ normal(0, 1);
  S_bSz ~ normal(0, 1);
  T_bSz ~ normal(0, 1);
  sd_cwtSre ~ normal(0, 1);
  sd_cwtTTre ~ normal(0, 1);
  Pro_sd ~ normal(0, 1);

  S_cwtre ~ normal(0, sd_cwtSre);
  TT_cwtre ~ normal(0, sd_cwtTTre);
  
  logit_pCap_Sim~normal(lt_mu_pCap, lt_sd_pCap);//from bt-spasx and read in as data

  // Survival likelihood
  for (i in 1:Ncwtgrp) {
    cwt_recaptures[i] ~ binomial(Nrel[i], pKL[i]);
  }

  // Travel time likelihood
  for (i in 1:Nobs) {
    ObsTT[i] ~ lognormal(log(pTT[TTind[i]]), Pro_sd);
  }
}

generated quantities{
  vector[Ncwtgrp] log_lik;

  for (i in 1:Ncwtgrp) {
    log_lik[i] = binomial_lpmf(cwt_recaptures[i] | Nrel[i], pKL[i]);
  }
}

