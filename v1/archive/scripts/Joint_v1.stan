// ---------------------------------------------------------------------------
// Joint_v1.stan (bp, 7 Oct 2026): joint telemetry + CWT model, version 1.
//
// Telemetry part = CovIndCont_ncRE.stan (03f) UNCHANGED: Flora's model with
// non-centered release-group effects. Everything below marked "CWT (joint v1)"
// is added; nothing above it is edited.
//
// CWT part = CWT_weekly.stan (06w) with ONE change: CWT fish no longer have
// their own survival and travel-time parameters. Following Perry et al.
// (exchangeability), they share the telemetry reach parameters:
//   survival to KL = prod over reaches 1-3 of
//     inv_logit(S_bReach[j] + S_bCov*flow + S_bSz*size + S_cwtre[g])^(km_gj/100)
//   median travel time to KL = sum over reaches 1-3 of
//     exp(T_bReach[j] + T_bSz*size + TT_bCov*flow + TT_cwtre[g]) * km_gj/100
// with km_gj the CWT group's distance in each telemetry reach (KL lies about
// 119 km into reach 3, Butte-Sac). Shared: S_bReach[1:3], S_bCov, S_bSz,
// T_bReach[1:3], T_bSz, TT_bCov (Flora's logit-linear size effect).
// CWT-only: release-group effects (S_cwtre, TT_cwtre; one per group, own SDs),
// observation SD of travel time (Pro_sd_cwt), and weekly KL capture
// probability with BT-SPAS-X priors, as in 06w.
// CWT flow and size are standardized on the telemetry scales (flow: A2
// Knights Landing monthly-table reference; size: musz, sdsz).
// log_lik (telemetry) is as Flora's, so 05 and LOO work unchanged;
// log_lik_cwt is separate.
// ---------------------------------------------------------------------------
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

//CJS model taken from stan example at https://mc-stan.org/docs/stan-users-guide/mark-recapture-models.html
//This version uses year-reach specific pCaps. 
// Survival influenced by fish size at release or 2 or 3-level water year type fixed effect or both cumulatively
//and random effect for year-release group
//
// ---------------------------------------------------------------------------
// CovIndCont_ncRE.stan (bp, 6 Oct 2026): Flora's CovIndCont.stan with ONE
// change: the release-group random effects are non-centered (S_RE, TT_RE,
// S_REt, TT_RET). Lines that changed are marked "CHANGE (non-centered)".
// The detection effects P_b stay exactly as hers (centered, bounded -10 to 10).
//
// Why: with release groups defined as study x release month (02, 6 Oct),
// some groups are small and some travel-time SDs are near zero
// (TTRE_sd[4] about 0.1). Her centered parameterization then forms a funnel:
// in the 6 Oct fit, chain 1 stuck there (39 divergences, E-BFMI 0.15,
// TTRE_sd[4] Rhat 1.19, ESS 13). Non-centering removes the funnel.
// Same model, same posterior, same priors (flat on the SDs, as hers), same
// parameter names, so 05 reads the fit unchanged. Her -10 to 10 bounds on the
// four random-effect arrays are dropped (they are now derived quantities);
// the effects are far inside those bounds in every fit so far.
// ---------------------------------------------------------------------------
data {
  int Nind;
  int Nreaches;
  int Ndetlocs;
  int Nyrs;
  int Nrg;
  int UseSizeEffect; // turn on and off the fish size effect
  int Ntribs;//# of tributaries where survival is modelled (Butte + Feather = 2)
  int NindT;// # of individuals released from tributaries summed across Butte and Feather
  int NreachesT;
  int NrgT;// # of release groups summed across Butte and Feather
  vector[NindT] RmultT;//The stream length multiplier for survival (based on distance from tributary release to sac station) for each individual 
  array[Nind,Nreaches] real Rmult; // changed this to work for new Rmult
  array[Nind,Ndetlocs] int CH;
  array[Nind] int yrind;
  array[Nind] int rgind;
  array[Nreaches] int rch_covind;
  array[Nind] int firstCap;
  array[Nind] int lastCap;
  array[NindT] int trib_ind; //The tributary index for each fish
  array[NindT] int firstCapT;// the first detection station of each tributary fish (including release, thus firtCapT = 1 in all cases)
  array[NindT] int lastCapT;// the last detection station for each tributary fish. Note for tributary fish this is 3 or less (Delta = 3 for trib fish)
  array[NindT] int yrindT;// The annual index for each tributary fish for mainstem station pCap mapping
  array[NindT] int rgindT;//The tributary release group index for each tributary fish
  array[NindT,3] int CHT;//The capture history for each tributary fish (digit 1 = release, digit 2 = sac station, digit 3 = delta stations)
  array[Nind] real Sz;
  vector[NindT] SzT;
  array[Nind,Nreaches] real CovX; 
  array[NindT,NreachesT] real CovXT;
  
  int Nobs;// # number of individuals*observed detection locations, which allows travel time from release to detection location
  vector[Nobs] ObsTT;
  array[Nind,Nreaches] int TTind; // Index which specifies the element of ObsTT given the 'i' individual in the jth reach.
  int NobsT;// # number of individuals*observed detection locations, which allows travel time from release to detection location
  vector[NobsT] ObsTTT;
  array[NindT,NreachesT] int TTindT; // Index which specifies the element of ObsTT given the 'i' individual in the jth reach.
  array[Nind,Nreaches] real ReachKM_ind;
  array[NindT,NreachesT] real ReachKMT_ind;

  // ---- CWT (joint v1) ----
  int Ncwtgrp;
  int Nobs_cwt;
  array[Ncwtgrp,3] real Rmult_cwt;   // km/100 in reaches 1-3, release to KL
  array[Ncwtgrp] real Sz_cwt;        // group mean fork length, telemetry size scale
  array[Ncwtgrp] real CovX_cwt;      // KL monthly max flow in release month, telemetry flow scale
  vector[Nobs_cwt] ObsTT_cwt;        // days from release to recovery at KL
  array[Nobs_cwt] int TTind_cwt;     // group of each recovery
  array[Ncwtgrp] int Nrel;
  array[Ncwtgrp] int cwt_recaptures;
  int Nwk;                           // KL year x run weeks with a capture probability
  vector[Nwk] lp_mu;                 // BT-SPAS-X posterior mean of logit pCap
  vector<lower=0>[Nwk] lp_sd;        // BT-SPAS-X posterior SD of logit pCap
  int Nseg;                          // group x week segments, days since release
  array[Nseg] real seg_lo;
  array[Nseg] real seg_hi;
  array[Nseg] int seg_wk;
  array[Ncwtgrp] int seg_start;
  array[Ncwtgrp] int seg_end;
  array[Nobs_cwt] int obs_wk;        // week of each recovery
}

parameters {
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] P_b;// reach -specific detection probability
  vector <lower=-10,upper=10> [Nreaches] muPb; // mean detection probability
  vector <lower=0.001> [Nreaches] sdPb; // standard deviation of detection probability
  vector <lower=-5,upper=5> [Nreaches] S_bReach; //survival by reach under baseline covariate
  vector <lower=-5,upper=5> [Nreaches] T_bReach; //survival by reach under baseline covariate
  real <lower=-10,upper=10>  S_bCov; //sac mainstem environmental covariate
  real <lower=-10,upper=10> TT_bCov; //sac mainstem environmental covariate for travel time
  matrix[Nrg, Nreaches] z_S; // CHANGE (non-centered): replaces S_RE, which is now a transformed parameter
  vector <lower=0.001> [Nreaches] RE_sd; //sac mainstem SD of random effect
  real <lower=-10,upper=10> S_bSz; // fish size covariate
  real <lower=-10,upper=10> T_bSz; // fish size covariate
  array[Ntribs] real <lower=-10,upper=10> S_bTrib;
  array[Ntribs] real <lower=-10,upper=10> T_bTrib;
  matrix[NrgT, NreachesT] z_St; // CHANGE (non-centered): replaces S_REt
  vector <lower=0.001> [NreachesT] RE_sdT; //tribSD of random effect
  array[Ntribs] real <lower=-10,upper=10> S_bCovT;//trib environmental covariate
  array[Ntribs] real <lower=-10,upper=10> TT_bCovT;//trib environmental covariate
   
  real <lower=0.001, upper=5> Pro_sd;
  matrix[Nrg, Nreaches] z_TT; // CHANGE (non-centered): replaces TT_RE
  vector <lower=0.001> [Nreaches] TTRE_sd;
  real <lower=0.001, upper=5> Pro_sdT;
  matrix[NrgT, NreachesT] z_TTt; // CHANGE (non-centered): replaces TT_RET
 vector <lower=0.001> [NreachesT] TTRE_sdT;

  // ---- CWT (joint v1) ----
  vector[Ncwtgrp] z_Scwt;
  real<lower=0.001> sd_cwtSre;
  vector[Ncwtgrp] z_TTcwt;
  real<lower=0.001> sd_cwtTTre;
  real<lower=0.001, upper=5> Pro_sd_cwt;
  vector[Nwk] z_pcap;
}

transformed parameters{
  // CHANGE (non-centered): rebuild the random effects from the standard normals.
  // Same names and indexing as Flora's parameters, so the code below and 05 are unchanged.
  matrix[Nrg, Nreaches] S_RE = diag_post_multiply(z_S, RE_sd);        // column j = RE_sd[j] * z_S[,j]
  matrix[Nrg, Nreaches] TT_RE = diag_post_multiply(z_TT, TTRE_sd);
  matrix[NrgT, NreachesT] S_REt = diag_post_multiply(z_St, RE_sdT);
  matrix[NrgT, NreachesT] TT_RET = diag_post_multiply(z_TTt, TTRE_sdT);
  // end CHANGE

  array[Nind] vector[Nreaches] surv; //survival by reach
  array[Nind] vector[Ndetlocs] Pcap; //detection probability by detection location
  array[Nind] vector[Ndetlocs] chi; //probability of NOT detecting individual again by station
  array[NindT] vector[NreachesT] survT; // two reaches release in trib - sac, sac-delta
  array[NindT] vector[3] PcapT;//3 locations: release (not calculated as assumed = 1, 2 = sac, 3 = delta)
  array[NindT] vector[3] chiT; //probability of not detecting individual by station
  
  array[Nind] vector[Nreaches] pTT; //Cummulative travel time from release to detection location
  vector[Nobs] lg_pTT;
  array[NindT] vector[NreachesT] pTTT;  //Cummulative travel time from release to detection location
  vector[NobsT] lg_pTTT;
  real TT;
  
  // Sacramento River fish
  for(i in 1:Nind){
    TT=exp(T_bReach[1] + T_bSz*Sz[i] + TT_bCov*CovX[i,1] + TT_RE[rgind[i],1]);
    pTT[i,1]=TT*ReachKM_ind[i,1]/100;//predicted travel time is just days/km * km from release to Woodson
    if(TTind[i,1]>0) lg_pTT[TTind[i,1]]=log(pTT[i,1]);//map log predictions to vector for lognormal likelihood calc

    for(j in 2:Nreaches){
      TT=exp(T_bReach[j] + T_bSz*Sz[i] + TT_bCov*CovX[i,j] + TT_RE[rgind[i],j]);
      pTT[i,j] = pTT[i,j-1] + TT*ReachKM_ind[i,j]/100;//accumulate travel times in downstream direction
      if(TTind[i,j]>0)lg_pTT[TTind[i,j]]=log(pTT[i,j]);
    }
  }
  
  for(i in 1:Nind){
    for(j in 1:Ndetlocs){  
      if(j<Ndetlocs) {
      surv[i,j]=inv_logit(S_bReach[rch_covind[j]] + S_bCov*CovX[i,j] + UseSizeEffect*S_bSz*Sz[i] + S_RE[rgind[i],j])^Rmult[i,j]; 
     }
      if(j>1) Pcap[i,j]=inv_logit(P_b[yrind[i],j-1]);
    }
    chi[i,Ndetlocs]=1.0;
    for(j in 1:Nreaches){
      int r_curr; int r_next;
      r_curr=Ndetlocs-j;  r_next=r_curr+1;
      //            not surviving   or  surviving but not decteced        Accumulate across reaches
      chi[i,r_curr]=(1-surv[i,r_curr]) + surv[i,r_curr]*(1-Pcap[i,r_next]) * chi[i,r_next];        
    }
  }//Nind loop
  
  //Tributary fish (Butte and Feather)
  real TTT;

  for(i in 1:NindT){
    //survival and travel time rates in tributary until SAC detection location
    survT[i,1]=inv_logit(S_bTrib[trib_ind[i]] + S_bCovT[trib_ind[i]]*CovXT[i,1] + UseSizeEffect*S_bSz*SzT[i]+ S_REt[rgindT[i],1])^RmultT[i];//RmultT is distance from release point to Sac station
    TTT = exp(T_bTrib[trib_ind[i]] + T_bSz*SzT[i] + TT_bCovT[trib_ind[i]] *CovXT[i,1] + TT_RET[rgindT[i],1]);
    pTTT[i,1]=TTT*ReachKMT_ind[i,1]/100;//predicted travel time is just days/km * km from release to Sacramento
    if(TTindT[i,1]>0) lg_pTTT[TTindT[i,1]]=log(pTTT[i,1]);//map log predictions to vector for lognormal likelihood calc
    
    //survival and travel time rates in Sac-Delta reach of mainstem
    int j=NreachesT;
    survT[i,2]=inv_logit(S_bReach[Nreaches] + S_bCov*CovXT[i,j] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i],j])^Rmult[1,Nreaches];//Same fixed effect for mainstem and tributary fish in mainstem
    TTT=exp(T_bReach[Nreaches] + T_bSz*SzT[i] + TT_bCov*CovXT[i,j] + TT_RET[rgindT[i],j]);
    pTTT[i,2] = pTTT[i,1] + TTT*ReachKMT_ind[i,2]/100;//accumulate travel times in downstream direction
    if(TTindT[i,2]>0)lg_pTTT[TTindT[i,2]]=log(pTTT[i,2]);
      
    j=4;//pCap at Sac station for tributary fish. Same as for mainstem fish. yrindT assigns the correct year for the trib fish 'i'
    PcapT[i,2]=inv_logit(P_b[yrindT[i],j-1]);//here j points to the mainstem station index, while 2 refers to the digit in the CHT sequence (Sac)
    j=5;//pCap at Delta station for tributary fish
    PcapT[i,3]=inv_logit(P_b[yrindT[i],j-1]);
    
    //Cummulative probability of not being detected by station for tributary fish
    chiT[i,3]=1.0;//probability of not being detected after delta station
    //probability of not being detected after sac station = not surviving from Sac-Delta or suriving Sac-Delta but not detected at Delta station
    chiT[i,2]=(1-survT[i,2]) + survT[i,2]*(1-PcapT[i,3])*chiT[i,3];
    //probability of not being detected after release = not surviving from release in trib to sac station or surviving release-sac but not detected at sac station * probability for next station computed above
    chiT[i,1]=(1-survT[i,1]) + survT[i,1]*(1-PcapT[i,2])*chiT[i,2];
  }

  // ---- CWT (joint v1) ----
  vector[Ncwtgrp] S_cwtre = sd_cwtSre * z_Scwt;
  vector[Ncwtgrp] TT_cwtre = sd_cwtTTre * z_TTcwt;
  vector[Nwk] logit_pcap = lp_mu + lp_sd .* z_pcap;
  vector[Ncwtgrp] surv_cwt;   // survival release to KL
  vector[Ncwtgrp] mu_TTcwt;   // log median travel time release to KL
  vector[Ncwtgrp] log_prec;   // log P(recovered at KL | survived)
  vector[Ncwtgrp] pKL;
  for (g in 1:Ncwtgrp) {
    real ls = 0;
    real tt = 0;
    for (j in 1:3) {
      ls += Rmult_cwt[g,j] * log_inv_logit(S_bReach[rch_covind[j]] + S_bCov*CovX_cwt[g] +
                                           UseSizeEffect*S_bSz*Sz_cwt[g] + S_cwtre[g]);
      tt += exp(T_bReach[j] + T_bSz*Sz_cwt[g] + TT_bCov*CovX_cwt[g] + TT_cwtre[g]) * Rmult_cwt[g,j];
    }
    surv_cwt[g] = exp(ls);
    mu_TTcwt[g] = log(tt);
    {
      int n = seg_end[g] - seg_start[g] + 1;
      vector[n] lw;
      for (k in 1:n) {
        int s = seg_start[g] + k - 1;
        lw[k] = log_lognormal_interval_prob(seg_lo[s], seg_hi[s], mu_TTcwt[g], Pro_sd_cwt) +
                log_inv_logit(logit_pcap[seg_wk[s]]);
      }
      log_prec[g] = log_sum_exp(lw);
    }
    pKL[g] = surv_cwt[g] * exp(log_prec[g]);
  }
}

model {
    ObsTT[]~lognormal(lg_pTT[],Pro_sd);
   to_vector(z_TT) ~ std_normal(); // CHANGE (non-centered): was for(j in 1:Nreaches){TT_RE[,j]~normal(0,TTRE_sd[j]);}

   ObsTTT[]~lognormal(lg_pTTT[],Pro_sdT);
   to_vector(z_TTt) ~ std_normal(); // CHANGE (non-centered): was for(j in 1:NreachesT){ TT_RET[,j] ~normal(0,TTRE_sdT[j]);}
    
  for(i in 1:Nind){
    for(j in (firstCap[i]+1):lastCap[i]){  //Loop through all detection stations after first to last station individual was detected at
    1~bernoulli(surv[i,j-1]);//had to be alive prior to last detection station
    CH[i,j]~bernoulli(Pcap[i,j]);
    }
    1~bernoulli(chi[i,lastCap[i]]);//probability of individual never beeing seen again after last detection
  }
  

  // CHANGE (non-centered): was for (j in 1:Nreaches){ S_RE[,j] ~ normal(0,RE_sd[j]); }
  to_vector(z_S) ~ std_normal();

  for(j in 1:Nreaches){
    P_b[1:Nyrs,j]~normal(muPb[j],sdPb[j]);
    }
  
  for(i in 1:NindT){
    for(j in (firstCapT[i]+1):lastCapT[i]){//firstCapT is always 1, so +1 means loop starts at 2, and can end at 2 (if lastCapT=2 = Sac) or 3 (delta) 
    1~bernoulli(survT[i,j-1]);//had to be alive prior to last detection station
    CHT[i,j]~bernoulli(PcapT[i,j]);//note j can be 2 or 3
    }
    1~bernoulli(chiT[i,lastCapT[i]]);//probability of individual never beeing seen again after last detection
  }
  
  // CHANGE (non-centered): was for (j in 1:NreachesT){ S_REt[,j] ~ normal(0,RE_sdT[j]); }
  to_vector(z_St) ~ std_normal();

  // ---- CWT (joint v1) ----
  sd_cwtSre ~ normal(0, 1);
  sd_cwtTTre ~ normal(0, 1);
  Pro_sd_cwt ~ normal(0, 1);
  z_Scwt ~ std_normal();
  z_TTcwt ~ std_normal();
  z_pcap ~ std_normal();                       // logit_pcap ~ normal(lp_mu, lp_sd)
  cwt_recaptures ~ binomial(Nrel, pKL);
  ObsTT_cwt ~ lognormal(mu_TTcwt[TTind_cwt], Pro_sd_cwt);
  target += sum(log_inv_logit(logit_pcap[obs_wk]) - log_prec[TTind_cwt]);   // travel times given recovery
}


generated quantities{

  vector[Nind + NindT] log_lik;
  array[Nyrs] vector[Nreaches] pred_pcap;
  array[Nind] vector[Nreaches] pred_surv;
  array[NindT] vector[NreachesT] pred_survT;
  vector[Nind] SurvRelSac;
  vector[Nind] TT_RelSac;
  vector[NindT] TT_RelSacT;

  for (iyr in 1:Nyrs) {
    for (j in 1:Nreaches) {
      pred_pcap[iyr,j] = inv_logit(P_b[iyr,j]);
    }
  }

  for (i in 1:Nind) {
    pred_surv[i] = surv[i];
    SurvRelSac[i] = prod(pred_surv[i,1:3]);
    TT_RelSac[i] = pTT[i,3];
  }

  for (i in 1:NindT) {
    pred_survT[i] = survT[i];
    TT_RelSacT[i] = pTTT[i,NreachesT];
  }

  for(i in 1:Nind){
    log_lik[i]=0;
    for(j in (firstCap[i]+1):lastCap[i]){  //Loop through all detection stations after first to last station individual was detected at
    log_lik[i] += log(surv[i,j-1]);//had to be alive prior to last detection station
    if (CH[i,j]==1){
      log_lik[i] += log(Pcap[i,j]);
    } else {
      log_lik[i] += log1m(Pcap[i,j]);
    }
    }
    log_lik[i] += log(chi[i,lastCap[i]]);//probability of individual never beeing seen again after last detection
  }

  int k = Nind;
  for(i in 1:NindT){
        k=k+1;
        log_lik[k]=0;
    for(j in (firstCapT[i]+1):lastCapT[i]){  //Loop through all detection stations after first to last station individual was detected at
    log_lik[k] += log(survT[i,j-1]);//had to be alive prior to last detection station
    if (CHT[i,j]==1){
      log_lik[k] += log(PcapT[i,j]);
    } else {
      log_lik[k] += log1m(PcapT[i,j]);
    }
    }
     log_lik[k] += log(chiT[i,lastCapT[i]]);//probability of individual never beeing seen again after last detection
  }
  

  // ---- CWT (joint v1) ----
  vector[Ncwtgrp] log_lik_cwt;
  vector[Ncwtgrp] p_rec_given_surv = exp(log_prec);
  vector[Ncwtgrp] pTT_cwt = exp(mu_TTcwt);   // median travel time to KL (days)
  for (g in 1:Ncwtgrp) {
    log_lik_cwt[g] = binomial_lpmf(cwt_recaptures[g] | Nrel[g], pKL[g]);
  }
}

