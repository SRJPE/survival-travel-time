// CovIndCont_sim.stan: the production model CovIndCont.stan, rewritten ONLY
// so that per-fish quantities (surv, Pcap, chi, travel-time predictions) are
// local variables in the model block instead of transformed parameters, and
// the per-fish generated quantities are dropped. The parameters, their
// bounds, the (absent) priors and the likelihood are identical to
// CovIndCont.stan, so the posterior is the same. This keeps each simulation
// fit's output small (CmdStan writes every transformed parameter to disk).
data {
  int Nind;
  int Nreaches;
  int Ndetlocs;
  int Nyrs;
  int Nrg;
  int UseSizeEffect;
  int Ntribs;
  int NindT;
  int NreachesT;
  int NrgT;
  vector[NindT] RmultT;
  array[Nind, Nreaches] real Rmult;
  array[Nind, Ndetlocs] int CH;
  array[Nind] int yrind;
  array[Nind] int rgind;
  array[Nreaches] int rch_covind;
  array[Nind] int firstCap;
  array[Nind] int lastCap;
  array[NindT] int trib_ind;
  array[NindT] int firstCapT;
  array[NindT] int lastCapT;
  array[NindT] int yrindT;
  array[NindT] int rgindT;
  array[NindT, 3] int CHT;
  array[Nind] real Sz;
  vector[NindT] SzT;
  array[Nind, Nreaches] real CovX;
  array[NindT, NreachesT] real CovXT;

  int Nobs;
  vector[Nobs] ObsTT;
  array[Nind, Nreaches] int TTind;
  int NobsT;
  vector[NobsT] ObsTTT;
  array[NindT, NreachesT] int TTindT;
  array[Nind, Nreaches] real ReachKM_ind;
  array[NindT, NreachesT] real ReachKMT_ind;
}

parameters {
  array[Nyrs] vector<lower=-10, upper=10>[Nreaches] P_b;
  vector<lower=-10, upper=10>[Nreaches] muPb;
  vector<lower=0.001>[Nreaches] sdPb;
  vector<lower=-5, upper=5>[Nreaches] S_bReach;
  vector<lower=-5, upper=5>[Nreaches] T_bReach;
  real<lower=-10, upper=10> S_bCov;
  real<lower=-10, upper=10> TT_bCov;
  array[Nrg, Nreaches] real<lower=-10, upper=10> S_RE;
  vector<lower=0.001>[Nreaches] RE_sd;
  real<lower=-10, upper=10> S_bSz;
  real<lower=-10, upper=10> T_bSz;
  array[Ntribs] real<lower=-10, upper=10> S_bTrib;
  array[Ntribs] real<lower=-10, upper=10> T_bTrib;
  array[NrgT, NreachesT] real<lower=-10, upper=10> S_REt;
  vector<lower=0.001>[NreachesT] RE_sdT;
  array[Ntribs] real<lower=-10, upper=10> S_bCovT;
  array[Ntribs] real<lower=-10, upper=10> TT_bCovT;

  real<lower=0.001, upper=5> Pro_sd;
  array[Nrg, Nreaches] real<lower=-10, upper=10> TT_RE;
  vector<lower=0.001>[Nreaches] TTRE_sd;
  real<lower=0.001, upper=5> Pro_sdT;
  array[NrgT, NreachesT] real<lower=-10, upper=10> TT_RET;
  vector<lower=0.001>[NreachesT] TTRE_sdT;
}

model {
  array[Nind] vector[Nreaches] surv;
  array[Nind] vector[Ndetlocs] Pcap;
  array[Nind] vector[Ndetlocs] chi;
  array[NindT] vector[NreachesT] survT;
  array[NindT] vector[3] PcapT;
  array[NindT] vector[3] chiT;
  array[Nind] vector[Nreaches] pTT;
  vector[Nobs] lg_pTT;
  array[NindT] vector[NreachesT] pTTT;
  vector[NobsT] lg_pTTT;

  // Sacramento River fish: travel time
  for (i in 1:Nind) {
    real TT = exp(T_bReach[1] + T_bSz * Sz[i] + TT_bCov * CovX[i, 1] + TT_RE[rgind[i], 1]);
    pTT[i, 1] = TT * ReachKM_ind[i, 1] / 100;
    if (TTind[i, 1] > 0) lg_pTT[TTind[i, 1]] = log(pTT[i, 1]);
    for (j in 2:Nreaches) {
      TT = exp(T_bReach[j] + T_bSz * Sz[i] + TT_bCov * CovX[i, j] + TT_RE[rgind[i], j]);
      pTT[i, j] = pTT[i, j - 1] + TT * ReachKM_ind[i, j] / 100;
      if (TTind[i, j] > 0) lg_pTT[TTind[i, j]] = log(pTT[i, j]);
    }
  }

  // Sacramento River fish: survival and detection
  for (i in 1:Nind) {
    for (j in 1:Ndetlocs) {
      if (j < Ndetlocs) {
        surv[i, j] = inv_logit(S_bReach[rch_covind[j]] + S_bCov * CovX[i, j]
                               + UseSizeEffect * S_bSz * Sz[i] + S_RE[rgind[i], j]) ^ Rmult[i, j];
      }
      if (j > 1) Pcap[i, j] = inv_logit(P_b[yrind[i], j - 1]);
    }
    chi[i, Ndetlocs] = 1.0;
    for (j in 1:Nreaches) {
      int r_curr = Ndetlocs - j;
      int r_next = r_curr + 1;
      chi[i, r_curr] = (1 - surv[i, r_curr]) + surv[i, r_curr] * (1 - Pcap[i, r_next]) * chi[i, r_next];
    }
  }

  // Butte and Feather fish
  for (i in 1:NindT) {
    real TTT;
    int j;
    survT[i, 1] = inv_logit(S_bTrib[trib_ind[i]] + S_bCovT[trib_ind[i]] * CovXT[i, 1]
                            + UseSizeEffect * S_bSz * SzT[i] + S_REt[rgindT[i], 1]) ^ RmultT[i];
    TTT = exp(T_bTrib[trib_ind[i]] + T_bSz * SzT[i] + TT_bCovT[trib_ind[i]] * CovXT[i, 1] + TT_RET[rgindT[i], 1]);
    pTTT[i, 1] = TTT * ReachKMT_ind[i, 1] / 100;
    if (TTindT[i, 1] > 0) lg_pTTT[TTindT[i, 1]] = log(pTTT[i, 1]);

    j = NreachesT;
    survT[i, 2] = inv_logit(S_bReach[Nreaches] + S_bCov * CovXT[i, j]
                            + UseSizeEffect * S_bSz * SzT[i] + S_REt[rgindT[i], j]) ^ Rmult[1, Nreaches];
    TTT = exp(T_bReach[Nreaches] + T_bSz * SzT[i] + TT_bCov * CovXT[i, j] + TT_RET[rgindT[i], j]);
    pTTT[i, 2] = pTTT[i, 1] + TTT * ReachKMT_ind[i, 2] / 100;
    if (TTindT[i, 2] > 0) lg_pTTT[TTindT[i, 2]] = log(pTTT[i, 2]);

    PcapT[i, 2] = inv_logit(P_b[yrindT[i], 3]);
    PcapT[i, 3] = inv_logit(P_b[yrindT[i], 4]);
    chiT[i, 3] = 1.0;
    chiT[i, 2] = (1 - survT[i, 2]) + survT[i, 2] * (1 - PcapT[i, 3]) * chiT[i, 3];
    chiT[i, 1] = (1 - survT[i, 1]) + survT[i, 1] * (1 - PcapT[i, 2]) * chiT[i, 2];
  }

  // likelihood (and hierarchical terms) exactly as in CovIndCont.stan
  ObsTT ~ lognormal(lg_pTT, Pro_sd);
  for (j in 1:Nreaches) TT_RE[, j] ~ normal(0, TTRE_sd[j]);

  ObsTTT ~ lognormal(lg_pTTT, Pro_sdT);
  for (j in 1:NreachesT) TT_RET[, j] ~ normal(0, TTRE_sdT[j]);

  for (i in 1:Nind) {
    for (j in (firstCap[i] + 1):lastCap[i]) {
      1 ~ bernoulli(surv[i, j - 1]);
      CH[i, j] ~ bernoulli(Pcap[i, j]);
    }
    1 ~ bernoulli(chi[i, lastCap[i]]);
  }

  for (j in 1:Nreaches) S_RE[, j] ~ normal(0, RE_sd[j]);
  for (j in 1:Nreaches) P_b[1:Nyrs, j] ~ normal(muPb[j], sdPb[j]);

  for (i in 1:NindT) {
    for (j in (firstCapT[i] + 1):lastCapT[i]) {
      1 ~ bernoulli(survT[i, j - 1]);
      CHT[i, j] ~ bernoulli(PcapT[i, j]);
    }
    1 ~ bernoulli(chiT[i, lastCapT[i]]);
  }

  for (j in 1:NreachesT) S_REt[, j] ~ normal(0, RE_sdT[j]);
}
