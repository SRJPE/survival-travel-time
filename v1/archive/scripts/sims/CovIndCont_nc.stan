// CovIndCont_nc.stan: candidate revision of CovIndCont.stan (2026-09).
// Same likelihood as CovIndCont.stan (mainstem CJS + Butte/Feather CJS +
// lognormal cumulative travel time). Two changes:
//   1. Non-centered random effects: detection by year (P_b), release-group
//      survival (S_RE, S_REt) and travel time (TT_RE, TT_RET). The 03 fit
//      showed centered-parameterization trouble in the terminal reach
//      (TTRE_sd[4] Rhat 1.12, ESS 21; RE_sd[4] Rhat 1.03, ESS 82).
//   2. Explicit weakly informative priors. CovIndCont.stan has none, so its
//      SD parameters have improper flat priors.
// Parameter names match CovIndCont.stan so outputs are directly comparable.
// Per-fish quantities are model-block locals (not saved) to keep output small.
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
  array[Nyrs] vector[Nreaches] z_Pb;
  vector[Nreaches] muPb;
  vector<lower=0.001>[Nreaches] sdPb;
  vector<lower=-5, upper=5>[Nreaches] S_bReach;
  vector<lower=-5, upper=5>[Nreaches] T_bReach;
  real S_bCov;
  real TT_bCov;
  real S_bSz;
  real T_bSz;
  array[Ntribs] real S_bTrib;
  array[Ntribs] real T_bTrib;
  array[Ntribs] real S_bCovT;
  array[Ntribs] real TT_bCovT;

  matrix[Nrg, Nreaches] z_S;
  vector<lower=0>[Nreaches] RE_sd;
  matrix[Nrg, Nreaches] z_TT;
  vector<lower=0>[Nreaches] TTRE_sd;
  matrix[NrgT, NreachesT] z_St;
  vector<lower=0>[NreachesT] RE_sdT;
  matrix[NrgT, NreachesT] z_TTt;
  vector<lower=0>[NreachesT] TTRE_sdT;

  real<lower=0.001, upper=5> Pro_sd;
  real<lower=0.001, upper=5> Pro_sdT;
}

transformed parameters {
  array[Nyrs] vector[Nreaches] P_b;
  matrix[Nrg, Nreaches] S_RE = diag_post_multiply(z_S, RE_sd);
  matrix[Nrg, Nreaches] TT_RE = diag_post_multiply(z_TT, TTRE_sd);
  matrix[NrgT, NreachesT] S_REt = diag_post_multiply(z_St, RE_sdT);
  matrix[NrgT, NreachesT] TT_RET = diag_post_multiply(z_TTt, TTRE_sdT);
  for (y in 1:Nyrs) P_b[y] = muPb + sdPb .* z_Pb[y];
}

model {
  array[Nind] vector[Nreaches] surv;
  array[Nind] vector[Ndetlocs] Pcap;
  array[Nind] vector[Ndetlocs] chi;
  array[NindT] vector[NreachesT] survT;
  array[NindT] vector[3] PcapT;
  array[NindT] vector[3] chiT;
  vector[Nobs] lg_pTT;
  vector[NobsT] lg_pTTT;

  // Sacramento mainstem fish: travel time
  for (i in 1:Nind) {
    real cum = 0;
    for (j in 1:Nreaches) {
      cum += exp(T_bReach[j] + T_bSz * Sz[i] + TT_bCov * CovX[i, j] + TT_RE[rgind[i], j])
             * ReachKM_ind[i, j] / 100;
      if (TTind[i, j] > 0) lg_pTT[TTind[i, j]] = log(cum);
    }
  }

  // Sacramento mainstem fish: survival and detection
  for (i in 1:Nind) {
    for (j in 1:Ndetlocs) {
      if (j < Ndetlocs) {
        surv[i, j] = inv_logit(S_bReach[rch_covind[j]] + S_bCov * CovX[i, j]
                               + UseSizeEffect * S_bSz * Sz[i]
                               + S_RE[rgind[i], j]) ^ Rmult[i, j];
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
    real tt1;
    real tt2;
    survT[i, 1] = inv_logit(S_bTrib[trib_ind[i]] + S_bCovT[trib_ind[i]] * CovXT[i, 1]
                            + UseSizeEffect * S_bSz * SzT[i] + S_REt[rgindT[i], 1]) ^ RmultT[i];
    tt1 = exp(T_bTrib[trib_ind[i]] + T_bSz * SzT[i] + TT_bCovT[trib_ind[i]] * CovXT[i, 1]
              + TT_RET[rgindT[i], 1]) * ReachKMT_ind[i, 1] / 100;
    if (TTindT[i, 1] > 0) lg_pTTT[TTindT[i, 1]] = log(tt1);

    survT[i, 2] = inv_logit(S_bReach[Nreaches] + S_bCov * CovXT[i, NreachesT]
                            + UseSizeEffect * S_bSz * SzT[i]
                            + S_REt[rgindT[i], NreachesT]) ^ Rmult[1, Nreaches];
    tt2 = tt1 + exp(T_bReach[Nreaches] + T_bSz * SzT[i] + TT_bCov * CovXT[i, NreachesT]
                    + TT_RET[rgindT[i], NreachesT]) * ReachKMT_ind[i, 2] / 100;
    if (TTindT[i, 2] > 0) lg_pTTT[TTindT[i, 2]] = log(tt2);

    PcapT[i, 2] = inv_logit(P_b[yrindT[i], 3]);
    PcapT[i, 3] = inv_logit(P_b[yrindT[i], 4]);
    chiT[i, 3] = 1.0;
    chiT[i, 2] = (1 - survT[i, 2]) + survT[i, 2] * (1 - PcapT[i, 3]) * chiT[i, 3];
    chiT[i, 1] = (1 - survT[i, 1]) + survT[i, 1] * (1 - PcapT[i, 2]) * chiT[i, 2];
  }

  // priors
  muPb ~ normal(0, 2.5);
  sdPb ~ normal(0, 2);
  S_bReach ~ normal(0, 2.5);
  T_bReach ~ normal(0, 2.5);
  S_bTrib ~ normal(0, 2.5);
  T_bTrib ~ normal(0, 2.5);
  S_bCov ~ normal(0, 2);
  TT_bCov ~ normal(0, 2);
  S_bSz ~ normal(0, 2);
  T_bSz ~ normal(0, 2);
  S_bCovT ~ normal(0, 2);
  TT_bCovT ~ normal(0, 2);
  RE_sd ~ normal(0, 2);
  RE_sdT ~ normal(0, 2);
  TTRE_sd ~ normal(0, 1);
  TTRE_sdT ~ normal(0, 1);
  Pro_sd ~ normal(0, 1);
  Pro_sdT ~ normal(0, 1);

  for (y in 1:Nyrs) z_Pb[y] ~ std_normal();
  to_vector(z_S) ~ std_normal();
  to_vector(z_TT) ~ std_normal();
  to_vector(z_St) ~ std_normal();
  to_vector(z_TTt) ~ std_normal();

  // travel time
  ObsTT ~ lognormal(lg_pTT, Pro_sd);
  ObsTTT ~ lognormal(lg_pTTT, Pro_sdT);

  // CJS, mainstem
  for (i in 1:Nind) {
    for (j in (firstCap[i] + 1):lastCap[i]) {
      1 ~ bernoulli(surv[i, j - 1]);
      CH[i, j] ~ bernoulli(Pcap[i, j]);
    }
    1 ~ bernoulli(chi[i, lastCap[i]]);
  }

  // CJS, tributaries
  for (i in 1:NindT) {
    for (j in (firstCapT[i] + 1):lastCapT[i]) {
      1 ~ bernoulli(survT[i, j - 1]);
      CHT[i, j] ~ bernoulli(PcapT[i, j]);
    }
    1 ~ bernoulli(chiT[i, lastCapT[i]]);
  }
}
