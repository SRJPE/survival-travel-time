//CJS model taken from stan example at https://mc-stan.org/docs/stan-users-guide/mark-recapture-models.html
//This version uses year-reach specific pCaps. 
// Survival influenced by fish size at release or 2 or 3-level water year type fixed effect or both cumulatively
//and random effect for year-release group
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
}

parameters {
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] P_b;// reach -specific detection probability
  vector <lower=-10,upper=10> [Nreaches] muPb; // mean detection probability
  vector <lower=0.001> [Nreaches] sdPb; // standard deviation of detection probability
  vector <lower=-5,upper=5> [Nreaches] S_bReach; //survival by reach under baseline covariate
  vector <lower=-5,upper=5> [Nreaches] T_bReach; //survival by reach under baseline covariate
  real <lower=-10,upper=10>  S_bCov; //sac mainstem environmental covariate
  real <lower=-10,upper=10> TT_bCov; //sac mainstem environmental covariate for travel time
  array[Nrg, Nreaches] real <lower=-10, upper=10> S_RE; 
  vector <lower=0.001> [Nreaches] RE_sd; //sac mainstem SD of random effect
  real <lower=-10,upper=10> S_bSz; // fish size covariate
  real <lower=-10,upper=10> T_bSz; // fish size covariate
  array[Ntribs] real <lower=-10,upper=10> S_bTrib;
  array[Ntribs] real <lower=-10,upper=10> T_bTrib;
  array[NrgT, NreachesT] real  <lower=-10, upper=10> S_REt;
  vector <lower=0.001> [NreachesT] RE_sdT; //tribSD of random effect
  array[Ntribs] real <lower=-10,upper=10> S_bCovT;//trib environmental covariate
  array[Ntribs] real <lower=-10,upper=10> TT_bCovT;//trib environmental covariate
   
  real <lower=0.001, upper=5> Pro_sd;
  array[Nrg, Nreaches] real <lower=-10, upper=10> TT_RE;
  vector <lower=0.001> [Nreaches] TTRE_sd;
  real <lower=0.001, upper=5> Pro_sdT;
  array[NrgT, NreachesT] real  <lower=-10, upper=10> TT_RET;
 vector <lower=0.001> [NreachesT] TTRE_sdT;
}

transformed parameters{
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
}

model {
    ObsTT[]~lognormal(lg_pTT[],Pro_sd);
   for(j in 1:Nreaches){TT_RE[,j]~normal(0,TTRE_sd[j]);}

   ObsTTT[]~lognormal(lg_pTTT[],Pro_sdT);
    for(j in 1:NreachesT){ TT_RET[,j] ~normal(0,TTRE_sdT[j]);}
    
  for(i in 1:Nind){
    for(j in (firstCap[i]+1):lastCap[i]){  //Loop through all detection stations after first to last station individual was detected at
    1~bernoulli(surv[i,j-1]);//had to be alive prior to last detection station
    CH[i,j]~bernoulli(Pcap[i,j]);
    }
    1~bernoulli(chi[i,lastCap[i]]);//probability of individual never beeing seen again after last detection
  }
  

  for (j in 1:Nreaches){
      S_RE[,j] ~ normal(0,RE_sd[j]);
    }
  
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
  
  for (j in 1:NreachesT){
        S_REt[,j] ~ normal(0,RE_sdT[j]);
    }
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
  
}

