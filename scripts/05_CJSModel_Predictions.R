#### Generated quantities R code
#Clear the workspace
rm(list = ls())

library(rstan)
library(posterior)
library(tidyverse)
library(here)
library(loo)

source("scripts/02_GetData.R")

load(file="results/fit_CovIndCont_MaxFlow_FL.Rdata")
draws <- as_draws_array(fit)
draws_df <- as_draws_df(fit)

post <- rstan::extract(fit) 
nsim <- length(post$S_bCov) # number of posterior draws

# --- Extract posterior draws ---
P_b <- post$P_b
S_bReach <- post$S_bReach    
S_RE     <- post$S_RE         
T_bReach <- post$T_bReach
TT_RE    <- post$TT_RE
S_bCov   <- post$S_bCov
TT_bCov  <- post$TT_bCov
S_bSz    <- post$S_bSz 
T_bSz    <- post$T_bSz 
RE_sd    <- post$RE_sd
TTRE_sd  <- post$TTRE_sd
S_bTrib  <- post$S_bTrib
S_bCovT   <- post$S_bCovT
S_REt     <- post$S_REt   
T_bTrib  <- post$T_bTrib
TT_RET    <- post$TT_RET
TT_bCovT  <- post$TT_bCovT
RE_sdT    <- post$RE_sdT
TTRE_sdT  <- post$TTRE_sdT

# --- Define forecast variables ---
mux <- mean(MaxflowSac,na.rm=TRUE)
sdx <- sd(MaxflowSac,na.rm=TRUE)
CovX <- as.matrix(data.frame(cbind(MaxflowSac.z,MaxflowSac.z,MaxflowSac.z,MaxflowDelta.z)))
Xvec <- (seq(from =min(MaxflowSac,na.rm=TRUE),to=max(MaxflowSac,na.rm=TRUE),length.out=NsX)-mux)/sdx
musz <- mean(c(FL,FL_T),na.rm=TRUE)
sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
Xsz  <- (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
muxB <- mean(MaxflowB,na.rm=TRUE)
sdxB <- sd(MaxflowB,na.rm=TRUE)
muxF <- mean(MaxflowF,na.rm=TRUE)
sdxF <- sd(MaxflowF,na.rm=TRUE)
CovXT <-data.frame(cbind(MaxflowT.z,MaxflowDeltaT.z))
XvecT <- data.frame(cbind((seq(from =min(MaxflowB),to=max(MaxflowB),length.out=NsX)-muxB)/sdxB,
                          ((seq(from =min(MaxflowF),to=max(MaxflowF),length.out=NsX)-muxF)/sdxF)))

# ----- Detection Probability --------------
pred_pcap <- 1 / (1 + exp(-P_b))

# --- Sac fish survival and travel time predictions ------
TT_reach <- array(NA, dim=c(nsim, Nind, Nreaches))
TT_RelSac  <- matrix(NA, nrow=nsim,ncol=Nind)
pred_surv <- array(NA, dim=c(nsim, Nind, Nreaches))
pred_surv_per100  <- array(NA, dim=c(nsim, Nind, Nreaches))
SurvRelSac <- matrix(NA,nrow=nsim,ncol=Nind)
SurvRelSacSz <- array(NA, dim=c(nsim, Nsz, NsX))

for(i in 1:Nind){
  TT_reach[,i,1]= exp(T_bReach[,rch_covind[1]] + TT_bCov*CovX[i,1] + TT_RE[,rgind[i],1])*ReachKM_ind[i,1]/100
  
  TT_RelSac[,i]= exp(T_bReach[,rch_covind[1]] + TT_bCov*CovX[i,1] + TT_RE[,rgind[i],1])*ReachKM_ind[i,1]/100 +
    exp(T_bReach[,rch_covind[2]] + TT_bCov*CovX[i,2] + TT_RE[,rgind[i],2])*ReachKM_ind[i,2]/100 +
    exp(T_bReach[,rch_covind[3]] + TT_bCov*CovX[i,3] + TT_RE[,rgind[i],3])*ReachKM_ind[i,3]/100
  
  pred_surv[,i,1]=plogis(S_bReach[,rch_covind[1]] + S_bCov*CovX[i,1]+ S_RE[,rgind[i],1])^RmultSac
  
  pred_surv_per100[,i,1]=plogis(S_bReach[,rch_covind[1]] +  S_bCov*CovX[i,1] + S_RE[,rgind[i],1])
  
  for(j in 2:Nreaches){ 
    pred_surv[,i,j]=plogis(S_bReach[,rch_covind[j]] +  S_bCov*CovX[i,j] + S_RE[,rgind[i],j])^Rmult[1,j]
    
    TT_reach[,i,j]= exp(T_bReach[,rch_covind[j]] + TT_bCov*CovX[i,j] + TT_RE[,rgind[i],j])*ReachKM_ind[i,j]/100
    
    pred_surv_per100[,i,j]=plogis(S_bReach[,rch_covind[j]] +  S_bCov*CovX[i,j] + S_RE[,rgind[i],j])
  }
  
  SurvRelSac[,i] = pred_surv[,i,1]*pred_surv[,i,2]*pred_surv[,i,3]
}


for(ix in 1:Nsz){
  for(j in 1:NsX){
    SurvRelSacSz[,ix,j]= plogis(S_bReach[,1] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,1] *  
      plogis(S_bReach[,2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,2] * 
      plogis(S_bReach[,3] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,3]
  }
}


#------- Survival and travel time forecasts for Sac fish  -------
SurvForecast <- matrix(NA, nrow=nsim, ncol = NsX)
SurvForecast_nore <-matrix(NA, nrow=nsim, ncol = NsX)
TTForecast <- matrix(NA, nrow=nsim, ncol = NsX)

SurvForecastSz <- array(NA, dim=c(nsim, Nsz, NsX))
SurvForecastSz_nore <- array(NA, dim=c(nsim, Nsz, NsX))
TTForecastSz   <- array(NA, dim=c(nsim, Nsz, NsX))

SurvForecastSz_rst <- array(NA, dim=c(nsim, Nsz, NsX, Nrst))
TTForecastSz_rst   <- array(NA, dim=c(nsim, Nsz, NsX, Nrst))

SurvForecastSz_rst_kni <- array(NA, dim=c(nsim, Nsz, NsX, Nrst))
TTForecastSz_rst_kni <- array(NA, dim=c(nsim, Nsz, NsX, Nrst))
SurvForecastSz_rst_tis <- array(NA, dim=c(nsim, Nsz, NsX, Nrst))
TTForecastSz_rst_tis <- array(NA, dim=c(nsim, Nsz, NsX, Nrst))

s_re1 <- rnorm(nsim, 0, RE_sd[,1])
s_re2 <- rnorm(nsim, 0, RE_sd[,2])
s_re3 <- rnorm(nsim, 0, RE_sd[,3])
tt_re1 <- rnorm(nsim, 0, TTRE_sd[,1])
tt_re2 <- rnorm(nsim, 0, TTRE_sd[,2])
tt_re3 <- rnorm(nsim, 0, TTRE_sd[,3])

for(j in 1:NsX){
  SurvForecast[,j]= plogis(S_bReach[,1] + S_bCov*Xvec[j] + s_re1)^RmultSac * 
    plogis(S_bReach[,2] + S_bCov*Xvec[j] + s_re2)^Rmult[1,2] *
    plogis(S_bReach[,3] + S_bCov*Xvec[j] + s_re3)^Rmult[1,3]
  
  SurvForecast_nore[,j]= plogis(S_bReach[,1] + S_bCov*Xvec[j])^RmultSac* 
    plogis(S_bReach[,2] + S_bCov*Xvec[j])^Rmult[1,2] *
    plogis(S_bReach[,3] + S_bCov*Xvec[j])^Rmult[1,3]
  
  TTForecast[,j] = exp(T_bReach[,1] + TT_bCov *Xvec[j] + tt_re1)*ReachKM[1]/100 +
    exp(T_bReach[,2] + TT_bCov *Xvec[j] + tt_re2)*ReachKM[2]/100+
    exp(T_bReach[,3] + TT_bCov *Xvec[j] + tt_re3)*ReachKM[3]/100
}

for(ix in 1:Nsz){
  for(j in 1:NsX){
    SurvForecastSz[,ix,j]= plogis(S_bReach[,1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re1)^RmultSac * 
      plogis(S_bReach[,2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
      plogis(S_bReach[,3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re3)^Rmult[1,3]
    
    SurvForecastSz_nore[,ix,j]= plogis(S_bReach[,1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix])^RmultSac * 
      plogis(S_bReach[,2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,2] *
      plogis(S_bReach[,3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,3]
    
    TTForecastSz[,ix,j]= exp(T_bReach[,1]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re1)*ReachKM[1]/100 +
      exp(T_bReach[,2]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re2)*ReachKM[2]/100 +
      exp(T_bReach[,3]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re3)*ReachKM[3]/100
    
    for(i in 1:Nrst){
      SurvForecastSz_rst[,ix,j,i]= plogis(S_bReach[,1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re1)^Rmultrst[i] *
        plogis(S_bReach[,2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
        plogis(S_bReach[,3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re3)^Rmult[1,3]
      
      TTForecastSz_rst[,ix,j,i]= exp(T_bReach[,1]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re1)*ReachKMrst[i]/100 +
        exp(T_bReach[,2]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re2)*ReachKM[2]/100 +
        exp(T_bReach[,3]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re3)*ReachKM[3]/100
      
      SurvForecastSz_rst_tis[,ix,j,i]= plogis(S_bReach[,1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re1)^Rmultrst[i] *
        plogis(S_bReach[,2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
        plogis(S_bReach[,3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re3)^Rmult_Tis
      
      TTForecastSz_rst_tis[,ix,j,i]= exp(T_bReach[,1]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re1)*ReachKMrst[i]/100 +
        exp(T_bReach[,2]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re2)*ReachKM[2]/100 +
        exp(T_bReach[,3]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re3)*ReachKM_Tis/100
      
      SurvForecastSz_rst_kni[,ix,j,i]= plogis(S_bReach[,1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re1)^Rmultrst[i] *
        plogis(S_bReach[,2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
        plogis(S_bReach[,3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re3)^Rmult_Kni
      
      TTForecastSz_rst_kni[,ix,j,i]= exp(T_bReach[,1]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re1)*ReachKMrst[i]/100 +
        exp(T_bReach[,2]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re2)*ReachKM[2]/100 +
        exp(T_bReach[,3]+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re3)*ReachKM_Kni/100
      
    }
  }
}


# --- Trib fish survival and travel time predictions ------
TT_reachT <- array(NA, dim=c(nsim, NindT, Ntribs))
pred_survT <- array(NA, dim=c(nsim, NindT, Ntribs))
pred_survT_per100  <- array(NA, dim=c(nsim, NindT, Ntribs))
pred_survTSz <- array(NA, dim=c(nsim,Nsz,NsX,Ntribs))

for(i in 1:NindT){
  pred_survT[,i,1]=plogis(S_bTrib[,trib_ind[i]] + S_bCovT[,trib_ind[i]]*CovXT[i,1] + S_REt[,rgindT[i],1])^RmultTrib[trib_ind[i]]
  pred_survT[,i,2]=plogis(S_bReach[,Nreaches]  + S_bCov*CovXT[i,2] + S_REt[,rgindT[i],2] )^Rmult[1,Nreaches]
  
  pred_survT_per100[,i,1]=plogis(S_bTrib[,trib_ind[i]] + S_bCovT[,trib_ind[i]]*CovXT[i,1] + S_REt[,rgindT[i],1])
  pred_survT_per100[,i,2]=plogis(S_bReach[,Nreaches] + S_bCov*CovXT[i,2] + S_REt[,rgindT[i],2] )
  
  TT_reachT[,i,1]= exp(T_bTrib[,trib_ind[i]] + TT_bCovT[,trib_ind[i]] *CovXT[i,1] +TT_RET[,rgindT[i],1])*ReachKMT_ind[i,1]/100
  TT_reachT[,i,2]= exp(T_bReach[,Nreaches]  + TT_bCov *CovXT[i,2] +TT_RET[,rgindT[i],2])*ReachKMT_ind[i,2]/100
}

for(ix in 1:Nsz){
  for(itrib in 1:Ntribs){
    for(j in 1:NsX){
      pred_survTSz[,ix,j,itrib]=plogis(S_bTrib[,itrib] + S_bCovT[,itrib]*XvecT[j,itrib] + S_bSz*Xsz[ix])^RmultTrib[itrib]
    }
  }
}

#------- Survival and travel time forecasts for Trib fish  -------
TribSurvForecast <- array(NA, dim=c(nsim,NsX,Ntribs))
TribSurvForecast_nore <-  array(NA, dim=c(nsim,NsX,Ntribs))
TTForecastT <- array(NA, dim=c(nsim,NsX,Ntribs))

TribSurvForecastSz <- array(NA, dim=c(nsim, Nsz, NsX,Ntribs))
TribSurvForecastSz_nore <- array(NA, dim=c(nsim, Nsz, NsX,Ntribs))
TTForecastTSz   <- array(NA, dim=c(nsim, Nsz, NsX,Ntribs))

TribSurvForecastSz_rst <- array(NA, dim=c(nsim, Nsz, NsX, NrstT))
TTForecastTSz_rst   <- array(NA, dim=c(nsim, Nsz, NsX, NrstT))

s_ret <- rnorm(nsim, 0, RE_sdT[,1])
tt_ret <- rnorm(nsim, 0, TTRE_sdT[,1])

for(itrib in 1:Ntribs){
  for(j in 1:NsX){
    TribSurvForecast[,j,itrib]= plogis(S_bTrib[,itrib] + S_bCovT[,itrib]*XvecT[j,itrib] + s_ret)^RmultTrib[itrib]
    
    TribSurvForecast_nore[,j,itrib]= plogis(S_bTrib[,itrib] + S_bCovT[,itrib]*XvecT[j,itrib] )^RmultTrib[itrib]
    
    TTForecastT[,j,itrib] = exp(T_bTrib[,itrib] + TT_bCovT[,itrib] *XvecT[j,itrib] + tt_ret)*ReachKMT[itrib,1]/100
  }
}

for(ix in 1:Nsz){
  for(j in 1:NsX){
    for(itrib in 1:Ntribs){
      TribSurvForecastSz[,ix,j,itrib]=plogis(S_bTrib[,itrib] + S_bCovT[,itrib]*XvecT[j,itrib] + S_bSz*Xsz[ix]+ s_ret)^RmultTrib[itrib]
      
      TribSurvForecastSz_nore[,ix,j,itrib]=plogis(S_bTrib[,itrib] + S_bCovT[,itrib]*XvecT[j,itrib] + S_bSz*Xsz[ix])^RmultTrib[itrib]
      
      TTForecastTSz[,ix,j,itrib] = exp(T_bTrib[,itrib] + TT_bCovT[,itrib] *XvecT[j,itrib] + T_bSz*Xsz[ix]+ tt_ret)*ReachKMT[itrib,1]/100
    }
    for(i in 1:NrstT){
      TTForecastTSz_rst[,ix,j,i] = exp(T_bTrib[,TribForRST[i]] + TT_bCovT[,TribForRST[i]] *XvecT[j,TribForRST[i]] + T_bSz*Xsz[ix]+ tt_ret)*ReachKMTrst[i]/100
      
      TribSurvForecastSz_rst[,ix,j,i]=plogis(S_bTrib[,TribForRST[i]] + S_bCovT[,TribForRST[i]]*XvecT[j,TribForRST[i]] + S_bSz*Xsz[ix]+ s_ret)^RmultTrst[i]
    }
  }
}



vars_to_save <- c("P_b","S_bReach","S_RE","T_bReach","TT_RE","S_bSz","S_bCov","TT_bCov","T_bSz","RE_sd","TTRE_sd",
                  "S_REt","TT_RET", "TT_bCovT","S_bTrib","S_bCovT","RE_sdT","T_bTrib","TT_bCovT","T_bSz","TTRE_sdT",
                  "mux","sdx","CovX","Xvec","musz","sdsz","Xsz", "muxB","sdxB","muxF","sdxF","CovXT","XvecT",
                  "pred_pcap","pred_surv","pred_surv_per100","TT_reach", "TT_RelSac", "SurvRelSac", "SurvRelSacSz",
                  "pred_survT","pred_survT_per100", "TT_reachT", "pred_survTSz",
                  "SurvForecast","SurvForecast_nore","TTForecast",
                  "SurvForecastSz","TTForecastSz","SurvForecastSz_nore","SurvForecastSz_rst","TTForecastSz_rst",
                  "SurvForecastSz_rst_tis","TTForecastSz_rst_tis","SurvForecastSz_rst_kni","TTForecastSz_rst_kni",
                  "TribSurvForecast","TribSurvForecast_nore","TTForecastT", 
                  "TribSurvForecastSz","TribSurvForecastSz_nore","TTForecastTSz",
                  "TribSurvForecastSz_rst","TTForecastTSz_rst")  

stm <-mget(vars_to_save)

save(stm, file = here("results", paste0("stm_", output_tag, ".Rdata")))

# LOOIC calculation ----------------------------------------------------------
log_lik <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
(loo <- loo(log_lik, r_eff = r_eff))
(pointwise_elpd <- loo$pointwise[ , "elpd_loo"])
loopointwise_stats <- loo$pointwise 
write.csv(as.data.frame(loopointwise_stats), tagged_table_file("loo_pointwise"), row.names = FALSE)
write.csv(as.data.frame(loo$estimates), tagged_table_file("loo_summary"))

# Detection probability figure ---------------------------------------------------------------------------
years <- 2013:2024
locations <- c("Woodson", "Butte", "Sacramento", "Delta")

pred <- stm$pred_pcap

df <- expand.grid(
  iter = 1:dim(pred)[1],
  year = years,
  location = locations
)

df$value <- as.vector(pred)

summary_df <- df %>%
  group_by(year, location) %>%
  summarise(
    mean = mean(value),
    lwr = quantile(value, 0.025),
    upr = quantile(value, 0.975),
    .groups = "drop"
  )

(pcap_fig <- ggplot(summary_df,
                                   aes(x = location, y = mean, colour = location)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = location), width = 0.2) +
    facet_wrap(~year, nrow = 1) +
    ylab("Detection Probability") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1))+
    scale_colour_viridis_d(direction = 1,name = "location")
)

ggsave(plot=pcap_fig, filename = tagged_fig_file("pcap"), width =10, height = 4)

# survival figures --------------------------------------------------------------------------------
# plot upper sac reach specific survival for each release group
reach_names <-c('Rel-Wood','Wood-But','But-Sac','Sac-Delta')

pred <- stm$pred_surv_per100

df <- expand.grid(
  iter  = 1:dim(pred)[1],
  fish  = 1:dim(pred)[2],
  reach = reach_names
)

df$survival <- as.vector(pred)

study_names <- d_Sac_sort %>%
  select(StudyID,year) %>%
  distinct()

study_lookup <- data.frame(
  study_num = 1:dim(study_names)[1],
  study = study_names$StudyID,
  year = study_names$year)

fish_study <- d_Sac_sort$StudyID
fish_meta <- data.frame(
  fish = 1:dim(pred)[2],
  study_id = fish_study
)

fish_meta <- left_join(fish_meta, study_lookup, by=c("study_id"="study"))

df_surv <- left_join(df, fish_meta, by="fish")

summary_df_surv <- df_surv %>%
  group_by(study_id, year, reach) %>%
  summarise(
    mean = mean(survival),
    lwr  = quantile(survival, 0.025),
    upr  = quantile(survival, 0.975),
    .groups = "drop"
  )

(surv_perreach100_fig <- ggplot(summary_df_surv,
                                           aes(x = reach, y = mean, colour = reach)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = reach), width = 0.2) +
    facet_wrap(~ year + study_id, scales = "free_x") +
    ylab("Survival rate") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Survival rate")+
    scale_colour_viridis_d(direction = 1,name = "reach")
)

ggsave(plot=surv_perreach100_fig, filename = tagged_fig_file("surv_perreach100"), width =13, height = 10)

# plot trib reach specific survival for each release group
reach_namesT <-c('Rel-Sac','Sac-Delta')

predT <- stm$pred_survT_per100

dfT <- expand.grid(
  iter  = 1:dim(predT)[1],
  fish  = 1:dim(predT)[2],
  reach = reach_namesT
)

dfT$survival <- as.vector(predT)

study_namesT <- d_FeaBut_sort %>%
  select(StudyID,year) %>%
  distinct()

study_lookupT <- data.frame(
  study_num = 1:dim(study_namesT)[1],
  study = study_namesT$StudyID,
  year = study_namesT$year)

fish_studyT <- d_FeaBut_sort$StudyID
fish_metaT <- data.frame(
  fish = 1:dim(predT)[2],
  study_id = fish_studyT)

fish_metaT <- left_join(fish_metaT, study_lookupT, by=c("study_id"="study"))

df_survT <- left_join(dfT, fish_metaT, by="fish")

summary_df_survT <- df_survT %>%
  group_by(study_id, year, reach) %>%
  summarise(
    mean = mean(survival),
    lwr  = quantile(survival, 0.025),
    upr  = quantile(survival, 0.975),
    .groups = "drop"
  )

(surv_tribperreach100_fig <- ggplot(summary_df_survT,
                                                   aes(x = reach, y = mean, colour = reach)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = reach), width = 0.2) +
    facet_wrap(~ year + study_id, scales = "free_x") +
    ylab("Survival rate") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Survival rate")+
    scale_colour_viridis_d(direction = 1,name = "reach")
)

ggsave(plot=surv_tribperreach100_fig, filename = tagged_fig_file("surv_tribperreach100"), width =15, height = 12)

# plot survival from upper sac release to Sac for each release group
pred_relsac <- stm$SurvRelSac

df_relsac <- expand.grid(
  iter = 1:dim(pred_relsac)[1],
  fish = 1:dim(pred_relsac)[2]
)

df_relsac$survival <- as.vector(pred_relsac)

study_names <- d_Sac_sort %>%
  select(StudyID,year) %>%
  distinct()

study_lookup <- data.frame(
  study_num = 1:dim(study_names)[1],
  study = study_names$StudyID,
  year = study_names$year)

fish_study <- d_Sac_sort$StudyID
fish_meta <- data.frame(
  fish = 1:length(fish_study),
  study_id = fish_study
)

fish_meta <- left_join(fish_meta, study_lookup, by=c("study_id"="study"))

df_surv_relsac <- left_join(df_relsac, fish_meta, by="fish")

summary_df_surv_relsac <- df_surv_relsac %>%
  group_by(study_id,year) %>%
  summarise(
    mean = mean(survival),
    lwr  = quantile(survival, 0.025),
    upr  = quantile(survival, 0.975),
    .groups = "drop") %>% 
  mutate(WY = case_when(year %in% c(2013,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W')) %>% 
  arrange(year,study_id) %>% 
  mutate(study_id = factor(study_id, levels = unique(study_id)))

(surv_relsac_fig <- ggplot(summary_df_surv_relsac,
                                          aes(x = study_id, y = mean, colour = WY)) +
    geom_point(size = 4) +
    geom_errorbar(
      aes(ymin = lwr, ymax = upr, colour = WY),
      width = 0.2,
      size = 1 ) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 45, hjust = 1))+
    xlab("")+ylab("Rel - Sac survival rate")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=surv_relsac_fig, filename = tagged_fig_file("surv_relsac"), width =7, height = 6)

# plot survival from trib release to Sac for each release group
pred_relsacT <- stm$pred_survT
reach_namesT <-c('Rel-Sac','Sac-Delta')

df_relsacT <- expand.grid(
  iter  = 1:dim(pred_relsacT)[1],
  fish  = 1:dim(pred_relsacT)[2],
  reach = reach_namesT)

df_relsacT$survival <- as.vector(pred_relsacT)

study_namesT <- d_FeaBut_sort %>%
  select(StudyID,year) %>%
  distinct()

study_lookupT <- data.frame(
  study_num = 1:dim(study_namesT)[1],
  study = study_namesT$StudyID,
  year = study_namesT$year)

fish_studyT <- d_FeaBut_sort$StudyID
fish_metaT <- data.frame(
  fish = 1:dim(predT)[2],
  study_id = fish_studyT)

fish_metaT <- left_join(fish_metaT, study_lookupT, by=c("study_id"="study"))

df_surv_relsacT <- left_join(df_relsacT, fish_metaT, by="fish")

summary_df_surv_relsacT <- df_surv_relsacT %>% filter(reach == 'Rel-Sac') %>% 
  group_by(study_id, year) %>%
  summarise(mean = mean(survival),
            lwr  = quantile(survival, 0.025),
            upr  = quantile(survival, 0.975),
            .groups = "drop") %>% 
  mutate(WY = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W')) %>% 
  arrange(year,study_id) %>% 
  mutate(study_id = factor(study_id, levels = unique(study_id)))

(surv_tribrelsac_fig <- ggplot(summary_df_surv_relsacT,
                                              aes(x = study_id, y = mean, colour = WY)) +
    geom_point(size = 4) +
    geom_errorbar(
      aes(ymin = lwr, ymax = upr, colour = WY),
      width = 0.2,
      size = 1 ) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 45, hjust = 1))+
    xlab("")+ylab("Rel - Sac survival rate")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=surv_tribrelsac_fig, filename = tagged_fig_file("surv_tribrelsac"), width =7, height = 6)


# Model forecast figures ----------------------------------------------------------------------
# Define index ranges
i_vals <- 1:25   #flow range
j_vals <- 1:2   # trib locations
n_cov <- 25

########### Flow forecast
## Mainstem plot
# flow values
mux <- mean(MaxflowSac,na.rm=TRUE)
sdx <- sd(MaxflowSac,na.rm=TRUE)
Xvec <- seq(from =4000,to=40000,length.out=NsX)
flow_vec <- Xvec

survfor <- stm$SurvForecast

survfor_summary <- data.frame(
  Flow = flow_vec,
  mean  = apply(survfor, 2, mean),
  lwr   = apply(survfor, 2, quantile, 0.025),
  upr   = apply(survfor, 2, quantile, 0.975))

survfor_nore <- stm$SurvForecast_nore

survfor_nore_summary <- data.frame(
  Flow = flow_vec,
  mean  = apply(survfor_nore, 2, mean),
  lwr   = apply(survfor_nore, 2, quantile, 0.025),
  upr   = apply(survfor_nore, 2, quantile, 0.975))

(survforecast_Sac_fig <- ggplot() + 
    geom_ribbon(data = survfor_summary,aes(x = Flow, ymin = lwr, ymax = upr),fill = "grey70", alpha = 0.6) +
    geom_line(data = survfor_summary,aes(x = Flow, y = mean),color = "black", size = 1.2) +
    geom_rug(data=d_Sac_sort,aes(x = Maxflowsac), sides = "b", alpha = 0.5) + 
    geom_line(data=survfor_nore_summary,aes(x = Flow, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=survfor_nore_summary,aes(x =Flow, y = mean,ymin =lwr, ymax =upr),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Flow",
          y = "Release - Sacramento survival rate")+
    ggtitle('Sacramento')+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
  
) 

ggsave(plot=survforecast_Sac_fig, filename = tagged_fig_file("survforecast", "Sac"), width =10, height = 6)

##  Butte plot
# flow values
muxB <- mean(MaxflowB,na.rm=TRUE)
sdxB <- sd(MaxflowB,na.rm=TRUE)
muxF <- mean(MaxflowF,na.rm=TRUE)
sdxF <- sd(MaxflowF,na.rm=TRUE)

XvecB <- seq(from =min(MaxflowB),to=max(MaxflowB),length.out=NsX)
XvecF <- seq(from =min(MaxflowF),to=max(MaxflowF),length.out=NsX)

survforB <- stm$TribSurvForecast[,,1]

survforB_summary <- data.frame(
  Flow = XvecB,
  mean  = apply(survforB, 2, mean),
  lwr   = apply(survforB, 2, quantile, 0.025),
  upr   = apply(survforB, 2, quantile, 0.975))

survforB_nore <- stm$TribSurvForecast_nore[,,1]

survforB_nore_summary <- data.frame(
  Flow = XvecB,
  mean  = apply(survforB_nore, 2, mean),
  lwr   = apply(survforB_nore, 2, quantile, 0.025),
  upr   = apply(survforB_nore, 2, quantile, 0.975))

(survforecast_But_fig <- ggplot() + 
    geom_ribbon(data = survforB_summary,aes(x = Flow, ymin = lwr, ymax = upr),fill = "grey70", alpha = 0.6) +
    geom_line(data = survforB_summary,aes(x = Flow, y = mean),color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort,aes(x = MaxflowB), sides = "b", alpha = 0.5) + 
    geom_line(data=survforB_nore_summary,aes(x = Flow, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=survforB_nore_summary,aes(x =Flow, y = mean,ymin =lwr, ymax =upr),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Flow",
          y = "Release - Sacramento survival rate")+
    ggtitle('Butte')+
    scale_y_continuous(limits = c(0, 1))
  #  scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
  
) 

ggsave(plot=survforecast_But_fig, filename = tagged_fig_file("survforecast", "But"), width =10, height = 6)

## Feather River
survforF <- stm$TribSurvForecast[,,2]

survforF_summary <- data.frame(
  Flow = XvecF,
  mean  = apply(survforF, 2, mean),
  lwr   = apply(survforF, 2, quantile, 0.025),
  upr   = apply(survforF, 2, quantile, 0.975))

survforF_nore <- stm$TribSurvForecast_nore[,,2]

survforF_nore_summary <- data.frame(
  Flow = XvecF,
  mean  = apply(survforF_nore, 2, mean),
  lwr   = apply(survforF_nore, 2, quantile, 0.025),
  upr   = apply(survforF_nore, 2, quantile, 0.975))

(survforecast_Fea_fig <- ggplot() + 
    geom_ribbon(data = survforF_summary,aes(x = Flow, ymin = lwr, ymax = upr),fill = "grey70", alpha = 0.6) +
    geom_line(data = survforF_summary,aes(x = Flow, y = mean),color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort,aes(x = MaxflowF), sides = "b", alpha = 0.5) + 
    geom_line(data=survforF_nore_summary,aes(x = Flow, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=survforF_nore_summary,aes(x =Flow, y = mean,ymin =lwr, ymax =upr),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Flow",
          y = "Release - Sacramento survival rate")+
    ggtitle('Feather')+
    scale_y_continuous(limits = c(0, 1))
  #  scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
  
) 

ggsave(plot=survforecast_Fea_fig, filename = tagged_fig_file("survforecast", "Fea"), width =10, height = 6)

##### Size forecast
k_vals <- 1:25  # sizes range
size_vec <- seq(from=10,to=150,length.out=Nsz)

## Sacramento plot
survforSz <- stm$SurvForecastSz

survforSz_summary <- data.frame(
  Size = size_vec,
  mean  = apply(survforSz, 2, mean),
  lwr   = apply(survforSz, 2, quantile, 0.025),
  upr   = apply(survforSz, 2, quantile, 0.975))

survforSz_nore <- stm$SurvForecastSz_nore

survforSz_nore_summary <- data.frame(
  Size = size_vec,
  mean  = apply(survforSz_nore, 2, mean),
  lwr   = apply(survforSz_nore, 2, quantile, 0.025),
  upr   = apply(survforSz_nore, 2, quantile, 0.975))

(survforecastSz_Sac_fig <- ggplot() + 
    geom_ribbon(data = survforSz_summary,aes(x = Size, ymin = lwr, ymax = upr),fill = "grey70", alpha = 0.6) +
    geom_line(data = survforSz_summary,aes(x = Size, y = mean),color = "black", size = 1.2) +
    geom_rug(data=d_Sac_sort,aes(x = fish_length), sides = "b", alpha = 0.5) + 
    geom_line(data=survforSz_nore_summary,aes(x = Size, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=survforSz_nore_summary,aes(x = Size, y = mean,ymin =lwr, ymax =upr),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Size",
          y = "Release - Sacramento survival rate")+
    ggtitle('Sacramento')+
    scale_y_continuous(limits = c(0, 1))
  #    scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
  
) 

ggsave(plot=survforecastSz_Sac_fig, filename = tagged_fig_file("survforecastSz", "Sac"), width =10, height = 6)

## Butte plot
survforSz_B <- stm$TribSurvForecastSz[,,12,1]

survforSz_B_summary <- data.frame(
  Size = size_vec,
  mean  = apply(survforSz_B, 2, mean),
  lwr   = apply(survforSz_B, 2, quantile, 0.025),
  upr   = apply(survforSz_B, 2, quantile, 0.975))

survforSz_B_nore <- stm$TribSurvForecastSz_nore[,,12,1]

survforSz_B_nore_summary <- data.frame(
  Size = size_vec,
  mean  = apply(survforSz_B_nore, 2, mean),
  lwr   = apply(survforSz_B_nore, 2, quantile, 0.025),
  upr   = apply(survforSz_B_nore, 2, quantile, 0.975))

(survforecastSz_But_fig <- ggplot() + 
    geom_ribbon(data = survforSz_B_summary,aes(x = Size, ymin = lwr, ymax = upr),fill = "grey70", alpha = 0.6) +
    geom_line(data = survforSz_B_summary,aes(x = Size, y = mean),color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort %>% filter(!rl %in% c("1F","2F")),aes(x = fish_length), sides = "b", alpha = 0.5) + 
    geom_line(data=survforSz_B_nore_summary,aes(x = Size, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=survforSz_B_nore_summary,aes(x = Size, y = mean,ymin =lwr, ymax =upr),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Size",
          y = "Release - Sacramento survival rate")+
    ggtitle('Butte')+
    scale_y_continuous(limits = c(0, 1))
  #    scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
  
) 

ggsave(plot=survforecastSz_But_fig, filename = tagged_fig_file("survforecastSz", "But"), width =10, height = 6)


## Feather plot
survforSz_F <- stm$TribSurvForecastSz[,,12,2]

survforSz_F_summary <- data.frame(
  Size = size_vec,
  mean  = apply(survforSz_F, 2, mean),
  lwr   = apply(survforSz_F, 2, quantile, 0.025),
  upr   = apply(survforSz_F, 2, quantile, 0.975))

survforSz_F_nore <- stm$TribSurvForecastSz_nore[,,12,2]

survforSz_F_nore_summary <- data.frame(
  Size = size_vec,
  mean  = apply(survforSz_F_nore, 2, mean),
  lwr   = apply(survforSz_F_nore, 2, quantile, 0.025),
  upr   = apply(survforSz_F_nore, 2, quantile, 0.975))

(survforecastSz_Fea_fig <- ggplot() + 
    geom_ribbon(data = survforSz_F_summary,aes(x = Size, ymin = lwr, ymax = upr),fill = "grey70", alpha = 0.6) +
    geom_line(data = survforSz_F_summary,aes(x = Size, y = mean),color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort  %>% filter(rl %in% c("1F","2F")),aes(x = fish_length), sides = "b", alpha = 0.5) + 
    geom_line(data=survforSz_F_nore_summary,aes(x = Size, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=survforSz_F_nore_summary,aes(x = Size, y = mean,ymin =lwr, ymax =upr),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Size",
          y = "Release - Sacramento survival rate")+
    ggtitle('Feather')+
    scale_y_continuous(limits = c(0, 1))
  #    scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
  
) 

ggsave(plot=survforecastSz_Fea_fig, filename = tagged_fig_file("survforecastSz", "Fea"), width =10, height = 6)

# travel time figures -------------------------------------------------------------------------------
summarise_draw_matrix <- function(draw_matrix, x, x_name) {
  data.frame(
    x = x,
    mean = apply(draw_matrix, 2, mean, na.rm = TRUE),
    lwr = apply(draw_matrix, 2, quantile, 0.025, na.rm = TRUE),
    upr = apply(draw_matrix, 2, quantile, 0.975, na.rm = TRUE)
  ) %>% 
    rename(!!x_name := x)
}

# plot upper sac reach specific travel time for each release group
tt_reach_names <- c('Rel-Wood','Wood-But','But-Sac','Sac-Delta')
tt_pred <- stm$TT_reach

tt_df <- expand.grid(
  iter = 1:dim(tt_pred)[1],
  fish = 1:dim(tt_pred)[2],
  reach = tt_reach_names
)

tt_df$travel_time <- as.vector(tt_pred)

tt_fish_meta <- data.frame(
  fish = 1:dim(tt_pred)[2],
  study_id = d_Sac_sort$StudyID
) %>% 
  left_join(study_lookup, by = c("study_id" = "study"))

tt_summary_df <- tt_df %>% 
  left_join(tt_fish_meta, by = "fish") %>% 
  group_by(study_id, year, reach) %>% 
  summarise(
    mean = mean(travel_time, na.rm = TRUE),
    lwr = quantile(travel_time, 0.025, na.rm = TRUE),
    upr = quantile(travel_time, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

(tt_perreach_fig <- ggplot(tt_summary_df,
                                      aes(x = reach, y = mean, colour = reach)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = reach), width = 0.2) +
    facet_wrap(~ year + study_id, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Travel time (days)")+
    scale_colour_viridis_d(direction = 1,name = "reach")
)

ggsave(plot=tt_perreach_fig, filename = tagged_fig_file("tt_perreach"), width =13, height = 10)

# plot tributary reach specific travel time for each release group
tt_reach_namesT <- c('Rel-Sac','Sac-Delta')
tt_predT <- stm$TT_reachT

tt_dfT <- expand.grid(
  iter = 1:dim(tt_predT)[1],
  fish = 1:dim(tt_predT)[2],
  reach = tt_reach_namesT
)

tt_dfT$travel_time <- as.vector(tt_predT)

tt_fish_metaT <- data.frame(
  fish = 1:dim(tt_predT)[2],
  study_id = d_FeaBut_sort$StudyID
) %>% 
  left_join(study_lookupT, by = c("study_id" = "study"))

tt_summary_dfT <- tt_dfT %>% 
  left_join(tt_fish_metaT, by = "fish") %>% 
  group_by(study_id, year, reach) %>% 
  summarise(
    mean = mean(travel_time, na.rm = TRUE),
    lwr = quantile(travel_time, 0.025, na.rm = TRUE),
    upr = quantile(travel_time, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

(tt_tribperreach_fig <- ggplot(tt_summary_dfT,
                                          aes(x = reach, y = mean, colour = reach)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = reach), width = 0.2) +
    facet_wrap(~ year + study_id, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Travel time (days)")+
    scale_colour_viridis_d(direction = 1,name = "reach")
)

ggsave(plot=tt_tribperreach_fig, filename = tagged_fig_file("tt_tribperreach"), width =15, height = 12)

# plot travel time from upper sac release to Sacramento for each release group
tt_relsac <- stm$TT_RelSac

tt_relsac_df <- expand.grid(
  iter = 1:dim(tt_relsac)[1],
  fish = 1:dim(tt_relsac)[2]
)

tt_relsac_df$travel_time <- as.vector(tt_relsac)

tt_relsac_summary <- tt_relsac_df %>% 
  left_join(fish_meta, by = "fish") %>% 
  group_by(study_id, year) %>% 
  summarise(
    mean = mean(travel_time, na.rm = TRUE),
    lwr = quantile(travel_time, 0.025, na.rm = TRUE),
    upr = quantile(travel_time, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(WY = case_when(year %in% c(2013,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W')) %>% 
  arrange(year,study_id) %>% 
  mutate(study_id = factor(study_id, levels = unique(study_id)))

(tt_relsac_fig <- ggplot(tt_relsac_summary,
                                    aes(x = study_id, y = mean, colour = WY)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = WY), width = 0.2, size = 1) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 45, hjust = 1))+
    xlab("")+ylab("Rel - Sac travel time (days)")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=tt_relsac_fig, filename = tagged_fig_file("tt_relsac"), width =7, height = 6)

# flow forecast travel time plots
ttfor_summary <- summarise_draw_matrix(stm$TTForecast, flow_vec, "Flow")

(ttforecast_Sac_fig <- ggplot(ttfor_summary, aes(x = Flow, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
    geom_line(color = "black", size = 1.2) +
    geom_rug(data=d_Sac_sort,aes(x = Maxflowsac), sides = "b", alpha = 0.5) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    labs(x = "Flow", y = "Release - Sacramento travel time (days)")+
    ggtitle('Sacramento')+
    scale_x_continuous(breaks = seq(5000, 40000, by = 5000),limits = c(4000, 40000))
)

ggsave(plot=ttforecast_Sac_fig, filename = tagged_fig_file("ttforecast", "Sac"), width =10, height = 6)

ttforB_summary <- summarise_draw_matrix(stm$TTForecastT[,,1], XvecB, "Flow")

(ttforecast_But_fig <- ggplot(ttforB_summary, aes(x = Flow, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
    geom_line(color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort,aes(x = MaxflowB), sides = "b", alpha = 0.5) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    labs(x = "Flow", y = "Release - Sacramento travel time (days)")+
    ggtitle('Butte')
)

ggsave(plot=ttforecast_But_fig, filename = tagged_fig_file("ttforecast", "But"), width =10, height = 6)

ttforF_summary <- summarise_draw_matrix(stm$TTForecastT[,,2], XvecF, "Flow")

(ttforecast_Fea_fig <- ggplot(ttforF_summary, aes(x = Flow, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
    geom_line(color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort,aes(x = MaxflowF), sides = "b", alpha = 0.5) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    labs(x = "Flow", y = "Release - Sacramento travel time (days)")+
    ggtitle('Feather')
)

ggsave(plot=ttforecast_Fea_fig, filename = tagged_fig_file("ttforecast", "Fea"), width =10, height = 6)

# size forecast travel time plots
ttforSz_summary <- summarise_draw_matrix(stm$TTForecastSz[,,12], size_vec, "Size")

(ttforecastSz_Sac_fig <- ggplot(ttforSz_summary, aes(x = Size, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
    geom_line(color = "black", size = 1.2) +
    geom_rug(data=d_Sac_sort,aes(x = fish_length), sides = "b", alpha = 0.5) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    labs(x = "Size", y = "Release - Sacramento travel time (days)")+
    ggtitle('Sacramento')
)

ggsave(plot=ttforecastSz_Sac_fig, filename = tagged_fig_file("ttforecastSz", "Sac"), width =10, height = 6)

ttforSz_B_summary <- summarise_draw_matrix(stm$TTForecastTSz[,,12,1], size_vec, "Size")

(ttforecastSz_But_fig <- ggplot(ttforSz_B_summary, aes(x = Size, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
    geom_line(color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort %>% filter(!rl %in% c("1F","2F")),aes(x = fish_length), sides = "b", alpha = 0.5) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    labs(x = "Size", y = "Release - Sacramento travel time (days)")+
    ggtitle('Butte')
)

ggsave(plot=ttforecastSz_But_fig, filename = tagged_fig_file("ttforecastSz", "But"), width =10, height = 6)

ttforSz_F_summary <- summarise_draw_matrix(stm$TTForecastTSz[,,12,2], size_vec, "Size")

(ttforecastSz_Fea_fig <- ggplot(ttforSz_F_summary, aes(x = Size, y = mean)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.6) +
    geom_line(color = "black", size = 1.2) +
    geom_rug(data=d_FeaBut_sort %>% filter(rl %in% c("1F","2F")),aes(x = fish_length), sides = "b", alpha = 0.5) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    labs(x = "Size", y = "Release - Sacramento travel time (days)")+
    ggtitle('Feather')
)

ggsave(plot=ttforecastSz_Fea_fig, filename = tagged_fig_file("ttforecastSz", "Fea"), width =10, height = 6)


