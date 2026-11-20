# ============================================================
# 03c_Run_CJS_Model_centered_cmdstan.R
#
# Flora's ORIGINAL centered model (scripts/CovIndCont.stan,
# unchanged) refit to the CURRENT data, whose travel times are
# correctly aligned (A1). Her 25 Sep fit used misaligned travel
# times, so it cannot be compared directly with the non-centered
# refit. This fit separates the two changes:
#   her fit (25 Sep)  vs  this fit         = effect of aligned data
#   this fit          vs  non-centered fit = effect of non-centering
#
# Same script as 03_Run_CJS_Model_noncentered.R, with only these
# differences: ModNm = "CovIndCont" (her Stan file in scripts/),
# her inits (P_b, not z_Pb), and different output names. Still
# fitted with cmdstanr (same settings as hers: 3 chains x 500/500,
# seed 1234).
#
# Output: results/fit_CovIndCont_centered_alignedTT_MaxFlow_FL.Rdata
# Raw CSVs: results/cmdstan_03_centered/ (delete once saved).
# Nothing existing is overwritten. Expect about 5 hours.
# ============================================================

library("rstan")      # still needed: read_stan_csv() builds the stanfit that 05 expects
library("cmdstanr")   # CHANGE 2
library("purrr")

#Clear the workspace
rm(list = ls())

source("scripts/02_GetData.R")

nchains <- 3
niter <- 2000

ModNm <- "CovIndCont"   # her original (centered) Stan file in scripts/
fnext <- "_MaxFlow_FL"

# Same list as her CovIndCont branch
ParSaveList <- c("TT_bCov","TT_bCovT","TT_RE","TT_RET","TTRE_sd","TTRE_sdT","Pro_sd","Pro_sdT","T_bSz",
                 "P_b","muPb","sdPb","S_bReach","S_bTrib","T_bReach","T_bTrib","S_bCov","S_bCovT","S_bSz",
                 "RE_sd","RE_sdT","S_RE","S_REt","pTT","pTTT","pred_pcap","pred_surv","pred_survT",
                 "SurvRelSac","TT_RelSac","TT_RelSacT",'log_lik')

# her inits, unchanged
inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
               S_bReach=c(1.5,0.5,0.5,0.5) )

# ---- Her irun 15 block, unchanged ----
#Flow at RBDD for each sac fish and at release for Butte and Feather fish
UseSizeEffect <- 1
CovX <- data.frame(cbind(MaxflowSac.z,MaxflowSac.z,MaxflowSac.z,MaxflowDelta.z))
CovXT <-data.frame(cbind(MaxflowT.z,MaxflowDeltaT.z))
mux <- mean(MaxflowSac,na.rm=TRUE)
sdx <- sd(MaxflowSac,na.rm=TRUE)
muxB <- mean(MaxflowB,na.rm=TRUE)
sdxB <- sd(MaxflowB,na.rm=TRUE)
muxF <- mean(MaxflowF,na.rm=TRUE)
sdxF <- sd(MaxflowF,na.rm=TRUE)
rgwy_ind <- rep(0,Nrg)
rgwy_indT <- rep(0,NrgT)
musz <- mean(c(FL,FL_T),na.rm=TRUE)
sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
Sz <- (FL-musz)/sdsz
SzT <- (FL_T-musz)/sdsz
Xsz  <- (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
Xvec <- (seq(from =min(MaxflowSac,na.rm=TRUE),to=max(MaxflowSac,na.rm=TRUE),length.out=NsX)-mux)/sdx
XvecT <- data.frame(cbind((seq(from =min(MaxflowB),to=max(MaxflowB),length.out=NsX)-muxB)/sdxB,
                          ((seq(from =min(MaxflowF),to=max(MaxflowF),length.out=NsX)-muxF)/sdxF)))

model_data<-list(Nind=Nind,Nreaches=Nreaches,Ndetlocs=Ndetlocs,Rmult=Rmult,
                 Nyrs=Nyrs,NyrsT=NyrsT,Nrg=Nrg,Nsites=Nsites,CH=CH,
                 yrind=yrind, rgind=rgind,indsite=indsite, indsitet=indsitet,
                 rch_covind=rch_covind,
                 UseSizeEffect=UseSizeEffect,CovX=CovX,firstCap=firstCap,lastCap=lastCap,
                 Ntribs=Ntribs,NindT=NindT,NreachesT=NreachesT,NrgT=NrgT,trib_ind=trib_ind,firstCapT=firstCapT,
                 lastCapT=lastCapT,RmultT=RmultT,yrindT=yrindT,rgindT=rgindT,CHT=CHT,
                 CovXT=CovXT,SzT=SzT,NobsT=NobsT,ObsTTT=ObsTTT, TTindT=TTindT,
                 Sz=Sz,Nobs=Nobs, ObsTT=ObsTT,TTind=TTind,ReachKMT_ind=ReachKMT_ind,ReachKM_ind=ReachKM_ind)

inits <- list(inits1,inits1,inits1)


# ============================================================
# CHANGE 2: fit with cmdstanr (was rstan::stan()).
# Same settings as her call: 3 chains, iter = 1000 (rstan's default
# warmup is half, so 500 warmup + 500 sampling), seed 1234, default
# adapt_delta (0.8) and max_treedepth (10).
# ============================================================

mod <- cmdstan_model(paste0("scripts/",ModNm,".stan"))

# cmdstanr is stricter than rstan about data: pass only the variables
# the model declares (her list also carries NyrsT, Nsites, indsite,
# indsitet, which the model does not use), and data frames as matrices.
model_data_cs <- model_data[names(mod$variables()$data)] |>
  map(\(x) if (is.data.frame(x)) as.matrix(x) else x)

cs_dir <- "results/cmdstan_03_centered_long"
dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)

fit_cs <- mod$sample(data = model_data_cs, init = inits, chains = nchains,
                     parallel_chains = nchains,
                     iter_warmup = niter/2, iter_sampling = niter/2,
                     seed = 1234, refresh = 50,
                     output_dir = cs_dir)   # CSVs kept on disk, not in temp

fit_cs$diagnostic_summary()   # divergences, treedepth hits, E-BFMI per chain

# ---- Convert to an rstan stanfit (as in 06 and 09) ----
# CmdStan writes every parameter, transformed parameter and generated
# quantity, including the per-fish ones her stan() call left out
# (surv, Pcap, chi, ...). Each CSV is copied with only her ParSaveList
# columns (plus lp__ and the sampler columns), and with the header flags
# true/false rewritten as 1/0, which rstan::read_stan_csv() needs.
slim_csv <- function(f, keep) {
  x <- readLines(f)
  com <- startsWith(x, "#")
  x[com] <- gsub("= true", "= 1", gsub("= false", "= 0", x[com]))
  rows <- which(!com)                     # header row, then one row per draw
  cols <- strsplit(x[rows[1]], ",", fixed = TRUE)[[1]]
  keep_i <- which(endsWith(cols, "__") | sub("\\..*$", "", cols) %in% keep)
  x[rows] <- map_chr(x[rows], \(r) paste(strsplit(r, ",", fixed = TRUE)[[1]][keep_i], collapse = ","))
  out <- sub("\\.csv$", "_rstan.csv", f)
  writeLines(x, out)
  out
}
csvs <- map_chr(fit_cs$output_files(), slim_csv, keep = ParSaveList)
fit <- rstan::read_stan_csv(csvs)

fitnm <- "results/fit_CovIndCont_centered_alignedTT_long_MaxFlow_FL.Rdata"
save(fit,file=fitnm)
