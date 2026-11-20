# ============================================================
# 03y_flow_year_core.R
#
# Flow and year-effect comparison (team request, 8 Oct 2026).
# Fits scripts/CovIndCont_ncRE_yr.stan on the production data (02 with
# the Knights Landing Sac-Delta flow), as one of:
#   M1 = no flow, year effect      (UseFlow 0, UseYearRE 1)
#   M2 = flow + year effect        (UseFlow 1, UseYearRE 1)
# M0 = flow, no year effect = the production fit 03g (not refitted here).
#
# Run through the wrappers 03y_M1.R and 03y_M2.R (each sets model_tag),
# e.g. as two RStudio Background Jobs alongside 03g.
# Saves only population-level parameters, the year effects and log_lik,
# so the conversion is small (these fits are for comparison, not for 05).
# 4 chains x 1000/1000, seed 1234. Expect about 4 to 5 h per chain.
# Output: results/fit_flowyear_<model_tag>.Rdata (object fit),
#         results/cmdstan_flowyear_<model_tag>/ (raw CSVs).
# ============================================================

library("rstan")
library("cmdstanr")
library("purrr")

stopifnot(exists("model_tag"), model_tag %in% c("M1", "M2"))
# 02 starts with rm(list = ls()), so keep the model tag outside the workspace
Sys.setenv(FLOWYEAR_TAG = model_tag)

source("scripts/02_GetData.R")

library("rstan")
library("cmdstanr")
library("purrr")
model_tag <- Sys.getenv("FLOWYEAR_TAG")
switches <- list(M1 = c(UseFlow = 0L, UseYearRE = 1L),
                 M2 = c(UseFlow = 1L, UseYearRE = 1L))[[model_tag]]
message("Fitting ", model_tag, ": UseFlow = ", switches[["UseFlow"]], ", UseYearRE = ", switches[["UseYearRE"]])

nchains <- 4
niter <- 2000   # 1000 warmup + 1000 sampling per chain

# ---- Data exactly as 03f / 03g ----
UseSizeEffect <- 1
CovX <- data.frame(cbind(MaxflowSac.z,MaxflowSac.z,MaxflowSac.z,MaxflowDelta.z))
CovXT <-data.frame(cbind(MaxflowT.z,MaxflowDeltaT.z))
musz <- mean(c(FL,FL_T),na.rm=TRUE)
sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
Sz <- (FL-musz)/sdsz
SzT <- (FL_T-musz)/sdsz

model_data<-list(Nind=Nind,Nreaches=Nreaches,Ndetlocs=Ndetlocs,Rmult=Rmult,
                 Nyrs=Nyrs,NyrsT=NyrsT,Nrg=Nrg,Nsites=Nsites,CH=CH,
                 yrind=yrind, rgind=rgind,indsite=indsite, indsitet=indsitet,
                 rch_covind=rch_covind,
                 UseSizeEffect=UseSizeEffect,CovX=CovX,firstCap=firstCap,lastCap=lastCap,
                 Ntribs=Ntribs,NindT=NindT,NreachesT=NreachesT,NrgT=NrgT,trib_ind=trib_ind,firstCapT=firstCapT,
                 lastCapT=lastCapT,RmultT=RmultT,yrindT=yrindT,rgindT=rgindT,CHT=CHT,
                 CovXT=CovXT,SzT=SzT,NobsT=NobsT,ObsTTT=ObsTTT, TTindT=TTindT,
                 Sz=Sz,Nobs=Nobs, ObsTT=ObsTT,TTind=TTind,ReachKMT_ind=ReachKMT_ind,ReachKM_ind=ReachKM_ind,
                 UseFlow = switches[["UseFlow"]], UseYearRE = switches[["UseYearRE"]])

# her inits (as 03f), plus the year effect starting at zero
inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
               S_bReach=c(1.5,0.5,0.5,0.5), z_yr = rep(0, Nyrs), sd_yr = 0.5)
inits <- rep(list(inits1), nchains)

mod <- cmdstan_model("scripts/CovIndCont_ncRE_yr.stan")
model_data_cs <- model_data[names(mod$variables()$data)] |>
  map(\(x) if (is.data.frame(x)) as.matrix(x) else x)

cs_dir <- paste0("results/cmdstan_flowyear_", model_tag)
dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
fit_cs <- mod$sample(data = model_data_cs, init = inits, chains = nchains,
                     parallel_chains = nchains,
                     iter_warmup = niter/2, iter_sampling = niter/2,
                     seed = 1234, refresh = 50, output_dir = cs_dir)
print(fit_cs$diagnostic_summary())

# ---- Slim conversion: population-level parameters, year effects, log_lik ----
ParSaveList <- c("TT_bCov","TT_bCovT","TTRE_sd","TTRE_sdT","Pro_sd","Pro_sdT","T_bSz",
                 "muPb","sdPb","S_bReach","S_bTrib","T_bReach","T_bTrib","S_bCov","S_bCovT","S_bSz",
                 "RE_sd","RE_sdT","S_yr","sd_yr","log_lik")
slim_csv <- function(f, keep) {
  x <- readLines(f)
  com <- startsWith(x, "#")
  x[com] <- gsub("= true", "= 1", gsub("= false", "= 0", x[com]))
  rows <- which(!com)
  cols <- strsplit(x[rows[1]], ",", fixed = TRUE)[[1]]
  keep_i <- which(endsWith(cols, "__") | sub("\\..*$", "", cols) %in% keep)
  x[rows] <- map_chr(x[rows], \(r) paste(strsplit(r, ",", fixed = TRUE)[[1]][keep_i], collapse = ","))
  out <- sub("\\.csv$", "_rstan.csv", f)
  writeLines(x, out)
  out
}
fit <- rstan::read_stan_csv(map_chr(fit_cs$output_files(), slim_csv, keep = ParSaveList))
save(fit, Year_all, file = paste0("results/fit_flowyear_", model_tag, ".Rdata"))
message("Saved results/fit_flowyear_", model_tag, ".Rdata")
