# ============================================================
# 03k_Run_Joint_v2_twoway.R
#
# TWO-WAY joint telemetry + CWT model with the Lorenzen link
# (scripts/Joint_v2_twoway.stan): the sensitivity check for the one-way
# model (06o). Same data and same CWT model as 06o; the only difference is
# that telemetry mortality per day M_r is computed inside the model, so the
# CWT recoveries can pull on the telemetry estimates.
#
# Run AFTER 03g (production telemetry fit) and 06o (one-way) have saved:
# starting values come from both, and the comparison at the end reads them.
# 4 chains x 1000/1000, adapt_delta 0.9, seed 1234: about the telemetry
# fit's time (4 to 6 h per chain). Saves population-level parameters,
# log_lik and the CWT terms only (small conversion).
# Output: results/fit_Joint_v2_twoway.Rdata, results/cmdstan_joint_v2_twoway/,
#         results/tables/joint_v2_*.csv
# ============================================================

rm(list = ls())

library(tidyverse)
library(here)
library(lubridate)
library(rstan)
library(cmdstanr)

source("scripts/02_GetData.R")

nchains <- 4
niter <- 2000            # 1000 warmup + 1000 sampling per chain
adapt_delta <- 0.9
seed <- 1234
L_r <- 70
c_prior <- c(mu = 1, sd = 0.25)
use_offset <- 1L
fit03g_file <- here("results", "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata")
oneway_file <- here("results", "fit_CWT_oneway_Hmed.Rdata")
stopifnot(file.exists(fit03g_file), file.exists(oneway_file))

# ---- Telemetry data: exactly as 03f / 03g -------------------------------------------
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
                 Sz=Sz,Nobs=Nobs, ObsTT=ObsTT,TTind=TTind,ReachKMT_ind=ReachKMT_ind,ReachKM_ind=ReachKM_ind)
stopifnot(all(is.finite(FL)))

# ---- CWT data: the same objects the one-way run saved -------------------------------
e <- new.env(); load(oneway_file, envir = e)
bd <- e$base_data; cwt_groups <- e$cwt_groups
ow_fit <- e$fit
rm(e)
cwt_data <- list(
  Ncwtgrp = bd$Ncwtgrp, Nobs_cwt = bd$Nobs,
  Sz_cwt = bd$Sz_cwt, CovX_cwt = bd$CovX_cwt, ReachKM_cwt = bd$ReachKM_cwt,
  ObsTT_cwt = bd$ObsTT, TTind_cwt = bd$TTind,
  Nrel = bd$Nrel, cwt_recaptures = bd$cwt_recaptures,
  Nwk = bd$Nwk, lp_mu = bd$lp_mu, lp_sd = bd$lp_sd,
  Nseg = bd$Nseg, seg_lo = bd$seg_lo, seg_hi = bd$seg_hi, seg_wk = bd$seg_wk,
  seg_start = bd$seg_start, seg_end = bd$seg_end, obs_wk = bd$obs_wk,
  log_Lratio_tel = log(L_r / FL), log_Lratio_cwt = bd$log_Lratio_cwt,
  c_mu = c_prior[["mu"]], c_sd = c_prior[["sd"]], UseOffset = use_offset)

# ---- Starting values: telemetry from 03g, CWT from the one-way fit ------------------
pmed <- function(f, p) {
  a <- rstan::extract(f, p)[[1]]
  if (length(dim(a)) <= 1) median(a) else apply(a, 2:length(dim(a)), median)
}
e <- new.env(); load(fit03g_file, envir = e); f <- e$fit
tel <- c("P_b","muPb","sdPb","S_bReach","T_bReach","S_bCov","TT_bCov","S_bSz","T_bSz",
         "S_bTrib","T_bTrib","S_bCovT","TT_bCovT","Pro_sd","Pro_sdT",
         "RE_sd","RE_sdT","TTRE_sd","TTRE_sdT","S_RE","S_REt","TT_RE","TT_RET")
tm <- set_names(map(tel, \(p) pmed(f, p)), tel)
rm(e, f); gc()
tel_init <- c(tm[c("P_b","muPb","sdPb","S_bReach","T_bReach","S_bCov","TT_bCov","S_bSz","T_bSz",
                   "S_bTrib","T_bTrib","S_bCovT","TT_bCovT","Pro_sd","Pro_sdT",
                   "RE_sd","RE_sdT","TTRE_sd","TTRE_sdT")],
              list(z_S  = sweep(tm$S_RE,   2, tm$RE_sd,    "/"),
                   z_St = sweep(tm$S_REt,  2, tm$RE_sdT,   "/"),
                   z_TT = sweep(tm$TT_RE,  2, tm$TTRE_sd,  "/"),
                   z_TTt = sweep(tm$TT_RET, 2, tm$TTRE_sdT, "/")))
cwt_init <- list(T_bReach_cwt = pmed(ow_fit, "T_bReach"), TT_bCov_cwt = pmed(ow_fit, "TT_bCov"),
                 T_bSz_cwt = pmed(ow_fit, "T_bSz"), z_TTcwt = pmed(ow_fit, "z_TT"),
                 sd_cwtTTre = pmed(ow_fit, "sd_cwtTTre"), Pro_sd_cwt = pmed(ow_fit, "Pro_sd"),
                 z_pcap = pmed(ow_fit, "z_pcap"), c = pmed(ow_fit, "c"), delta = pmed(ow_fit, "delta"),
                 z_h = pmed(ow_fit, "z_h"), sd_h = pmed(ow_fit, "sd_h"))
set.seed(seed)
inits <- map(seq_len(nchains), \(ch) {
  jit <- \(x) if (is.matrix(x)) x + matrix(rnorm(length(x), 0, 0.02), nrow(x)) else x + rnorm(length(x), 0, 0.02)
  c(map(tel_init, \(x) x), map(cwt_init, \(x) x)) |>
    imap(\(x, nm) if (nm %in% c("sdPb","RE_sd","RE_sdT","TTRE_sd","TTRE_sdT","Pro_sd","Pro_sdT",
                                "sd_cwtTTre","Pro_sd_cwt","sd_h","c")) x * exp(rnorm(length(x), 0, 0.02)) else jit(x))
})

# ---- Fit ---------------------------------------------------------------------------------
mod <- cmdstan_model("scripts/Joint_v2_twoway.stan")
model_data_cs <- c(model_data, cwt_data)[names(mod$variables()$data)] |>
  map(\(x) if (is.data.frame(x)) as.matrix(x) else x)
cs_dir <- "results/cmdstan_joint_v2_twoway"
dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
fit_cs <- mod$sample(data = model_data_cs, init = inits, chains = nchains, parallel_chains = nchains,
                     iter_warmup = niter/2, iter_sampling = niter/2, seed = seed,
                     adapt_delta = adapt_delta, refresh = 50, output_dir = cs_dir)
print(fit_cs$diagnostic_summary())

ParSaveList <- c("TT_bCov","TT_bCovT","TTRE_sd","TTRE_sdT","Pro_sd","Pro_sdT","T_bSz",
                 "muPb","sdPb","S_bReach","S_bTrib","T_bReach","T_bTrib","S_bCov","S_bCovT","S_bSz",
                 "RE_sd","RE_sdT","log_lik",
                 "c","delta","sd_h","M_r","H_tel","T_bReach_cwt","T_bSz_cwt","TT_bCov_cwt",
                 "sd_cwtTTre","Pro_sd_cwt","surv_cwt","pTT_cwt","p_rec_given_surv","pKL","log_lik_cwt")
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
save(fit, cwt_groups, file = here("results", "fit_Joint_v2_twoway.Rdata"))

# ---- Compare: telemetry vs 03g, CWT link vs one-way -------------------------------------
summ <- \(f, pars, nm) summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |>
  as.data.frame() |> rownames_to_column("parameter") |>
  transmute(model = nm, parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`, rhat = Rhat, ess = n_eff)
tel_pars <- c("S_bReach", "T_bReach", "S_bCov", "S_bSz", "T_bSz", "TT_bCov", "RE_sd", "TTRE_sd")
e <- new.env(); load(fit03g_file, envir = e)
tel_cmp <- bind_rows(summ(fit, tel_pars, "two_way"), summ(e$fit, tel_pars, "telemetry_03g"))
rm(e); gc()
link_cmp <- bind_rows(
  summ(fit, c("c", "delta", "sd_h", "M_r", "T_bReach_cwt", "T_bSz_cwt"), "two_way"),
  summ(ow_fit, c("c", "delta", "sd_h", "M_r", "T_bReach", "T_bSz"), "one_way") |>
    mutate(parameter = recode(parameter, T_bReach = "T_bReach_cwt", T_bSz = "T_bSz_cwt")))
cmp <- bind_rows(tel_cmp, link_cmp) |> mutate(across(where(is.numeric), \(x) signif(x, 3)))
write.csv(cmp, here("results", "tables", "joint_v2_twoway_vs_oneway_03g.csv"), row.names = FALSE)
cmp |> select(model, parameter, median, q2.5, q97.5, rhat) |>
  pivot_wider(names_from = model, values_from = c(median, q2.5, q97.5, rhat)) |>
  as.data.frame() |> print()

yr <- \(f, nm) tibble(model = nm, year = cwt_groups$year,
                      surv = apply(rstan::extract(f, "surv_cwt")[[1]], 2, median)) |>
  group_by(model, year) |> summarise(mean_surv = mean(surv), .groups = "drop")
by_year <- bind_rows(yr(fit, "two_way"), yr(ow_fit, "one_way")) |>
  pivot_wider(names_from = model, values_from = mean_surv)
write.csv(by_year, here("results", "tables", "joint_v2_twoway_cwt_survival_by_year.csv"), row.names = FALSE)
print(as.data.frame(by_year), digits = 3)
