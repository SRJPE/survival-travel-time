# ============================================================
# 03j_Run_Joint_v1.R
#
# Joint telemetry + CWT model, version 1 (for the 15 Nov handoff).
# Stan model: scripts/Joint_v1.stan =
#   telemetry: CovIndCont_ncRE.stan (03f), unchanged
#   + CWT: CWT_weekly.stan (06w), sharing the telemetry reach survival and
#     travel-time parameters (Perry-style exchangeability; Flora's
#     logit-linear size effect). KL lies about 119 km into reach 3
#     (Butte-Sac), so CWT fish use reaches 1, 2 and part of 3.
#
# Data:
#   telemetry exactly as 03f (02 with study x month release groups);
#   CWT exactly as 06w (same groups, recoveries, weeks, BT-SPAS-X
#   priors), except:
#     - flow on the telemetry reference scale (flow_mu_sac, flow_sd_sac
#       from 02, fix A2), not the CWT group mean/sd;
#     - size on the telemetry scale (musz, sdsz), so S_bSz and T_bSz
#       mean the same thing for both tag types;
#     - distance split by telemetry reach: reach 1 = RST site to Woodson
#       (02's dist_rstwoodson), reach 2 = Woodson-Butte (telemetry
#       median, 88 km), reach 3 = the rest to KL (06's distances).
#
# Settings: 4 chains x 1000/1000 (as 03f), adapt_delta 0.9, seed 1234.
# Expect about the 03f time (3.5 h per chain) or a bit more.
# Output: results/fit_Joint_v1_MaxFlow_FL.Rdata (object `fit`, readable
# by 05 like the 03f fit), results/cmdstan_joint_v1/ (raw CSVs),
# results/tables/joint_v1_*.csv.
# Needs: data/kl_weekly_pcap.csv; for the comparison at the end,
# results/fit_CovIndCont_studymonth_ncRE_MaxFlow_FL.Rdata (03f) and
# results/fit_CWT_weekly.Rdata (06w).
# ============================================================

library(rstan)
library(cmdstanr)
library(tidyverse)
library(here)
library(lubridate)

rm(list = ls())

source("scripts/02_GetData.R")

nchains <- 4
niter <- 2000            # 1000 warmup + 1000 sampling per chain
adapt_delta <- 0.9
seed <- 1234

# ---- Telemetry: exactly as 03f ---------------------------------------------------------
UseSizeEffect <- 1
CovX <- data.frame(cbind(MaxflowSac.z,MaxflowSac.z,MaxflowSac.z,MaxflowDelta.z))
CovXT <-data.frame(cbind(MaxflowT.z,MaxflowDeltaT.z))
mux <- mean(MaxflowSac,na.rm=TRUE)
sdx <- sd(MaxflowSac,na.rm=TRUE)
muxB <- mean(MaxflowB,na.rm=TRUE)
sdxB <- sd(MaxflowB,na.rm=TRUE)
muxF <- mean(MaxflowF,na.rm=TRUE)
sdxF <- sd(MaxflowF,na.rm=TRUE)
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

# ---- CWT: as 06w -------------------------------------------------------------------------
release_distance_km <- c(Battle = 300.443064082, RBDD = 240.3631)   # to KL (06; = 02's dist_rstknights)
kl_file <- here("data", "kl_weekly_pcap.csv")
if (!file.exists(kl_file)) stop("Copy kl_weekly_pcap.csv into data/")

first_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

cwt_groups_06 <- drerelrec %>%
  filter(!is.na(relloc_area), relloc_area %in% names(release_distance_km),
         release_group_id != 873) %>%
  group_by(release_group_id) %>%
  summarise(
    release_location_name = first(release_location_name),
    relloc_area = first(relloc_area),
    year = first(year),
    month = first(month),
    mid_release_date = first(mid_release_date),
    Nrel = first(group_total_marked_N),
    avg_length = first(avg_length),
    Maxflowsac = first_finite(monthly_max_flow),
    cwt_recaptures = n_distinct(
      if_else(is.na(tag_code), NA_character_, paste(tag_code, date, forklength, sep = "_")),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ReachKM_cwt = unname(release_distance_km[relloc_area])) %>%
  filter(is.finite(Nrel), Nrel > 0, is.finite(avg_length),
         is.finite(Maxflowsac), is.finite(ReachKM_cwt)) %>%
  arrange(year, month, release_group_id)

kl_weekly <- read.csv(kl_file) |>
  transmute(kl_year = as.integer(str_sub(fit, -4)), run_week = as.integer(week), lp_mean, lp_sd)

kl_week_of <- function(d) {
  jw <- as.integer(lubridate::week(d))
  yr <- as.integer(lubridate::year(d))
  tibble(kl_year = if_else(jw >= 45L, yr + 1L, yr),
         run_week = case_when(jw >= 45L ~ jw - 44L, jw <= 22L ~ jw + 9L, .default = NA_integer_))
}
season_end <- function(kl_year) make_date(kl_year, 1, 1) + days(7 * 22 - 1)

grp_season <- cwt_groups_06 |>
  mutate(rel_date = as_date(mid_release_date), kl_week_of(rel_date),
         keep = !is.na(run_week) & kl_year %in% kl_weekly$kl_year)

recs <- drerelrec |>
  filter(release_group_id %in% grp_season$release_group_id, !is.na(tag_code)) |>
  distinct(release_group_id, tag_code, date, forklength, .keep_all = TRUE) |>
  select(release_group_id, tag_code, date, forklength, cwt_tt) |>
  left_join(grp_season |> select(release_group_id, rel_kl_year = kl_year), by = "release_group_id") |>
  mutate(date = as_date(date), in_season = is.na(date) | date <= season_end(rel_kl_year))
stopifnot(nrow(recs) == sum(cwt_groups_06$cwt_recaptures))

cwt_groups <- grp_season |>
  filter(keep) |>
  select(-cwt_recaptures) |>
  left_join(recs |> filter(in_season) |> count(release_group_id, name = "cwt_recaptures"),
            by = "release_group_id") |>
  mutate(cwt_recaptures = replace_na(cwt_recaptures, 0L)) |>
  arrange(year, month, release_group_id) |>
  mutate(grp = row_number())

# JOINT CHANGE: telemetry scales for flow and size; distance split by telemetry reach
km_reach1 <- c(Battle = dist_rstwoodson[2], RBDD = dist_rstwoodson[6])   # RST site to Woodson (02)
km_reach2 <- median(ReachKM_ind[, 2])                                       # Woodson-Butte, telemetry
cwt_groups <- cwt_groups |>
  mutate(Sz_cwt = (avg_length - musz) / sdsz,
         CovX_cwt = (Maxflowsac - flow_mu_sac) / flow_sd_sac,
         km1 = unname(km_reach1[relloc_area]),
         km2 = km_reach2,
         km3 = ReachKM_cwt - km1 - km2)
print(distinct(cwt_groups, relloc_area, km1, km2, km3, ReachKM_cwt))
stopifnot(all(cwt_groups$km3 > 100 & cwt_groups$km3 < 170))   # KL inside reach 3 (170 km)
message("CWT size range (mm): ", paste(range(cwt_groups$avg_length), collapse = " to "),
        "; telemetry: ", paste(round(range(c(FL, FL_T), na.rm = TRUE)), collapse = " to "))

make_segments <- function(grp, rel_date, kl_year) {
  tibble(day = seq(rel_date, season_end(kl_year), by = "day")) |>
    mutate(off = as.numeric(day - rel_date), kl_week_of(day)) |>
    group_by(kl_year, run_week) |>
    summarise(lo = pmax(min(off) - 0.5, 0), hi = max(off) + 0.5, .groups = "drop") |>
    arrange(lo) |>
    mutate(grp = grp)
}
segs <- cwt_groups |> select(grp, rel_date, kl_year) |> pmap(make_segments) |> list_rbind()
wk_tbl <- segs |>
  distinct(kl_year, run_week) |>
  arrange(kl_year, run_week) |>
  left_join(kl_weekly, by = c("kl_year", "run_week")) |>
  mutate(wk = row_number())
stopifnot(!anyNA(wk_tbl$lp_mean))
segs <- segs |>
  left_join(wk_tbl |> select(kl_year, run_week, wk), by = c("kl_year", "run_week")) |>
  arrange(grp, lo) |>
  mutate(seg = row_number())
seg_idx <- segs |> group_by(grp) |> summarise(seg_start = min(seg), seg_end = max(seg))
stopifnot(identical(seg_idx$grp, cwt_groups$grp))

cwt_tt_obs <- recs |>
  filter(in_season, is.finite(cwt_tt), cwt_tt > 0) |>
  distinct(release_group_id, tag_code, date, forklength, cwt_tt) |>
  inner_join(cwt_groups |> select(release_group_id, grp), by = "release_group_id") |>
  inner_join(segs |> select(grp, lo, hi, wk), by = join_by(grp, cwt_tt >= lo, cwt_tt < hi)) |>
  arrange(grp, date, tag_code)
message("CWT groups: ", nrow(cwt_groups), "; recoveries: ", sum(cwt_groups$cwt_recaptures),
        "; travel times: ", nrow(cwt_tt_obs))   # 06w: 357, 3013, 3013

cwt_data <- list(
  Ncwtgrp = nrow(cwt_groups), Nobs_cwt = nrow(cwt_tt_obs),
  Rmult_cwt = as.matrix(cwt_groups[, c("km1", "km2", "km3")]) / 100,
  Sz_cwt = cwt_groups$Sz_cwt, CovX_cwt = cwt_groups$CovX_cwt,
  ObsTT_cwt = cwt_tt_obs$cwt_tt, TTind_cwt = cwt_tt_obs$grp,
  Nrel = as.integer(round(cwt_groups$Nrel)), cwt_recaptures = as.integer(cwt_groups$cwt_recaptures),
  Nwk = nrow(wk_tbl), lp_mu = wk_tbl$lp_mean, lp_sd = wk_tbl$lp_sd,
  Nseg = nrow(segs), seg_lo = segs$lo, seg_hi = segs$hi, seg_wk = segs$wk,
  seg_start = seg_idx$seg_start, seg_end = seg_idx$seg_end, obs_wk = cwt_tt_obs$wk)

save(cwt_data, cwt_groups, cwt_tt_obs, segs, wk_tbl, musz, sdsz,
     file = here("results", "joint_v1_cwt_data.Rdata"))

# ---- Fit ----------------------------------------------------------------------------------
mod <- cmdstan_model("scripts/Joint_v1.stan")
model_data_cs <- c(model_data, cwt_data)[names(mod$variables()$data)] |>
  map(\(x) if (is.data.frame(x)) as.matrix(x) else x)

# her telemetry inits, plus sensible starts for the travel-time and CWT terms
# (random starts can put CWT travel times far outside the KL season weeks);
# small jitter per chain so chains still start apart
set.seed(seed)
inits <- map(seq_len(nchains), \(ch) list(
  P_b = matrix(2.2, nrow = Nyrs, ncol = Nreaches), muPb = rep(0, Nreaches), sdPb = rep(1, Nreaches),
  S_bReach = c(1.5, 0.5, 0.5, 0.5) + rnorm(4, 0, 0.1),
  T_bReach = rep(1, Nreaches) + rnorm(Nreaches, 0, 0.1),
  TT_bCov = 0, T_bSz = 0, S_bCov = 0, S_bSz = 0,
  z_Scwt = rep(0, cwt_data$Ncwtgrp), sd_cwtSre = 0.5,
  z_TTcwt = rep(0, cwt_data$Ncwtgrp), sd_cwtTTre = 0.2,
  Pro_sd_cwt = 0.45, z_pcap = rep(0, cwt_data$Nwk)))

cs_dir <- "results/cmdstan_joint_v1"
dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
fit_cs <- mod$sample(data = model_data_cs, init = inits, chains = nchains, parallel_chains = nchains,
                     iter_warmup = niter/2, iter_sampling = niter/2, seed = seed,
                     adapt_delta = adapt_delta, refresh = 50, output_dir = cs_dir)
fit_cs$diagnostic_summary()

# ---- Convert to an rstan stanfit (as 03f), keeping her list plus the CWT terms ----------
ParSaveList <- c("TT_bCov","TT_bCovT","TT_RE","TT_RET","TTRE_sd","TTRE_sdT","Pro_sd","Pro_sdT","T_bSz",
                 "P_b","muPb","sdPb","S_bReach","S_bTrib","T_bReach","T_bTrib","S_bCov","S_bCovT","S_bSz",
                 "RE_sd","RE_sdT","S_RE","S_REt","pTT","pTTT","pred_pcap","pred_surv","pred_survT",
                 "SurvRelSac","TT_RelSac","TT_RelSacT",'log_lik',
                 # CWT (joint v1)
                 "sd_cwtSre","sd_cwtTTre","Pro_sd_cwt","S_cwtre","TT_cwtre","logit_pcap",
                 "surv_cwt","pTT_cwt","p_rec_given_surv","pKL","log_lik_cwt")
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
save(fit, file = "results/fit_Joint_v1_MaxFlow_FL.Rdata")

# ---- Compare with telemetry only (03f) and CWT only (06w) ------------------------------
shared <- c("S_bReach", "T_bReach", "S_bCov", "S_bSz", "T_bSz", "TT_bCov",
            "RE_sd", "TTRE_sd", "Pro_sd", "sd_cwtSre", "sd_cwtTTre", "Pro_sd_cwt")
summ <- \(f, pars, nm) summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |>
  as.data.frame() |> rownames_to_column("parameter") |>
  transmute(model = nm, parameter, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`, rhat = Rhat, ess = n_eff)
e <- new.env(); load(here("results", "fit_CovIndCont_studymonth_ncRE_MaxFlow_FL.Rdata"), envir = e)
cmp <- bind_rows(summ(fit, shared, "joint_v1"),
                 summ(e$fit, intersect(shared, names(e$fit@par_dims)), "telemetry_03f")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
rm(e); gc()
write.csv(cmp, here("results", "tables", "joint_v1_vs_03f_parameters.csv"), row.names = FALSE)
cmp |> select(model, parameter, median, q2.5, q97.5, rhat) |>
  pivot_wider(names_from = model, values_from = c(median, q2.5, q97.5, rhat)) |>
  as.data.frame() |> print()

# CWT survival to KL by year: joint vs 06w weekly (CWT only)
e <- new.env(); load(here("results", "fit_CWT_weekly.Rdata"), envir = e)
yr_surv <- \(f, nm) tibble(model = nm, year = cwt_groups$year,
                           surv = apply(rstan::extract(f, "surv_cwt")[[1]], 2, median)) |>
  group_by(model, year) |> summarise(mean_surv = mean(surv), .groups = "drop")
by_year <- bind_rows(yr_surv(fit, "joint_v1"), yr_surv(e$fit_cwt_weekly, "cwt_only_06w")) |>
  pivot_wider(names_from = model, values_from = mean_surv) |>
  left_join(cwt_groups |> group_by(year) |>
              summarise(groups = n(), recoveries = sum(cwt_recaptures)), by = "year")
rm(e); gc()
write.csv(by_year, here("results", "tables", "joint_v1_cwt_survival_by_year.csv"), row.names = FALSE)
print(as.data.frame(by_year), digits = 3)
