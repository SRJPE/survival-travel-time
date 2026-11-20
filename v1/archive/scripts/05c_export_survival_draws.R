# ============================================================
# 05c_export_survival_draws.R
#
# Exports the survival model as parameter draws plus a site table, so the
# JPE Delta step can compute survival at any fork length, discharge and
# release-group effect (by week, month and site). No refit: parameters come
# from the stm saved by 05, distances from 02.
#
# Writes results/survival_export/:
#   survival_pars.csv.gz   one row per posterior draw (same draw index for all sites)
#   survival_sites.csv     one row per RST: curve, distances, gauge, scaling, fitted range
# Functions to use them: scripts/survival_predict_functions.R
#
# Endpoint: RST to the Sacramento station (end of reach 3), as in 05's
# SurvForecastSz_rst / TribSurvForecastSz_rst. Sac-Delta reach not included.
#
# Check at the end: the exported pieces must reproduce 05's _nore forecasts.
# ============================================================

rm(list = ls())
library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)

source("scripts/02_GetData.R")  # Rmultrst, Rmult, RmultSac, RmultTrst, RmultTrib, TribForRST, Nrst, NrstT, FL, FL_T
source(here("scripts", "survival_predict_functions.R"))

stm_file <- here("results", "stm_MaxFlowFL_KLdelta.Rdata")
out_dir  <- here("results", "survival_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# RST index -> site, same codes and order as Liz's extract_stm_inputs.R and 11_export.
# 02 has Nrst = 6 mainstem RSTs (5 Redding, 6 Red Bluff, not used) and NrstT = 3
# tributary RSTs with TribForRST = 1, 2, 2 (Hallwood and LFC both on the Feather curve).
rst_map <- tribble(
  ~site,        ~curve,    ~rst_index,
  "ucc",        "sac",     1,
  "ubc",        "sac",     2,
  "mill creek", "sac",     3,
  "deer creek", "sac",     4,
  "okie dam",   "butte",   1,
  "hallwood",   "feather", 2,
  "lfc",        "feather", 3
)
trib_curves <- c("butte", "feather")   # S_bTrib / S_bCovT column order (trib_ind 1 = Butte, 2 = Feather)

load(stm_file)
cat("Nrst =", Nrst, "; NrstT =", NrstT, "; TribForRST =", TribForRST, "\n")
trib_rows <- filter(rst_map, curve != "sac")
stopifnot(Nrst >= 4, NrstT == 3,
          identical(trib_curves[TribForRST[trib_rows$rst_index]], trib_rows$curve))

# ---- Parameter draws ----------------------------------------------------------
pars <- tibble(
  draw           = seq_len(nrow(stm$S_bReach)),
  b_sac1         = stm$S_bReach[, 1],
  b_sac2         = stm$S_bReach[, 2],
  b_sac3         = stm$S_bReach[, 3],
  b_cov          = as.vector(stm$S_bCov),
  b_sz           = as.vector(stm$S_bSz),
  re_sd1         = stm$RE_sd[, 1],
  re_sd2         = stm$RE_sd[, 2],
  re_sd3         = stm$RE_sd[, 3],
  b_trib_butte   = stm$S_bTrib[, 1],
  b_trib_feather = stm$S_bTrib[, 2],
  b_covT_butte   = stm$S_bCovT[, 1],
  b_covT_feather = stm$S_bCovT[, 2],
  re_sdT         = stm$RE_sdT[, 1]
)

# ---- Site table -----------------------------------------------------------------
gauge_tab <- tribble(
  ~curve,    ~gauge,                                                     ~mux,      ~sdx,
  "sac",     "USGS 11390500 Wilkins Slough (SRJPEdata 'knights landing')", stm$mux,  stm$sdx,
  "butte",   "USGS 11390000 Butte Creek near Chico",                      stm$muxB, stm$sdxB,
  "feather", "SRJPEdata Feather composite (max of reporting gauges)",     stm$muxF, stm$sdxF
) |>
  # fitted range = ends of each forecast grid, back on the cfs scale
  mutate(xz = list(stm$Xvec, stm$XvecT[[1]], stm$XvecT[[2]]),
         cfs_lo = map_dbl(xz, min) * sdx + mux,
         cfs_hi = map_dbl(xz, max) * sdx + mux) |>
  select(-xz)
print(select(gauge_tab, curve, cfs_lo, cfs_hi))

sites <- rst_map |>
  mutate(
    re_group = curve,   # mainstem RSTs share; Hallwood and LFC share the Feather group
    r1 = if_else(curve == "sac", Rmultrst[rst_index], RmultTrst[rst_index]),
    r2 = if_else(curve == "sac", Rmult[1, 2], 0),
    r3 = if_else(curve == "sac", Rmult[1, 3], 0),
    # Hallwood (Yuba): Josh's For_Input used 0.5 x Feather discharge; Liz's site table
    # uses Feather discharge as is. 1 here; set to 0.5 if the team keeps Josh's rule.
    cfs_mult = 1
  ) |>
  left_join(gauge_tab, by = "curve") |>
  mutate(musz = stm$musz, sdsz = stm$sdsz, endpoint = "Sacramento station")

print(as.data.frame(sites))
cat("Fork-length range of tagged fish (mm):", range(c(FL, FL_T), na.rm = TRUE), "\n")
cat("Release-group SDs (median): Sac reaches",
    round(map_dbl(select(pars, re_sd1:re_sd3), median), 2),
    "; tributary", round(median(pars$re_sdT), 2), "\n")

# ---- Check: reproduce 05's _nore forecasts from the exported pieces ---------------
check_sites <- tibble(
  site = c("chk_sac", "chk_butte", "chk_feather"),
  curve = c("sac", "butte", "feather"), re_group = curve,
  r1 = c(RmultSac, RmultTrib), r2 = c(Rmult[1, 2], 0, 0), r3 = c(Rmult[1, 3], 0, 0),
  cfs_mult = 1
) |>
  left_join(gauge_tab, by = "curve") |>
  mutate(musz = stm$musz, sdsz = stm$sdsz)

grid_cfs <- list(sac     = stm$Xvec       * stm$sdx  + stm$mux,
                 butte   = stm$XvecT[[1]] * stm$sdxB + stm$muxB,
                 feather = stm$XvecT[[2]] * stm$sdxF + stm$muxF)
ref <- list(sac     = stm$SurvForecastSz_nore,
            butte   = stm$TribSurvForecastSz_nore[, , , 1],
            feather = stm$TribSurvForecastSz_nore[, , , 2])
fl_grid <- stm$Xsz * stm$sdsz + stm$musz

chk <- expand_grid(curve = names(grid_cfs), ix = c(7, 12, 18), jx = c(4, 13, 22),
                   draw = pars$draw) |>
  mutate(.by = curve,
         site  = paste0("chk_", curve),
         cfs   = grid_cfs[[curve[1]]][jx],
         fl_mm = fl_grid[ix],
         ref   = ref[[curve[1]]][cbind(draw, ix, jx)],
         e1 = 0, e2 = 0, e3 = 0)

chk_out <- predict_surv(select(chk, -ref), pars, check_sites) |>
  mutate(ref = chk$ref)
max_diff <- chk_out |> summarise(max_abs_diff = max(abs(S - ref)), .by = curve)
print(max_diff)
stopifnot(all(max_diff$max_abs_diff < 1e-8))
cat("Check passed: exported draws reproduce 05's _nore forecasts.\n")

# ---- Write ----------------------------------------------------------------------
write_csv(pars,  file.path(out_dir, "survival_pars.csv.gz"))
write_csv(sites, file.path(out_dir, "survival_sites.csv"))
cat("Wrote", nrow(pars), "draws and", nrow(sites), "sites to", out_dir, "\n")
