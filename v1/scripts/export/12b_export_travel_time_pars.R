# ============================================================
# scripts/export/12b_export_travel_time_pars.R      (new, 18 Nov 2026)
#
# Exports the travel-time part of the telemetry model as parameter draws plus
# a site table, so runSRJPE can compute travel time with predict_travel_time()
# instead of reading the precomputed grids in the stm. Same draw index as
# survival_pars.csv.gz (both are rows of the stm, in the same order).
#
# Inputs : results/stm_MaxFlowFL_KLdelta.Rdata (10, with 11's Pro_sd, Pro_sdT)
#          results/survival_export/survival_sites.csv (12): gauges, scaling, ranges
#          scripts/refit/02_get_data.R: reach distances in km
# Writes results/survival_export/:
#   travel_time_pars.csv.gz   one row per posterior draw
#   travel_time_sites.csv     survival_sites.csv plus km1-km3 (reach lengths, km)
#
# Endpoint: RST to the Sacramento station, as in the stm's TTForecastSz_rst
# (mainstem: RST-Woodson, Woodson-Butte, Butte-Sacramento) and TTForecastTSz_rst
# (tributary: RST-Sacramento).
#
# Check: the exported pieces must reproduce the stm's TTForecastSz_rst and
# TTForecastTSz_rst (and, for the mainstem, the _kni and _tis variants that end
# reach 3 at Knights Landing and Tisdale). Those include release-group draws made
# in 10 that the stm does not store, so the check backs them out for each
# posterior draw: tributary, the one effect; mainstem, reach 1 from two RSTs at
# different distances and reach 3 from the Sacramento vs Knights Landing endpoints.
# predict_travel_time() must then match every other grid point, RST and endpoint
# to 1e-7 (relative). That tests the discharge and size terms and their scaling,
# reach 1 and reach 3 distances (Tisdale), and the sum over reaches. The reach 2
# distance always appears times its own effect, so it is checked against the
# survival distance (r2) instead.
#
# Memory: loads the stm (about 11 GB). Restart R first.
# ============================================================

rm(list = ls())
library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)

source(here("scripts", "refit", "02_get_data.R"))  # ReachKM, ReachKMrst, ReachKMTrst, ReachKM_Kni, ReachKM_Tis
source(here("R", "survival_predict_functions.R"))

stm_file <- here("results", "stm_MaxFlowFL_KLdelta.Rdata")
out_dir  <- here("results", "survival_export")
sites    <- read_csv(file.path(out_dir, "survival_sites.csv"), show_col_types = FALSE)

load(stm_file)
need <- c("T_bReach", "TT_bCov", "T_bSz", "TTRE_sd", "T_bTrib", "TT_bCovT", "TTRE_sdT",
          "Pro_sd", "Pro_sdT", "TTForecastSz_rst", "TTForecastTSz_rst")
miss <- setdiff(need, names(stm))
if (length(miss) > 0) stop("Not in the stm: ", paste(miss, collapse = ", "),
                           if (any(c("Pro_sd", "Pro_sdT") %in% miss)) " (run 11 first)")

# ---- Parameter draws ----------------------------------------------------------
tt_pars <- tibble(
  draw           = seq_len(nrow(stm$T_bReach)),
  t_sac1         = stm$T_bReach[, 1],
  t_sac2         = stm$T_bReach[, 2],
  t_sac3         = stm$T_bReach[, 3],
  t_cov          = as.vector(stm$TT_bCov),
  t_sz           = as.vector(stm$T_bSz),
  tt_re_sd1      = stm$TTRE_sd[, 1],
  tt_re_sd2      = stm$TTRE_sd[, 2],
  tt_re_sd3      = stm$TTRE_sd[, 3],
  t_trib_butte   = stm$T_bTrib[, 1],
  t_trib_feather = stm$T_bTrib[, 2],
  t_covT_butte   = stm$TT_bCovT[, 1],
  t_covT_feather = stm$TT_bCovT[, 2],
  tt_re_sdT      = stm$TTRE_sdT[, 1],
  pro_sd         = as.vector(stm$Pro_sd),
  pro_sdT        = as.vector(stm$Pro_sdT)
)
stopifnot("survival and travel-time draws differ in number" = nrow(tt_pars) == nrow(stm$S_bReach))

# ---- Site table: survival_sites plus reach lengths (km), as 10 uses them -----------
tt_sites <- sites |>
  mutate(km1 = if_else(curve == "sac", ReachKMrst[rst_index], ReachKMTrst[rst_index]),
         km2 = if_else(curve == "sac", ReachKM[2], 0),
         km3 = if_else(curve == "sac", ReachKM[3], 0),
         .after = r3)
print(select(tt_sites, site, curve, rst_index, r1, r2, r3, km1, km2, km3) |> as.data.frame())
cat("Knights Landing and Tisdale endpoints (reach 3, km):", ReachKM_Kni, ReachKM_Tis, "\n")
# survival uses distance / 100 km (r1-r3); they should agree with km / 100
km_vs_r <- tt_sites |> mutate(d1 = r1 * 100 - km1, d2 = r2 * 100 - km2, d3 = r3 * 100 - km3)
if (any(abs(c(km_vs_r$d1, km_vs_r$d2, km_vs_r$d3)) > 1))
  warning("Reach km and survival distances (r x 100) differ by more than 1 km; see table above")

# ---- Check: reproduce the stm's travel-time forecasts ------------------------------------
fl_grid <- stm$Xsz * stm$sdsz + stm$musz
cfs_grid <- list(sac     = stm$Xvec       * stm$sdx  + stm$mux,
                 butte   = stm$XvecT[[1]] * stm$sdxB + stm$muxB,
                 feather = stm$XvecT[[2]] * stm$sdxF + stm$muxF)
zgrid <- list(sac = stm$Xvec, butte = stm$XvecT[[1]], feather = stm$XvecT[[2]])
ix0 <- 1; jx0 <- 1      # grid point used to back out the release-group draws

# Tributary: one effect per draw, shared by Butte, Feather and every tributary RST (as in 10)
tr1 <- filter(tt_sites, curve != "sac") |> slice(1)
c1  <- tr1$curve; k1 <- match(c1, c("butte", "feather"))
re_T <- log(stm$TTForecastTSz_rst[, ix0, jx0, tr1$rst_index] / (tr1$km1 / 100)) -
  (stm$T_bTrib[, k1] + stm$TT_bCovT[, k1] * zgrid[[c1]][jx0] + stm$T_bSz * stm$Xsz[ix0])

# Mainstem: TT / F = g1 * km1 / 100 + g2 * km2 / 100 + g3 * km3 / 100,
# F = exp(TT_bCov * z + T_bSz * zs), gk = exp(T_bReach_k + re_k).
# Two RSTs at different km1 give g1 and the reach 2 + 3 term (beta); the Knights
# Landing endpoint (reach 3 = ReachKM_Kni instead of ReachKM[3]) separates g3.
ms <- filter(tt_sites, curve == "sac") |> arrange(km1)
ra <- ms[1, ]; rb <- ms[nrow(ms), ]
stopifnot("need two mainstem RSTs at different distances" = rb$km1 - ra$km1 > 1)
has_kt <- all(c("TTForecastSz_rst_kni", "TTForecastSz_rst_tis") %in% names(stm))
if (!has_kt) stop("stm lacks TTForecastSz_rst_kni / _tis (made by 10); rerun 10 and 11")
B <- \(arr, r) arr[, ix0, jx0, r$rst_index] /
  exp(stm$TT_bCov * stm$Xvec[jx0] + stm$T_bSz * stm$Xsz[ix0])
Ba  <- B(stm$TTForecastSz_rst, ra)
Bb  <- B(stm$TTForecastSz_rst, rb)
Bak <- B(stm$TTForecastSz_rst_kni, ra)
g1   <- (Bb - Ba) / ((rb$km1 - ra$km1) / 100)
g3   <- (Ba - Bak) / ((ReachKM[3] - ReachKM_Kni) / 100)
g2   <- (Ba - g1 * ra$km1 / 100 - g3 * ReachKM[3] / 100) / (ReachKM[2] / 100)
stopifnot("backed-out mainstem terms not positive" = all(g1 > 0), all(g2 > 0), all(g3 > 0))
re_M <- tibble(draw = tt_pars$draw,
               t1 = log(g1) - stm$T_bReach[, 1],
               t2 = log(g2) - stm$T_bReach[, 2],
               t3 = log(g3) - stm$T_bReach[, 3])
re_all <- bind_rows(mutate(re_M, sac = TRUE),
                    tibble(draw = tt_pars$draw, t1 = re_T, t2 = 0, t3 = 0, sac = FALSE))

# Sites to check: every RST to Sacramento, plus mainstem RSTs to Knights Landing and Tisdale
chk_sites <- bind_rows(
  mutate(tt_sites, endpoint_chk = "sac"),
  filter(tt_sites, curve == "sac") |> mutate(km3 = ReachKM_Kni, endpoint_chk = "kni"),
  filter(tt_sites, curve == "sac") |> mutate(km3 = ReachKM_Tis, endpoint_chk = "tis")) |>
  mutate(rst = site, site = paste(site, endpoint_chk, sep = "_to_"))
ref_arr <- list(sac = stm$TTForecastSz_rst, kni = stm$TTForecastSz_rst_kni,
                tis = stm$TTForecastSz_rst_tis, trib = stm$TTForecastTSz_rst)
nsz <- length(fl_grid); nsx <- length(stm$Xvec)

chk <- expand_grid(site = chk_sites$site, ix = c(4, 7, 12, 18, 23), jx = c(2, 4, 13, 22, 24),
                   draw = tt_pars$draw) |>
  filter(ix <= nsz, jx <= nsx) |>
  left_join(select(chk_sites, site, curve, rst_index, endpoint_chk), by = "site") |>
  mutate(sac   = curve == "sac",
         arr   = if_else(sac, endpoint_chk, "trib"),
         fl_mm = fl_grid[ix],
         cfs   = map2_dbl(curve, jx, \(cv, j) cfs_grid[[cv]][j])) |>
  mutate(.by = arr, ref = ref_arr[[arr[1]]][cbind(draw, ix, jx, rst_index)]) |>
  left_join(re_all, by = c("draw", "sac"))

chk_out <- predict_travel_time(select(chk, draw, site, fl_mm, cfs, t1, t2, t3), tt_pars, chk_sites) |>
  mutate(ref = chk$ref)
stopifnot("discharge grid points should not be clamped" = !any(chk_out$clamped))
max_diff <- chk_out |>
  summarise(max_rel_diff = max(abs(tt / ref - 1)), .by = site)
print(as.data.frame(max_diff), row.names = FALSE)
stopifnot(all(max_diff$max_rel_diff < 1e-7))
cat("Check passed: exported draws reproduce the stm's travel-time forecasts at every RST.\n")

# ---- Summary at the median discharge and a 70 mm fish, no release-group effect ----------
med <- tt_sites |>
  transmute(site, fl_mm = 70, cfs = (cfs_lo + cfs_hi) / 2) |>
  cross_join(select(tt_pars, draw)) |>
  mutate(t1 = 0, t2 = 0, t3 = 0) |>
  predict_travel_time(tt_pars, tt_sites) |>
  summarise(tt_median_days = median(tt), tt_lo95 = quantile(tt, 0.025), tt_hi95 = quantile(tt, 0.975),
            sdlog = median(tt_sdlog), .by = site)
cat("\nTravel time to Sacramento, 70 mm fish, mid-range discharge, no release-group effect:\n")
print(mutate(med, across(where(is.numeric), \(v) round(v, 2))) |> as.data.frame())

# ---- Write ----------------------------------------------------------------------
write_csv(tt_pars,  file.path(out_dir, "travel_time_pars.csv.gz"))
write_csv(tt_sites, file.path(out_dir, "travel_time_sites.csv"))
cat("Wrote", nrow(tt_pars), "draws and", nrow(tt_sites), "sites to", out_dir, "\n")
