# ============================================================
# survival_predict_functions.R
#
# Survival from RST to the Sacramento station (end of reach 3) using the
# draws exported by 05c. Same formulas as 05's SurvForecastSz_rst and
# TribSurvForecastSz_rst, but at any fork length, discharge and
# release-group effect.
#
# predict_surv(x, pars, sites)
#   x: tibble with draw, site, fl_mm, cfs, e1, e2, e3
#      cfs   = monthly max discharge at the site's gauge, raw cfs
#              (Yuba: give Feather cfs; the 0.5 factor is applied here)
#      e1-e3 = release-group effects on the logit scale, one per reach
#              (mainstem: 3 reaches; tributary fish: e1 only). From draw_re().
#   Returns x plus cfs_used (after the Yuba factor and clamping to the
#   fitted range), clamped (TRUE if clamping changed it) and S.
#   Fork length is NOT limited: below the tagged fish (about 70 mm) the
#   telemetry size slope is extended in a straight line (logit scale).
#
# predict_surv_fl(x, pars, sites, lz, small_fish, L_ref)   (added 18 Nov 2026)
#   predict_surv() with an explicit rule for fish smaller than L_ref (70 mm,
#   the smallest tagged fish). At L_ref and above it is predict_surv() exactly.
#   small_fish:
#     "lorenzen"    S(L) = S(L_ref)^((L_ref / L)^c). The cumulative hazard over
#                   the route scales with (L_ref / L)^c, c from the CWT gamma
#                   fit (G1, 07_Run_CWT_gamma.R), one c per draw from lz
#                   (07b_export_lorenzen_draws.R). Default.
#     "clamp"       S(L) = S(L_ref): small fish get 70 mm survival.
#     "extrapolate" predict_surv() as is (straight-line telemetry slope).
#   Not applied, by design: the CWT offset delta (hatchery / tag / capture
#   differences of CWT fish), the CWT group effects, and the exposure term
#   T^gamma (gamma about 0.02, so (T_L / T_ref)^gamma is about 1).
#   Returns predict_surv()'s columns plus size_rule ("telemetry" or the
#   small-fish rule used for that row).
#
# draw_re(pars, sites, months, mode, seed)
#   Release-group effects for a forecast season. One effect is drawn per:
#   "none"       no effect (05's _nore curves)
#   "season"     draw x re_group: one effect for the whole season. This is
#                what 05's stm does (one s_re per draw for every site,
#                size and discharge).
#   "month"      draw x month x re_group: mainstem RSTs share a month's
#                effect; Butte and Feather (+ Yuba) have their own
#   "month_site" draw x month x site: independent for every site
# ============================================================

predict_surv <- function(x, pars, sites) {
  stopifnot(all(x$site %in% sites$site))
  x |>
    # .curve so a 'curve' column already in x can't collide with the site table's
    left_join(select(sites, site, .curve = curve, r1, r2, r3, cfs_mult, mux, sdx,
                     cfs_lo, cfs_hi, musz, sdsz), by = "site") |>
    left_join(pars, by = "draw") |>
    mutate(
      cfs_site = cfs * cfs_mult,
      cfs_used = pmin(pmax(cfs_site, cfs_lo), cfs_hi),
      clamped  = cfs_used != cfs_site,
      z  = (cfs_used - mux) / sdx,
      zs = (fl_mm - musz) / sdsz,
      b1 = case_when(.curve == "sac"     ~ b_sac1,
                     .curve == "butte"   ~ b_trib_butte,
                     .curve == "feather" ~ b_trib_feather),
      bc = case_when(.curve == "sac"     ~ b_cov,
                     .curve == "butte"   ~ b_covT_butte,
                     .curve == "feather" ~ b_covT_feather),
      # tributary sites have r2 = r3 = 0, so those terms are 1
      S = plogis(b1     + bc * z + b_sz * zs + e1)^r1 *
          plogis(b_sac2 + bc * z + b_sz * zs + e2)^r2 *
          plogis(b_sac3 + bc * z + b_sz * zs + e3)^r3
    ) |>
    select(all_of(names(x)), cfs_used, clamped, S)
}

predict_surv_fl <- function(x, pars, sites, lz = NULL,
                            small_fish = c("lorenzen", "clamp", "extrapolate"),
                            L_ref = 70) {
  small_fish <- match.arg(small_fish)
  stopifnot(!".fl_orig" %in% names(x), !".c_lz" %in% names(x))
  small <- x$fl_mm < L_ref

  if (small_fish == "extrapolate") {
    return(predict_surv(x, pars, sites) |>
             mutate(size_rule = if_else(small, "extrapolate", "telemetry")))
  }

  # survival at max(fl, L_ref); for fish at or above L_ref this is predict_surv() exactly
  out <- x |>
    mutate(.fl_orig = fl_mm, fl_mm = pmax(fl_mm, L_ref)) |>
    predict_surv(pars, sites) |>
    mutate(fl_mm = .fl_orig, size_rule = if_else(fl_mm < L_ref, small_fish, "telemetry"))

  if (small_fish == "lorenzen") {
    stopifnot("lz (from 07b) is needed for small_fish = 'lorenzen'" = !is.null(lz),
              all(c("draw", "c") %in% names(lz)),
              all(unique(x$draw[small]) %in% lz$draw))
    out <- out |>
      left_join(select(lz, draw, .c_lz = c), by = "draw") |>
      mutate(S = if_else(fl_mm < L_ref, S^((L_ref / fl_mm)^.c_lz), S)) |>
      select(-.c_lz)
  }
  select(out, -.fl_orig)
}

draw_re <- function(pars, sites, months,
                    mode = c("none", "season", "month", "month_site"), seed = 1) {
  mode <- match.arg(mode)
  grid <- expand_grid(draw = pars$draw, month = months, site = sites$site) |>
    left_join(select(sites, site, curve, re_group), by = "site")
  if (mode == "none") {
    return(transmute(grid, draw, month, site, e1 = 0, e2 = 0, e3 = 0))
  }
  unit <- switch(mode,
                 season     = "re_group",
                 month      = c("month", "re_group"),
                 month_site = c("month", "site"))
  set.seed(seed)
  keys <- grid |>
    distinct(across(all_of(c("draw", unit)))) |>
    mutate(u1 = rnorm(n()), u2 = rnorm(n()), u3 = rnorm(n()))
  grid |>
    left_join(keys, by = c("draw", unit)) |>
    left_join(select(pars, draw, re_sd1, re_sd2, re_sd3, re_sdT), by = "draw") |>
    transmute(draw, month, site,
              e1 = u1 * if_else(curve == "sac", re_sd1, re_sdT),
              e2 = if_else(curve == "sac", u2 * re_sd2, 0),
              e3 = if_else(curve == "sac", u3 * re_sd3, 0))
}

# ---- Example use in the Delta step ---------------------------------------------
# ex    <- here("results", "survival_export")
# pars  <- read_csv(file.path(ex, "survival_pars.csv.gz"))
# sites <- read_csv(file.path(ex, "survival_sites.csv"))
# lz    <- read_csv(file.path(ex, "lorenzen_pars.csv.gz"))
# re    <- draw_re(pars, sites, months = c(11, 12, 1:5), mode = "month")
# wk: one row per site x week, with month, fl_mm (PLAD) and cfs (that month's
#     max: observed for past months, the sampled analog year for future ones;
#     if the analog year is sampled per draw, wk also has draw and you join by it)
# S  <- wk |>
#   inner_join(re, by = c("site", "month"), relationship = "many-to-many") |>
#   predict_surv_fl(pars, sites, lz, small_fish = "lorenzen")
