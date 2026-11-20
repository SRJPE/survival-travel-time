# ============================================================
# R/survival_predict_functions.R          (srjpe_survival, 18 Nov 2026)
#
# Survival and travel time from an RST to the Sacramento station (Delta
# entry), from the posterior draws exported by
#   scripts/export/12_export_survival_pars.R  -> survival_pars.csv.gz, survival_sites.csv
#   scripts/export/12b_export_travel_time_pars.R -> travel_time_pars.csv.gz, travel_time_sites.csv
#   scripts/export/13_export_lorenzen_pars.R  -> lorenzen_pars.csv.gz
# Same formulas as the stm forecasts (SurvForecastSz_rst, TribSurvForecastSz_rst,
# TTForecastSz_rst, TTForecastTSz_rst), but at any fork length, discharge and
# release-group effect. Uses dplyr and tidyr.
#
# predict_survival(x, pars, sites, lz, small_fish, L_ref)
#   x: one row per draw x site x (week, month, ...), with columns
#      draw, site, fl_mm, cfs, e1, e2, e3
#      cfs   = monthly max discharge at the site's gauge, raw cfs
#      e1-e3 = release-group effects on the logit scale, one per reach
#              (mainstem: 3 reaches; tributary sites use e1 only). From draw_re().
#   Discharge outside the fitted range is clamped to its ends (cfs_used,
#   clamped = TRUE where that changed it).
#   small_fish: rule for fish smaller than L_ref (70 mm, the smallest tagged fish).
#     At L_ref and above all three rules are identical.
#     "extrapolate" (DEFAULT) telemetry size slope extended as a straight line
#                   on the logit scale.
#     "lorenzen"    S(L) = S(L_ref)^((L_ref / L)^c): cumulative hazard over the
#                   route scales with (L_ref / L)^c, c from the CWT gamma fit,
#                   one c per draw from lz (lorenzen_pars.csv.gz).
#     "clamp"       S(L) = S(L_ref): small fish get 70 mm survival. An upper bound.
#   TO USE LORENZEN: read lz and pass both arguments, e.g.
#     predict_survival(x, pars, sites, lz = lz, small_fish = "lorenzen")
#   To make Lorenzen the default everywhere instead, put "lorenzen" first in the
#   small_fish vector in the function definition below.
#   Not applied, by design: the CWT offset delta, the CWT group effects, and the
#   exposure term T^gamma (gamma about 0.02, so it is about 1).
#   Returns x plus cfs_used, clamped, S and size_rule ("telemetry" or the rule
#   used for that row).
#
# predict_travel_time(x, tt_pars, tt_sites, small_fish, L_ref)
#   x: draw, site, fl_mm, cfs, t1, t2, t3 (t1-t3 = release-group effects on the
#      log travel-rate scale, from draw_re(..., tt_pars = tt_pars)).
#   Travel time in days, RST to the Sacramento station: the sum over reaches of
#   exp(rate) * km / 100. Discharge is clamped as in predict_survival().
#   small_fish: rule for fish smaller than L_ref (70 mm). Identical at L_ref and above.
#     "extrapolate" (DEFAULT) telemetry size slope extended in a straight line
#                   on the log travel-rate scale.
#     "clamp"       small fish get the 70 mm travel time. A bound.
#     There is no Lorenzen-type rule for travel time.
#   Returns x plus cfs_used, clamped and
#     tt       model travel time (days): the median of individual fish
#     tt_sdlog fish-to-fish SD of log travel time (Pro_sd mainstem, Pro_sdT tributary)
#     tt_mean  mean of individual fish = tt * exp(tt_sdlog^2 / 2)
#     tt_size_rule "telemetry" or the small_fish rule used for that row
#   Individual fish travel times: tt * exp(rnorm(n(), 0, tt_sdlog)).
#
# draw_re(pars, sites, months, mode, seed, tt_pars)
#   Release-group effects for a forecast season. One effect is drawn per:
#   "none"       no effect (the stm's _nore curves)
#   "season"     draw x re_group: one effect for the whole season (what the stm does)
#   "month"      draw x month x re_group: mainstem RSTs share a month's effect;
#                Butte and Feather have their own (recommended)
#   "month_site" draw x month x site: independent for every site
#   Returns draw, month, site, e1, e2, e3 (survival). If tt_pars is given, also
#   t1, t2, t3 (travel time), drawn independently of e1-e3 (the model has no
#   correlation between them). Adding tt_pars does not change e1-e3 for a seed.
#
# predict_surv(x, pars, sites)   core survival formula, kept for scripts 12 and 16.
#   predict_survival(..., small_fish = "extrapolate") returns the same S.
# predict_surv_fl()              old name of predict_survival() (default "lorenzen"),
#   kept so earlier code still runs. New code: predict_survival().
# ============================================================

predict_survival <- function(x, pars, sites, lz = NULL,
                             small_fish = c("extrapolate", "lorenzen", "clamp"),
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
    stopifnot("lz (lorenzen_pars.csv.gz) is needed for small_fish = 'lorenzen'" = !is.null(lz),
              all(c("draw", "c") %in% names(lz)),
              all(unique(x$draw[small]) %in% lz$draw))
    out <- out |>
      left_join(select(lz, draw, .c_lz = c), by = "draw") |>
      mutate(S = if_else(fl_mm < L_ref, S^((L_ref / fl_mm)^.c_lz), S)) |>
      select(-.c_lz)
  }
  select(out, -.fl_orig)
}

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
                            small_fish = c("lorenzen", "clamp", "extrapolate"), L_ref = 70) {
  predict_survival(x, pars, sites, lz = lz, small_fish = match.arg(small_fish), L_ref = L_ref)
}

predict_travel_time <- function(x, tt_pars, tt_sites,
                                small_fish = c("extrapolate", "clamp"), L_ref = 70) {
  small_fish <- match.arg(small_fish)
  stopifnot(all(x$site %in% tt_sites$site),
            "x needs t1, t2, t3 (draw_re(..., tt_pars = tt_pars))" = all(c("t1", "t2", "t3") %in% names(x)))
  x |>
    select(-any_of(c("cfs_used", "clamped"))) |>   # recomputed below (same values if from predict_survival)
    left_join(select(tt_sites, site, .curve = curve, .km1 = km1, .km2 = km2, .km3 = km3,
                     .mult = cfs_mult, .mux = mux, .sdx = sdx, .lo = cfs_lo, .hi = cfs_hi,
                     .musz = musz, .sdsz = sdsz), by = "site") |>
    left_join(rename_with(tt_pars, \(n) paste0(".", n), -draw), by = "draw") |>
    mutate(
      .cfs     = cfs * .mult,
      cfs_used = pmin(pmax(.cfs, .lo), .hi),
      clamped  = cfs_used != .cfs,
      .z  = (cfs_used - .mux) / .sdx,
      .fl = if (small_fish == "clamp") pmax(fl_mm, L_ref) else fl_mm,
      .zs = (.fl - .musz) / .sdsz,
      .a1 = case_when(.curve == "sac"     ~ .t_sac1,
                      .curve == "butte"   ~ .t_trib_butte,
                      .curve == "feather" ~ .t_trib_feather),
      .c1 = case_when(.curve == "sac"     ~ .t_cov,
                      .curve == "butte"   ~ .t_covT_butte,
                      .curve == "feather" ~ .t_covT_feather),
      # reach 1: RST to Woodson Bridge (mainstem) or RST to Sacramento (tributary);
      # reaches 2 and 3 (Woodson-Butte, Butte-Sacramento) are mainstem only (km2 = km3 = 0 otherwise)
      tt = exp(.a1 + .c1 * .z + .t_sz * .zs + t1) * .km1 / 100 +
        if_else(.km2 > 0, exp(.t_sac2 + .t_cov * .z + .t_sz * .zs + t2) * .km2 / 100, 0) +
        if_else(.km3 > 0, exp(.t_sac3 + .t_cov * .z + .t_sz * .zs + t3) * .km3 / 100, 0),
      tt_sdlog = if_else(.curve == "sac", .pro_sd, .pro_sdT),
      tt_mean  = tt * exp(tt_sdlog^2 / 2),
      tt_size_rule = if_else(fl_mm < L_ref, small_fish, "telemetry")
    ) |>
    select(-starts_with("."))
}

draw_re <- function(pars, sites, months,
                    mode = c("none", "season", "month", "month_site"), seed = 1,
                    tt_pars = NULL) {
  mode <- match.arg(mode)
  grid <- expand_grid(draw = pars$draw, month = months, site = sites$site) |>
    left_join(select(sites, site, curve, re_group), by = "site")
  if (mode == "none") {
    out <- transmute(grid, draw, month, site, e1 = 0, e2 = 0, e3 = 0)
    if (!is.null(tt_pars)) out <- mutate(out, t1 = 0, t2 = 0, t3 = 0)
    return(out)
  }
  unit <- switch(mode,
                 season     = "re_group",
                 month      = c("month", "re_group"),
                 month_site = c("month", "site"))
  set.seed(seed)
  # u1-u3 first so e1-e3 are the same whether or not tt_pars is given
  keys <- grid |>
    distinct(across(all_of(c("draw", unit)))) |>
    mutate(u1 = rnorm(n()), u2 = rnorm(n()), u3 = rnorm(n()))
  if (!is.null(tt_pars)) {
    stopifnot("tt_pars and pars have different draws" = setequal(tt_pars$draw, pars$draw))
    keys <- mutate(keys, v1 = rnorm(n()), v2 = rnorm(n()), v3 = rnorm(n()))
  }
  out <- grid |>
    left_join(keys, by = c("draw", unit)) |>
    left_join(select(pars, draw, re_sd1, re_sd2, re_sd3, re_sdT), by = "draw")
  if (!is.null(tt_pars)) {
    out <- out |>
      left_join(select(tt_pars, draw, tt_re_sd1, tt_re_sd2, tt_re_sd3, tt_re_sdT), by = "draw") |>
      mutate(t1 = v1 * if_else(curve == "sac", tt_re_sd1, tt_re_sdT),
             t2 = if_else(curve == "sac", v2 * tt_re_sd2, 0),
             t3 = if_else(curve == "sac", v3 * tt_re_sd3, 0))
  }
  out |>
    mutate(e1 = u1 * if_else(curve == "sac", re_sd1, re_sdT),
           e2 = if_else(curve == "sac", u2 * re_sd2, 0),
           e3 = if_else(curve == "sac", u3 * re_sd3, 0)) |>
    select(draw, month, site, e1, e2, e3, any_of(c("t1", "t2", "t3")))
}

# ---- Example use in the Delta step ---------------------------------------------
# ex       <- "path/to/survival files"
# pars     <- readr::read_csv(file.path(ex, "survival_pars.csv.gz"))
# sites    <- readr::read_csv(file.path(ex, "survival_sites.csv"))
# tt_pars  <- readr::read_csv(file.path(ex, "travel_time_pars.csv.gz"))
# tt_sites <- readr::read_csv(file.path(ex, "travel_time_sites.csv"))
# lz       <- readr::read_csv(file.path(ex, "lorenzen_pars.csv.gz"))   # only for "lorenzen"
# re <- draw_re(pars, sites, months = c(11, 12, 1:5), mode = "month", tt_pars = tt_pars)
# wk: one row per site x week, with month, fl_mm (PLAD) and cfs (that month's max:
#     observed for past months, the sampled analog year for future ones; if the
#     analog year is sampled per draw, wk also has draw and you join by it)
# out <- wk |>
#   inner_join(re, by = c("site", "month"), relationship = "many-to-many") |>
#   predict_survival(pars, sites) |>                    # default: extrapolate
#   # predict_survival(pars, sites, lz = lz, small_fish = "lorenzen") |>   # Lorenzen instead
#   predict_travel_time(tt_pars, tt_sites)              # default: extrapolate
#   # predict_travel_time(tt_pars, tt_sites, small_fish = "clamp") # 70 mm bound instead
