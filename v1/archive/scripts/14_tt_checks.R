# ============================================================
# 14_tt_checks.R
# Travel-time checks before the 15 Nov handoff. No model fitting.
#
#  A. Data only: is the discharge effect on travel time (TT_bCov = +0.52 in
#     03g, higher discharge = slower) present WITHIN release sites, or only
#     BETWEEN them? Release-group medians of observed travel time, release to
#     Sacramento station, regressed on the same z-scored discharge the model
#     uses (MaxflowSac.z = CovX), so slopes are on the TT_bCov scale.
#  B. Fit check: model-predicted vs observed median travel time by release
#     group (new stm, TT_reach), and residuals vs discharge and fork length.
#  C. Forecast curves: RST to Delta entry travel time vs discharge (at mean
#     size) and vs size (at mean discharge), new stm vs the team's old stm.
#  D. Fish-level spread (Pro_sd), only once Pro_sd has been added to the stm.
#
# Inputs : scripts/02_GetData.R (sourced; clears the workspace)
#          results/handoff/stm.Rdata  (new stm, thinned to 1,500 draws)
#          data/stm.Rdata             (team's previous stm)
# Outputs: results/figures/14_tt_checks/*.png
#          results/tables/14_tt_checks/*.csv
# Restart R first. The two stm objects are loaded one at a time.
# ============================================================

source(here::here("scripts", "02_GetData.R"))   # d_Sac_sort, TTfR, ReachKM_ind, MaxflowSac(.z), FL, rgind

library(tidyverse)
library(here)

fig_dir <- here("results", "figures", "14_tt_checks")
tab_dir <- here("results", "tables", "14_tt_checks")
walk(c(fig_dir, tab_dir), \(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))

stm_new_file <- here("results", "handoff", "stm.Rdata")
stm_old_file <- here("data", "stm.Rdata")

probs   <- c(0.1, 0.5, 0.9)                     # 80% intervals throughout
stm_col <- c(new = "#0072B2", old = "#E69F00")  # Okabe-Ito
stm_lab <- c(new = "New (03g, Wilkins Slough)", old = "Old (team stm)")
theme_set(theme_bw(base_size = 12))

qt <- function(x) {
  v <- quantile(x, probs, names = FALSE, na.rm = TRUE)
  tibble(lo = v[1], med = v[2], hi = v[3])
}

# ============================================================
# A. Data only: discharge vs observed travel time, within and between sites
# ============================================================
col_or_na <- \(nm) if (nm %in% names(d_Sac_sort)) d_Sac_sort[[nm]] else NA
site_col  <- if ("rel_site" %in% names(d_Sac_sort)) "rel_site" else "release_location"

fish <- tibble(
  row    = seq_len(Nind),
  study  = d_Sac_sort$StudyID,
  site   = d_Sac_sort[[site_col]],
  year   = col_or_na("year"),
  rg     = rgind,
  flow   = MaxflowSac,
  flow_z = as.numeric(MaxflowSac.z),
  fl     = FL,
  km     = rowSums(ReachKM_ind[, 1:3]),          # release to Sacramento station
  tt     = TTfR$TTfR3                             # observed days, release to Sacramento station
) |>
  filter(!is.na(tt), tt > 0, km > 0) |>
  mutate(log_rate = log(tt / km * 100))           # log days per 100 km, the TT model's scale

# One row per release group (study x month): discharge is constant within a group
grp <- fish |>
  summarise(n = n(), study = first(study), site = first(site), year = first(year),
            flow = first(flow), flow_z = first(flow_z), fl = mean(fl, na.rm = TRUE),
            obs_med = median(tt), log_rate = median(log_rate), .by = rg) |>
  mutate(fl_z = as.numeric(scale(fl)))

flow_slope <- function(fit, label) {
  ci <- confint(fit)["flow_z", ]
  tibble(model = label, slope = coef(fit)[["flow_z"]], lo95 = ci[[1]], hi95 = ci[[2]],
         n_groups = nobs(fit))
}

# Unweighted: group-to-group variation dominates, so each group counts once
forms <- list(
  "pooled (between + within sites)" = log_rate ~ flow_z,
  "within site"                     = log_rate ~ flow_z + site,
  "within site, adj. fork length"   = log_rate ~ flow_z + site + fl_z
)
if (n_distinct(grp$year, na.rm = TRUE) > 1) {
  forms[["within site and year, adj. fork length"]] <- log_rate ~ flow_z + site + factor(year) + fl_z
}
slopes <- imap(forms, \(f, nm) flow_slope(lm(f, data = grp), nm)) |> list_rbind()

# Separate slope for each site with at least 3 distinct discharges
site_slopes <- grp |>
  filter(n_distinct(flow_z) >= 3, .by = site) |>
  nest(.by = site) |>
  mutate(res = map2(data, site, \(d, s) flow_slope(lm(log_rate ~ flow_z, data = d), paste("site:", s)))) |>
  select(res) |>
  unnest(res)

# Figure A1: observed travel time vs discharge, one panel per release site
pA1 <- ggplot(grp, aes(flow / 1000, exp(log_rate))) +
  geom_point(data = fish, colour = "grey75", alpha = 0.25, size = 0.6) +
  geom_point(aes(size = n), shape = 21, fill = "#0072B2", colour = "black") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D55E00", linewidth = 1) +
  scale_y_log10() +
  facet_wrap(~ site) +
  labs(x = "Monthly max discharge in release month, Wilkins Slough (kcfs)",
       y = "Days per 100 km, release to Sacramento station (log scale)",
       size = "Fish in group",
       caption = "Grey: individual fish. Blue: release-group medians. Orange: within-site fit to group medians.")
ggsave(file.path(fig_dir, "A1_tt_vs_discharge_by_site.png"), pA1, width = 10, height = 7, dpi = 300)

# Figure A2: what discharge is confounded with (site and fish size)
pA2 <- ggplot(grp, aes(fl, flow / 1000)) +
  geom_point(aes(size = n), shape = 21, fill = "#0072B2", colour = "black") +
  facet_wrap(~ site) +
  labs(x = "Mean fork length of release group (mm)",
       y = "Monthly max discharge, Wilkins Slough (kcfs)", size = "Fish in group")
ggsave(file.path(fig_dir, "A2_discharge_vs_size_by_site.png"), pA2, width = 10, height = 7, dpi = 300)

# ============================================================
# Load the NEW stm once: parameters (A), fit check (B), curves (C), Pro_sd (D)
# ============================================================
e <- new.env(); load(stm_new_file, envir = e); s_new <- e$stm; rm(e); invisible(gc())
need <- c("TT_bCov", "T_bSz", "TT_bCovT", "TT_reach", "TTForecastSz_rst",
          "Xvec", "Xsz", "mux", "sdx", "musz", "sdsz")
stopifnot("new stm is missing objects" = all(need %in% names(s_new)))

summ <- function(x, nm) {
  v <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  tibble(model = nm, slope = v[2], lo95 = v[1], hi95 = v[3], n_groups = NA_integer_)
}
TT_bCovT <- as.matrix(s_new$TT_bCovT)
par_tbl <- bind_rows(
  summ(s_new$TT_bCov, "03g posterior TT_bCov (mainstem)"),
  map(seq_len(ncol(TT_bCovT)), \(k) summ(TT_bCovT[, k], paste0("03g posterior TT_bCovT[", k, "] (trib)"))) |> list_rbind(),
  summ(s_new$T_bSz, "03g posterior T_bSz (size)")
)

tabA <- bind_rows(par_tbl[1, ], slopes, site_slopes, par_tbl[-1, ]) |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(tabA, file.path(tab_dir, "A_discharge_slopes_data_vs_model.csv"))

# ============================================================
# B. Model-predicted vs observed median travel time by release group
# ============================================================
# TT_reach[draw, fish, reach] = predicted median travel time in each reach,
# with that fish's release-group effect. Reaches 1 to 3 end at Sacramento.
pred_draws <- s_new$TT_reach[, , 1] + s_new$TT_reach[, , 2] + s_new$TT_reach[, , 3]

pred_grp <- split(fish$row, fish$rg) |>
  imap(\(i, g) {
    d <- apply(pred_draws[, i, drop = FALSE], 1, median)   # group median, each draw
    qt(d) |> rename(pred_lo = lo, pred_med = med, pred_hi = hi) |> mutate(rg = as.integer(g))
  }) |>
  list_rbind()
rm(pred_draws); invisible(gc())

tabB <- grp |>
  left_join(pred_grp, by = "rg") |>
  mutate(log_resid = log(obs_med / pred_med)) |>
  select(rg, study, site, year, n, flow, fl, obs_med, pred_med, pred_lo, pred_hi, log_resid) |>
  arrange(site, flow)
write_csv(tabB |> mutate(across(where(is.double), \(x) round(x, 3))),
          file.path(tab_dir, "B_obs_vs_pred_by_release_group.csv"))

pB1 <- ggplot(tabB, aes(pred_med, obs_med)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(xmin = pred_lo, xmax = pred_hi), colour = "#56B4E9") +
  geom_point(aes(size = n), shape = 21, fill = "#0072B2", colour = "black") +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ site) +
  labs(x = "Predicted median travel time, release to Sacramento (days, log)",
       y = "Observed median (days, log)", size = "Fish in group",
       caption = "Dashed: 1:1. Bars: 80% interval of the predicted group median.")
ggsave(file.path(fig_dir, "B1_obs_vs_pred_by_group.png"), pB1, width = 10, height = 7, dpi = 300)

# Residuals vs discharge and size: a trend means the model misses that covariate
pB2 <- tabB |>
  transmute(n, log_resid, `Discharge (kcfs)` = flow / 1000, `Mean fork length (mm)` = fl) |>
  pivot_longer(-c(n, log_resid), names_to = "covariate", values_to = "x") |>
  ggplot(aes(x, log_resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(size = n), shape = 21, fill = "#0072B2", colour = "black") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#D55E00") +
  facet_wrap(~ covariate, scales = "free_x") +
  labs(x = NULL, y = "log(observed / predicted) group median", size = "Fish in group")
ggsave(file.path(fig_dir, "B2_residuals_vs_discharge_and_size.png"), pB2, width = 10, height = 4.5, dpi = 300)

# ============================================================
# C. RST to Delta entry travel-time curves, new vs old stm
# ============================================================
# Mainstem RST order from 02 (dist_rstwoodson): reach 1 = RST to Woodson,
# then Woodson-Butte and Butte-Sacramento station (= Delta entry)
rst_names <- c("Upper Clear Ck", "Battle Ck", "Mill Ck", "Deer Ck", "Redding", "Red Bluff")

tt_curves <- function(s, label) {
  stopifnot(dim(s$TTForecastSz_rst)[4] == length(rst_names))
  arr <- s$TTForecastSz_rst                       # [draw, size, discharge, RST]
  isz <- which.min(abs(s$Xsz))                    # grid point nearest mean size
  ix  <- which.min(abs(s$Xvec))                   # grid point nearest mean discharge
  nx  <- dim(arr)[3]
  by_flow <- expand_grid(rst = seq_len(dim(arr)[4]), j = seq_len(nx)) |>
    mutate(q = map2(rst, j, \(r, jj) qt(arr[, isz, jj, r])),
           x = s$Xvec[j] * s$sdx + s$mux, curve = "discharge",
           at = case_when(j == 1 ~ "min discharge", j == ix ~ "mean discharge",
                          j == nx ~ "max discharge", .default = NA_character_))
  by_size <- expand_grid(rst = seq_len(dim(arr)[4]), k = seq_len(dim(arr)[2])) |>
    mutate(q = map2(rst, k, \(r, kk) qt(arr[, kk, ix, r])),
           x = s$Xsz[k] * s$sdsz + s$musz, curve = "size")
  bind_rows(by_flow, by_size) |>
    unnest(q) |>
    select(-j, -k) |>
    mutate(rst = factor(rst, levels = seq_along(rst_names), labels = rst_names),   # 02's RST order
           stm = label,
           mean_size = s$Xsz[isz] * s$sdsz + s$musz,
           mean_flow = s$Xvec[ix] * s$sdx + s$mux)
}

curves_new <- tt_curves(s_new, "new")

# ============================================================
# D. Fish-level spread around the median (needs Pro_sd in the stm)
# ============================================================
tabD <- if (!is.null(s_new$Pro_sd)) {
  isz <- which.min(abs(s_new$Xsz)); ix <- which.min(abs(s_new$Xvec))
  map(seq_len(dim(s_new$TTForecastSz_rst)[4]), \(r) {
    med <- s_new$TTForecastSz_rst[, isz, ix, r]
    tibble(rst = rst_names[r],
           median_tt   = median(med),
           fish_10pct  = median(med * exp(qnorm(0.1) * s_new$Pro_sd)),
           fish_90pct  = median(med * exp(qnorm(0.9) * s_new$Pro_sd)),
           Pro_sd      = median(s_new$Pro_sd))
  }) |> list_rbind() |> mutate(across(where(is.double), \(x) round(x, 2)))
} else {
  message("Pro_sd is not in the stm yet: run 05d and 11, then rerun this script.")
  NULL
}
if (!is.null(tabD)) write_csv(tabD, file.path(tab_dir, "D_fish_level_spread_at_mean.csv"))

rm(s_new); invisible(gc())

# Old stm (previous team stm; its discharge grid may be on a different gauge)
e <- new.env(); load(stm_old_file, envir = e); curves_old <- tt_curves(e$stm, "old"); rm(e); invisible(gc())

curves <- bind_rows(curves_new, curves_old)

plot_curves <- function(which_curve, xlab, xdiv = 1, vline = NULL) {
  p <- curves |>
    filter(curve == which_curve) |>
    ggplot(aes(x / xdiv, med, colour = stm, fill = stm, linetype = stm)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ rst, scales = "free_y") +
    scale_colour_manual(values = stm_col, labels = stm_lab) +
    scale_fill_manual(values = stm_col, labels = stm_lab) +
    scale_linetype_manual(values = c(new = "solid", old = "longdash"), labels = stm_lab) +
    labs(x = xlab, y = "RST to Delta entry travel time (days)", colour = NULL, fill = NULL, linetype = NULL) +
    theme(legend.position = "bottom")
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline, linetype = "dotted")
  p
}

pC1 <- plot_curves("discharge", "Discharge (kcfs), at each stm's mean fork length", xdiv = 1000)
ggsave(file.path(fig_dir, "C1_rst_delta_tt_vs_discharge.png"), pC1, width = 10, height = 7, dpi = 300)

pC2 <- plot_curves("size", "Fork length (mm), at each stm's mean discharge",
                   vline = min(FL, na.rm = TRUE)) +
  labs(caption = "Dotted: smallest tagged fish in the telemetry data. Left of it is extrapolation.")
ggsave(file.path(fig_dir, "C2_rst_delta_tt_vs_size.png"), pC2, width = 10, height = 7, dpi = 300)

tabC <- curves |>
  filter(curve == "discharge", !is.na(at)) |>
  transmute(stm, rst, at, discharge_cfs = round(x), mean_size_mm = round(mean_size, 1),
            tt_med = round(med, 1), tt_lo80 = round(lo, 1), tt_hi80 = round(hi, 1)) |>
  arrange(rst, at, stm)
write_csv(tabC, file.path(tab_dir, "C_rst_delta_tt_new_vs_old.csv"))

# ============================================================
# Console summary
# ============================================================
cat("\n=== A. Discharge slope on log travel time (TT_bCov scale) ===\n"); print(as.data.frame(tabA))
cat("\n=== B. Worst-fitting release groups (|log resid| largest) ===\n")
print(as.data.frame(tabB |> slice_max(abs(log_resid), n = 10) |> mutate(across(where(is.double), \(x) round(x, 2)))))
cat("\n=== C. RST to Delta travel time, new vs old ===\n"); print(as.data.frame(tabC))
if (!is.null(tabD)) { cat("\n=== D. Fish-level spread at mean size and discharge ===\n"); print(as.data.frame(tabD)) }
cat("\nFigures:", fig_dir, "\nTables: ", tab_dir, "\n")
