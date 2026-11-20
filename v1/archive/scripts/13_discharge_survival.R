# 13_discharge_survival.R
# Survival under the dry / average / wet scenario seasons from 12, from the corrected telemetry stm.
#
# Survival = release to the Sacramento station (Delta Entry 1st Bridge), typical release group
# (no release-group effect, the _nore arrays), mean telemetry fork length, as defined in 05. These
# curves show the discharge effect alone; production survival should use the arrays WITH the
# release-group effect. Each posterior draw is interpolated on the logit scale across the stm
# discharge grid; discharge outside the grid is clamped to the grid end and flagged.
#
# Scenario survival = survival in the scenario season(s) from 12 (within each draw, the mean over
#   the scenario's seasons when seasons_per_scen > 1).
# Expected survival = within each draw, the weighted mean over all seasons of survival at that
#   season's actual monthly max (the mixture; equal weights unless season_weights is set). This is
#   the central value to hand on. The same draw is used for every gauge, month and season.
#
# Inputs:  results/stm_MaxFlowFL_KLdelta.Rdata (from 05)
#          results/tables/12_monthly_max_by_season.csv, results/tables/12_scenario_seasons.csv
# Outputs: results/tables/13_stm_grid_ranges.csv          discharge range covered by each stm grid
#          results/tables/13_survival_scenarios.csv        survival in the dry / average / wet seasons
#          results/tables/13_survival_by_season_month.csv  survival at each season's observed discharge
#          results/tables/13_survival_expected.csv         expected survival (mixture over seasons)
#          results/tables/13_clamp_share.csv               share of seasons clamped, by gauge x month
#          results/figures/13_survival_scenarios.png

library(tidyverse)
library(here)

# ---- Settings ----------------------------------------------------------------------------------
stm_file      <- here("results", "stm_MaxFlowFL_KLdelta.Rdata")
season_months <- c(11, 12, 1, 2, 3, 4, 5)
month_labels  <- month.abb[season_months]
scen_levels   <- c("dry", "average", "wet")

# Season weights for the expected-survival mixture: NULL = equal weights. Later this can come from
# a DWR forecast, e.g. tibble(season = ..., weight = ...); weights are rescaled to sum to 1.
season_weights <- NULL

# stm arrays: SurvForecast_nore is draws x 25 over Xvec (z scale, mux / sdx). TribSurvForecast_nore
# is draws x 25 x 2 over the columns of XvecT (Butte first, Feather second).
surv_spec <- tribble(
  ~gauge,       ~array,                  ~grid,   ~grid_col,   ~mu,    ~sd,    ~trib_index,
  "Sacramento", "SurvForecast_nore",     "Xvec",  NA_integer_, "mux",  "sdx",  NA_integer_,
  "Butte",      "TribSurvForecast_nore", "XvecT", 1L,          "muxB", "sdxB", 1L,
  "Feather",    "TribSurvForecast_nore", "XvecT", 2L,          "muxF", "sdxF", 2L
)
gauge_labels <- c(
  Sacramento = "Sacramento (release to Sacramento station)",
  Butte      = "Butte Creek (release to Sacramento station)",
  Feather    = "Feather River (release to Sacramento station)"
)

pal  <- c(dry = "#D55E00", average = "#000000", wet = "#0072B2")  # Okabe-Ito
ltys <- c(dry = "dashed", average = "solid", wet = "longdash")
col_expected <- "#CC79A7"                                       # Okabe-Ito reddish purple

# ---- stm ---------------------------------------------------------------------------------------
load(stm_file)
if (!exists("stm")) stop(stm_file, " did not create an object called 'stm'")

get_grid <- function(grid, grid_col) {
  g <- stm[[grid]]
  if (is.null(g)) stop("stm has no '", grid, "'")
  if (!is.na(grid_col)) g <- as.matrix(g)[, grid_col]
  as.numeric(g)
}

get_surv_matrix <- function(array, g, trib_index) {
  a <- stm[[array]]
  if (is.null(a)) stop("stm has no '", array, "'")
  d <- dim(a)
  if (length(d) == 3) {
    if (is.na(trib_index)) stop(array, " is 3-dimensional; set trib_index in surv_spec")
    tdim <- setdiff(seq_along(d), c(which(d == length(g))[1], which.max(d)))
    if (length(tdim) != 1) {
      stop("Cannot tell which dimension of ", array, " is the tributary (dims ",
           paste(d, collapse = " x "), ")")
    }
    idx <- rep(list(TRUE), 3)
    idx[[tdim]] <- trib_index
    a <- do.call(`[`, c(list(a), idx, list(drop = TRUE)))
  } else if (length(d) != 2) {
    stop(array, " has ", length(d), " dimensions; expected 2 or 3")
  }
  if (ncol(a) != length(g)) a <- t(a)
  stopifnot(ncol(a) == length(g))
  a
}

# Interpolate each draw on the logit scale; targets outside the grid are clamped to the grid ends
interp_logit <- function(m, gz, x) {
  lm <- qlogis(pmin(pmax(m, 1e-9), 1 - 1e-9))
  xc <- pmin(pmax(x, min(gz)), max(gz))
  i  <- findInterval(xc, gz, all.inside = TRUE)
  w  <- (xc - gz[i]) / (gz[i + 1] - gz[i])
  plogis(sweep(lm[, i, drop = FALSE], 2, 1 - w, `*`) + sweep(lm[, i + 1, drop = FALSE], 2, w, `*`))
}

draw_summary <- function(s) {
  q <- apply(s, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
  tibble(surv_mean = colMeans(s), surv_q025 = q[1, ], surv_q10 = q[2, ], surv_med = q[3, ],
         surv_q90 = q[4, ], surv_q975 = q[5, ])
}

# Within-draw weighted mean over a set of columns, summarized across draws
mix_summary <- function(s, cols, w) draw_summary(matrix(s[, cols, drop = FALSE] %*% (w / sum(w)), ncol = 1))

# ---- Inputs from 12 ----------------------------------------------------------------------------
mm <- read_csv(here("results", "tables", "12_monthly_max_by_season.csv"), show_col_types = FALSE) |>
  mutate(scenario = factor(scenario, levels = scen_levels))
scen_seasons <- read_csv(here("results", "tables", "12_scenario_seasons.csv"), show_col_types = FALSE)

weights <- (season_weights %||% tibble(season = unique(mm$season), weight = 1)) |>
  filter(season %in% mm$season)
if (nrow(weights) != n_distinct(mm$season)) stop("season_weights must cover every season in 12")
mm <- left_join(mm, weights, by = "season")

# ---- Survival by gauge -------------------------------------------------------------------------
survival_for_gauge <- function(gauge, array, grid, grid_col, mu, sd, trib_index) {
  gz  <- get_grid(grid, grid_col)
  m   <- get_surv_matrix(array, gz, trib_index)
  ord <- order(gz)
  gz  <- gz[ord]
  m   <- m[, ord, drop = FALSE]

  on_cfs   <- max(abs(gz)) > 50                       # grid stored in cfs rather than z
  to_x     <- \(cfs) if (on_cfs) cfs else (cfs - stm[[mu]]) / stm[[sd]]
  grid_cfs <- if (on_cfs) range(gz) else range(gz) * stm[[sd]] + stm[[mu]]

  mm_g <- filter(mm, gauge == .env$gauge) |>
    mutate(col = row_number(),
           outside_fit = discharge < grid_cfs[1] | discharge > grid_cfs[2])
  s_mm <- interp_logit(m, gz, to_x(mm_g$discharge))   # draws x season-months

  by_month <- \(d, ...) d |>
    summarise(cols = list(col), w = list(weight), n_seasons = n(),
              share_clamped = mean(outside_fit), .by = c(...)) |>
    mutate(stats = map2(cols, w, \(cc, ww) mix_summary(s_mm, cc, ww))) |>
    select(-cols, -w) |>
    unnest(stats) |>
    mutate(gauge = .env$gauge, .before = 1)

  list(
    grid_ranges   = tibble(gauge = gauge, n_draws = nrow(m), grid_lo_cfs = grid_cfs[1],
                           grid_hi_cfs = grid_cfs[2]),
    scenarios     = mm_g |> filter(!is.na(scenario)) |> mutate(weight = 1) |>
      by_month(scenario, month, mpos, month_lab),
    season_points = bind_cols(select(mm_g, -col), draw_summary(s_mm)),
    expected      = by_month(mm_g, month, mpos, month_lab)
  )
}

res <- surv_spec |>
  pmap(survival_for_gauge) |>
  list_transpose() |>
  map(bind_rows)

cat("\nDischarge range of each stm grid (compare with fit_range in 12):\n")
print(res$grid_ranges)

clamp_share <- res$expected |>
  select(gauge, month_lab, mpos, n_seasons, share_clamped) |>
  arrange(gauge, mpos)

write_csv(res$grid_ranges,   here("results", "tables", "13_stm_grid_ranges.csv"))
write_csv(res$scenarios,     here("results", "tables", "13_survival_scenarios.csv"))
write_csv(res$season_points, here("results", "tables", "13_survival_by_season_month.csv"))
write_csv(res$expected,      here("results", "tables", "13_survival_expected.csv"))
write_csv(clamp_share,       here("results", "tables", "13_clamp_share.csv"))

cat("\nScenario seasons (November to May): ",
    paste0(scen_seasons$scenario, " = ", scen_seasons$season - 1, "-", substr(scen_seasons$season, 3, 4),
           collapse = ", "), "\n", sep = "")
cat("Survival by scenario (posterior median) and expected survival (posterior mean of the mixture):\n")
bind_rows(
  select(res$scenarios, gauge, row = scenario, mpos, month_lab, value = surv_med),
  mutate(select(res$expected, gauge, mpos, month_lab, value = surv_mean), row = "expected")
) |>
  mutate(value = round(value, 3), row = factor(row, levels = c(scen_levels, "expected"))) |>
  arrange(gauge, row, mpos) |>
  select(-mpos) |>
  pivot_wider(names_from = month_lab, values_from = value) |>
  print(n = Inf)

cat("\nShare of seasons with discharge outside the fitted range (survival clamped):\n")
clamp_share |>
  mutate(share_clamped = round(share_clamped, 2)) |>
  select(gauge, month_lab, share_clamped) |>
  pivot_wider(names_from = month_lab, values_from = share_clamped) |>
  print()

# ---- Figure --------------------------------------------------------------------------------------
scen_labs <- scen_seasons |>
  mutate(season_label = paste0(season - 1, "-", substr(season, 3, 4))) |>   # 2022 -> "2021-22"
  summarise(lab = paste0(first(scenario), " (", paste(season_label, collapse = ", "), ")"), .by = scenario) |>
  (\(d) set_names(d$lab, d$scenario))()

prep <- \(d) d |>
  mutate(gauge = factor(gauge, levels = surv_spec$gauge),
         across(any_of("scenario"), \(x) factor(x, levels = scen_levels)),
         across(any_of("outside_fit"),
                \(x) if_else(x, "outside fitted discharge range (clamped)", "inside fitted discharge range")))

p13 <- ggplot() +
  geom_ribbon(data = prep(res$scenarios),
              aes(mpos, ymin = surv_q10, ymax = surv_q90, fill = scenario), alpha = 0.15) +
  geom_point(data = prep(res$season_points),
             aes(mpos, surv_med, shape = outside_fit),
             position = position_jitter(width = 0.15, height = 0, seed = 1),
             colour = "grey35", alpha = 0.6, size = 1.6) +
  geom_line(data = prep(res$scenarios),
            aes(mpos, surv_med, colour = scenario, linetype = scenario), linewidth = 1) +
  geom_line(data = prep(res$expected),
            aes(mpos, surv_mean), colour = col_expected, linetype = "dotted", linewidth = 1.2) +
  facet_wrap(~gauge, labeller = as_labeller(gauge_labels)) +
  scale_x_continuous(breaks = seq_along(month_labels), labels = month_labels) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = pal,    breaks = scen_levels, labels = scen_labs, name = NULL) +
  scale_fill_manual(values = pal,      breaks = scen_levels, labels = scen_labs, name = NULL) +
  scale_linetype_manual(values = ltys, breaks = scen_levels, labels = scen_labs, name = NULL) +
  scale_shape_manual(values = c("inside fitted discharge range" = 16,
                                "outside fitted discharge range (clamped)" = 1), name = NULL) +
  labs(x = NULL, y = "Survival, release to Sacramento station",
       title = "Survival under the discharge scenarios, by month",
       subtitle = "Scenarios are real seasons; typical release group (no release-group effect), mean telemetry fork length",
       caption = paste(
         "Lines: posterior median in the dry, average and wet seasons from 12. Band: 80% credible interval (parameter uncertainty only).",
         "Dots: posterior median at each season's observed monthly max; open dots fall outside the discharge range of the telemetry fit (clamped).",
         "Dotted purple: expected survival, the within-draw mean over all seasons (the mixture). Survival is curved in discharge, so this differs from the average season.",
         sep = "\n")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", plot.caption = element_text(hjust = 0))

ggsave(here("results", "figures", "13_survival_scenarios.png"), p13, width = 13, height = 5.5, dpi = 300)
message("Wrote results/figures/13_survival_scenarios.png")
