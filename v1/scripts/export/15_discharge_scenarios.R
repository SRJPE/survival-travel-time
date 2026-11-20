# ---- srjpe_survival: export step 15, dry / average / wet discharge seasons
# ---- made from scripts/12_discharge_classes.R (its header says 12_discharge_scenarios.R) by setup/build_from_old_project.R (exact edits listed there)

# 12_discharge_scenarios.R
# Discharge scenarios for the survival covariate, built from whole historical seasons.
#
# What: monthly max discharge (the covariate the telemetry survival model was fitted on) at the three
#       survival gauges, November to May, for every season in SRJPEdata. Each season gets one wetness
#       index; dry / average / wet are the real seasons nearest the 16th / 50th / 84th percentiles of
#       that index. A scenario therefore uses one actual season's discharge in every month and at
#       every gauge, so months stay consistent within the year and the three rivers stay consistent
#       with each other. (The previous version took the 16th / 50th / 84th percentile separately for
#       each gauge and month, which strings together a year that never happened.)
# Season: November (year - 1) to May (year). Stored as the May year ('season', = water year), shown
#         as both years ('season_label': season 2022 = "2021-22" = November 2021 to May 2022).
#
# Inputs:  SRJPEdata (pinned commit): forecast_covariates, flow_data
# Outputs: results/tables/12_year_alignment_check.csv   is monthly_max_flow 'year' calendar or water year?
#          results/tables/12_monthly_max_by_season.csv  one row per gauge x season x month
#          results/tables/12_season_index.csv           wetness index and percentile rank per season
#          results/tables/12_scenario_seasons.csv       which seasons define dry / average / wet
#          results/tables/12_discharge_scenarios.csv    scenario discharge by gauge x month (x season)
#          results/figures/12_discharge_scenarios.png

library(tidyverse)
library(here)
library(SRJPEdata)

# ---- Settings ----------------------------------------------------------------------------------
pinned_sha        <- "d3584db2a7230a85d3fdb06b5df1b9f794c3a197"
season_months     <- c(11, 12, 1, 2, 3, 4, 5)
month_labels      <- month.abb[season_months]
scenario_probs    <- c(dry = 0.16, average = 0.50, wet = 0.84)
seasons_per_scen  <- 1               # seasons nearest each percentile; 3 = average over 3 seasons
season_range      <- NULL            # e.g. c(1992, 2026); NULL = all seasons complete at every gauge
telemetry_seasons <- 2013:2024       # seasons with telemetry releases in the survival fit (dots filled)

# Wetness index used to rank seasons:
#   "sacramento": Nov-May mean of log monthly max at the Sacramento gauge
#   "all_gauges": Nov-May mean over the three gauges of log monthly max, each gauge standardized
#                 (z of log discharge by gauge and month), so all three rivers count equally
index_source <- "all_gauges"

# Survival gauges, same series as the telemetry fit (03g)
gauges <- tribble(
  ~gauge,       ~stream,            ~fc_site_group,    ~fd_site_group,
  "Sacramento", "sacramento river", "knights landing", "knights landing",
  "Butte",      "butte creek",      NA,                "butte creek",
  "Feather",    "feather river",    NA,                NA  # SRJPEdata's Feather series (max of available gauges)
)
gauge_labels <- c(
  Sacramento = "Sacramento (USGS 11390500, Wilkins Slough; 'knights landing')",
  Butte      = "Butte Creek (USGS 11390000)",
  Feather    = "Feather River (SRJPEdata composite of Feather gauges)"
)

# Discharge range of the telemetry releases in the survival fit (05 forecast grids,
# stm_MaxFlowFL_KLdelta). Survival outside it is clamped. 13 prints the exact grid ends from the stm.
fit_range <- tribble(
  ~gauge,       ~lo,  ~hi,
  "Sacramento", 4480, 24100,
  "Butte",      168,  3160,
  "Feather",    983,  30325
)

scen_levels <- names(scenario_probs)
season_lab  <- \(s) paste0(s - 1, "-", substr(s, 3, 4))   # 2022 -> "2021-22"
pal  <- c(dry = "#D55E00", average = "#000000", wet = "#0072B2")  # Okabe-Ito
ltys <- c(dry = "dashed", average = "solid", wet = "longdash")

dir.create(here("results", "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ---- Data --------------------------------------------------------------------------------------
sha <- packageDescription("SRJPEdata")$RemoteSha %||% "unknown"
if (!identical(sha, pinned_sha)) warning("SRJPEdata commit is ", sha, ", not the pinned ", pinned_sha)

load_obj <- function(nm) {
  e <- new.env()
  data(list = nm, package = "SRJPEdata", envir = e)
  e[[nm]]
}
fc <- load_obj("forecast_covariates")
fd <- load_obj("flow_data")

mm_raw <- fc |>
  filter(name == "monthly_max_flow") |>
  inner_join(select(gauges, gauge, stream, fc_site_group), by = "stream") |>
  filter(is.na(fc_site_group) | site_group == fc_site_group) |>
  select(gauge, year, month, discharge = value)

missing_g <- setdiff(gauges$gauge, mm_raw$gauge)
if (length(missing_g) > 0) stop("No monthly_max_flow rows for: ", paste(missing_g, collapse = ", "))
dups <- mm_raw |> count(gauge, year, month) |> filter(n > 1)
if (nrow(dups) > 0) {
  print(dups)
  stop("More than one monthly max per gauge x year x month (see above); check the site_group filter")
}

# ---- Check: is 'year' in monthly_max_flow the calendar year or the water year? -------------------
# Only Nov and Dec differ between the two. Rebuild an approximate monthly max from the weekly max in
# flow_data (week w of year y assigned to the month of its middle day) and see which alignment matches.
wk_month <- fd |>
  filter(statistic == "max", !is.na(value)) |>
  inner_join(gauges |> filter(!is.na(fd_site_group)) |> select(gauge, fd_site_group),
             by = c("site_group" = "fd_site_group")) |>
  mutate(month = month(make_date(year, 1, 1) + days(7 * (week - 1) + 3))) |>
  summarise(wk_max = max(value, na.rm = TRUE), .by = c(gauge, year, month))

align <- mm_raw |>
  filter(month %in% c(11, 12)) |>
  mutate(calendar = year, water_year = year - 1) |>
  pivot_longer(c(calendar, water_year), names_to = "if_year_is", values_to = "cal_year") |>
  inner_join(wk_month, by = c("gauge", "cal_year" = "year", "month")) |>
  filter(discharge > 0, wk_max > 0) |>
  summarise(n = n(), median_abs_log_ratio = median(abs(log(discharge / wk_max))),
            .by = c(gauge, if_year_is))
write_csv(align, here("results", "tables", "12_year_alignment_check.csv"))
print(align)

year_is <- align |>
  summarise(score = mean(median_abs_log_ratio), .by = if_year_is) |>
  slice_min(score, n = 1) |>
  pull(if_year_is)
message("monthly_max_flow 'year' matches the ", year_is, " year; seasons are built on that basis")

# ---- Monthly max by season -----------------------------------------------------------------------
mm_all <- mm_raw |>
  filter(month %in% season_months) |>
  mutate(season    = year + if (year_is == "calendar") as.integer(month >= 11) else 0L,
         mpos      = match(month, season_months),
         month_lab = month_labels[mpos],
         season_label = season_lab(season)) |>
  mutate(complete = n() == length(season_months) & all(!is.na(discharge)),
         .by = c(gauge, season))

# Seasons used everywhere below: complete at all three gauges (and inside season_range, if set),
# so the scenarios and the expected-survival mixture in 13 use the same seasons for every river.
common_seasons <- mm_all |>
  filter(complete) |>
  distinct(gauge, season) |>
  count(season) |>
  filter(n == nrow(gauges)) |>
  pull(season) |>
  (\(s) if (is.null(season_range)) s else s[between(s, season_range[1], season_range[2])])() |>
  sort()
message("Seasons complete at all three gauges: ", length(common_seasons), " (",
        season_lab(min(common_seasons)), " to ", season_lab(max(common_seasons)), ")")

mm <- mm_all |>
  filter(season %in% common_seasons) |>
  mutate(telemetry_season = season %in% telemetry_seasons)

# ---- Wetness index and scenario seasons ----------------------------------------------------------
index_by_gauge <- mm |>
  mutate(logq = log(pmax(discharge, 1)),
         z    = (logq - mean(logq)) / sd(logq), .by = c(gauge, month)) |>
  summarise(index = mean(z), .by = c(gauge, season))

cat("\nRank correlation of the season index between gauges (do the rivers agree on wet vs dry?):\n")
index_by_gauge |>
  pivot_wider(names_from = gauge, values_from = index) |>
  select(-season) |>
  cor(method = "spearman") |>
  round(2) |>
  print()

season_index <- (if (index_source == "sacramento")
  filter(index_by_gauge, gauge == "Sacramento") else index_by_gauge) |>
  summarise(index = mean(index), .by = season) |>
  mutate(pct = percent_rank(index), season_label = season_lab(season), .after = season) |>
  arrange(index)

scenario_seasons <- scenario_probs |>
  map(\(p) season_index |>
        mutate(target = p, dist = abs(pct - p)) |>
        slice_min(dist, n = seasons_per_scen, with_ties = FALSE)) |>
  list_rbind(names_to = "scenario") |>
  mutate(scenario = factor(scenario, levels = scen_levels)) |>
  arrange(scenario, season)

if (anyDuplicated(scenario_seasons$season)) {
  warning("A season is in more than one scenario; reduce seasons_per_scen")
}

season_index <- season_index |>
  left_join(select(scenario_seasons, season, scenario), by = "season")
write_csv(season_index,     here("results", "tables", "12_season_index.csv"))
write_csv(scenario_seasons, here("results", "tables", "12_scenario_seasons.csv"))

cat("\nScenario seasons (index = ", index_source, "):\n", sep = "")
print(select(scenario_seasons, scenario, season_label, index, pct))

mm <- mm |> left_join(select(scenario_seasons, season, scenario), by = "season")
write_csv(mm, here("results", "tables", "12_monthly_max_by_season.csv"))

# Scenario discharge: the scenario seasons' actual values. With seasons_per_scen > 1 the per-season
# rows are kept (13 averages survival over them); 'discharge_med' is only for display.
scen_discharge <- mm |>
  filter(!is.na(scenario)) |>
  mutate(discharge_med = median(discharge), n_seasons = n(), .by = c(gauge, scenario, month)) |>
  mutate(outside_fit = discharge < fit_range$lo[match(gauge, fit_range$gauge)] |
                       discharge > fit_range$hi[match(gauge, fit_range$gauge)]) |>
  arrange(gauge, scenario, season, mpos)
write_csv(scen_discharge, here("results", "tables", "12_discharge_scenarios.csv"))

cat("\nScenario discharge (cfs; median over the scenario seasons):\n")
scen_discharge |>
  distinct(gauge, scenario, month_lab, mpos, discharge_med) |>
  arrange(gauge, scenario, mpos) |>
  select(-mpos) |>
  pivot_wider(names_from = month_lab, values_from = discharge_med) |>
  print(n = Inf)

cat("\nShare of scenario months outside the fitted discharge range:\n")
scen_discharge |>
  summarise(share_outside = round(mean(outside_fit), 2), .by = c(gauge, scenario)) |>
  pivot_wider(names_from = scenario, values_from = share_outside) |>
  print()

# ---- Figure --------------------------------------------------------------------------------------
scen_labs <- scenario_seasons |>
  summarise(lab = paste0(first(scenario), " (", paste(season_label, collapse = ", "), ")"),
            .by = scenario) |>
  (\(d) set_names(d$lab, as.character(d$scenario)))()

as_gauge_factor <- \(d) mutate(d, gauge = factor(gauge, levels = gauges$gauge))

p12 <- ggplot() +
  geom_point(data = as_gauge_factor(mm),
             aes(mpos, discharge, shape = telemetry_season),
             position = position_jitter(width = 0.15, height = 0, seed = 1),
             colour = "grey35", alpha = 0.6, size = 1.6) +
  geom_line(data = as_gauge_factor(scen_discharge),
            aes(mpos, discharge, colour = scenario, linetype = scenario,
                group = interaction(scenario, season)),
            linewidth = 1) +
  geom_hline(data = fit_range |> pivot_longer(c(lo, hi), values_to = "cfs") |> as_gauge_factor(),
             aes(yintercept = cfs), linetype = "dotted", colour = "grey40") +
  facet_wrap(~gauge, scales = "free_y", labeller = as_labeller(gauge_labels)) +
  scale_x_continuous(breaks = seq_along(month_labels), labels = month_labels) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = pal,    breaks = scen_levels, labels = scen_labs, name = NULL) +
  scale_linetype_manual(values = ltys, breaks = scen_levels, labels = scen_labs, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "other seasons",
                                `TRUE`  = paste0("telemetry seasons (", season_lab(min(telemetry_seasons)),
                                                 " to ", season_lab(max(telemetry_seasons)), ")")),
                     name = NULL) +
  labs(x = NULL, y = "Monthly max discharge (cfs, log scale)",
       title = "Discharge scenarios for the survival covariate, November to May",
       subtitle = paste0("Seasons ", season_lab(min(common_seasons)), " to ", season_lab(max(common_seasons)),
                         " (November to May); scenarios are real seasons; SRJPEdata ", substr(sha, 1, 7)),
       caption = paste(
         paste0("Lines: the seasons nearest the 16th, 50th and 84th percentiles of a season wetness index (",
                if (index_source == "sacramento") "Sacramento gauge" else "all three gauges",
                ", mean of log monthly max, Nov to May)."),
         "Each scenario is one real season at every gauge and in every month. A season runs November to May: 2021-22 = November 2021 to May 2022.",
         "Dots: each season's monthly max (filled: seasons with telemetry releases in the survival fit).",
         "Dotted grey: discharge range of the telemetry releases in the survival fit; outside it survival is clamped.",
         sep = "\n")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", plot.caption = element_text(hjust = 0))

ggsave(here("results", "figures", "12_discharge_scenarios.png"), p12, width = 13, height = 5.5, dpi = 300)
message("Wrote results/figures/12_discharge_scenarios.png")
