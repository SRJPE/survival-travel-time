# 12_discharge_scenarios.R
# Discharge scenarios for the survival covariate.
#
# What: monthly max discharge, the covariate the telemetry survival model was fitted on, at the three
#       survival gauges, November to May, for every season in SRJPEdata. Dry / average / wet are the
#       16th / 50th / 84th percentiles of monthly max discharge for each gauge and month (the discharge
#       memo's -1 / 0 / +1 SD, without negative cfs at Butte and Feather). No year classification:
#       expected survival in 13 uses each season's actual monthly values instead.
# Season: November (year - 1) to May (year), labelled by the May year (= water year).
#
# Inputs:  SRJPEdata (pinned commit): forecast_covariates, flow_data
# Outputs: results/tables/12_year_alignment_check.csv   is monthly_max_flow 'year' calendar or water year?
#          results/tables/12_monthly_max_by_season.csv  one row per gauge x season x month
#          results/tables/12_discharge_scenarios.csv    dry / average / wet discharge by gauge x month
#          results/figures/12_discharge_scenarios.png

library(tidyverse)
library(here)
library(SRJPEdata)

# ---- Settings ----------------------------------------------------------------------------------
pinned_sha     <- "d3584db2a7230a85d3fdb06b5df1b9f794c3a197"
season_months  <- c(11, 12, 1, 2, 3, 4, 5)
month_labels   <- month.abb[season_months]
scenario_probs <- c(dry = 0.16, average = 0.50, wet = 0.84)
season_range   <- NULL           # e.g. c(1991, 2026) to restrict the record; NULL = all complete seasons
telemetry_seasons <- 2013:2024   # seasons with telemetry releases in the survival fit (dots filled)

# Survival gauges, same series as the telemetry fit (03g)
gauges <- tribble(
  ~gauge,       ~stream,            ~fc_site_group,    ~fd_site_group,
  "Sacramento", "sacramento river", "knights landing", "knights landing",
  "Butte",      "butte creek",      NA,                "butte creek",
  "Feather",    "feather river",    NA,                NA  # composite: max of the available Feather gauges
)
gauge_labels <- c(
  Sacramento = "Sacramento ('knights landing' = USGS 11390500, Wilkins Slough)",
  Butte      = "Butte Creek (USGS 11390000)",
  Feather    = "Feather River (composite of Feather gauges)"
)

# Monthly max discharge of the telemetry release months in the survival fit (05 forecast grids,
# stm_MaxFlowFL_KLdelta). Survival outside this range is extrapolated. 13 prints the exact grid ends.
fit_range <- tribble(
  ~gauge,       ~lo,  ~hi,
  "Sacramento", 4480, 24100,
  "Butte",      168,  3160,
  "Feather",    983,  30325
)

pal  <- c(dry = "#D55E00", average = "#000000", wet = "#0072B2")  # Okabe-Ito
ltys <- c(dry = "dashed", average = "solid", wet = "longdash")

dir.create(here("results", "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ---- Data --------------------------------------------------------------------------------------
sha <- packageDescription("SRJPEdata")$RemoteSha
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
  filter(statistic == "max") |>
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
mm <- mm_raw |>
  filter(month %in% season_months) |>
  mutate(season    = year + if (year_is == "calendar") as.integer(month >= 11) else 0L,
         mpos      = match(month, season_months),
         month_lab = month_labels[mpos]) |>
  mutate(complete = n() == length(season_months), .by = c(gauge, season)) |>
  mutate(in_range = if (is.null(season_range)) TRUE else between(season, season_range[1], season_range[2]),
         telemetry_season = season %in% telemetry_seasons)
write_csv(mm, here("results", "tables", "12_monthly_max_by_season.csv"))

use <- filter(mm, complete, in_range)
cat("\nComplete seasons used, by gauge:\n")
use |> summarise(n_seasons = n_distinct(season), first = min(season), last = max(season), .by = gauge) |>
  print()

# ---- Scenarios -----------------------------------------------------------------------------------
scenarios <- use |>
  reframe(scenario = names(scenario_probs),
          prob = unname(scenario_probs),
          discharge = quantile(discharge, scenario_probs, names = FALSE),
          n_seasons = n(),
          .by = c(gauge, month, mpos, month_lab)) |>
  arrange(gauge, mpos, prob)
write_csv(scenarios, here("results", "tables", "12_discharge_scenarios.csv"))

cat("\nDischarge scenarios, monthly max (cfs):\n")
scenarios |>
  mutate(discharge = round(discharge)) |>
  select(gauge, scenario, month_lab, discharge) |>
  pivot_wider(names_from = month_lab, values_from = discharge) |>
  print(n = Inf)

# ---- Figure --------------------------------------------------------------------------------------
prep <- \(d) d |>
  mutate(gauge = factor(gauge, levels = gauges$gauge),
         across(any_of("scenario"), \(x) factor(x, levels = names(scenario_probs))),
         across(any_of("telemetry_season"),
                \(x) if_else(x, "telemetry seasons (2013 to 2024)", "other seasons")))

scen_labs <- set_names(paste0(names(scenario_probs), " (", 100 * scenario_probs, "th pct)"),
                       names(scenario_probs))

p12 <- ggplot() +
  geom_point(data = prep(use), aes(mpos, discharge, shape = telemetry_season),
             position = position_jitter(width = 0.15, height = 0, seed = 1),
             colour = "grey35", alpha = 0.6, size = 1.6) +
  geom_line(data = prep(scenarios), aes(mpos, discharge, colour = scenario, linetype = scenario),
            linewidth = 1.1) +
  geom_hline(data = fit_range |> pivot_longer(c(lo, hi), values_to = "cfs") |> prep(),
             aes(yintercept = cfs), linetype = "dotted", colour = "grey40") +
  facet_wrap(~gauge, scales = "free_y", labeller = as_labeller(gauge_labels)) +
  scale_x_continuous(breaks = seq_along(month_labels), labels = month_labels) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = pal, labels = scen_labs, name = NULL) +
  scale_linetype_manual(values = ltys, labels = scen_labs, name = NULL) +
  scale_shape_manual(values = c("telemetry seasons (2013 to 2024)" = 16, "other seasons" = 1),
                     name = NULL) +
  labs(x = NULL, y = "Monthly max discharge (cfs, log scale)",
       title = "Discharge scenarios for the survival covariate, November to May",
       subtitle = paste0("Seasons ", min(use$season), " to ", max(use$season),
                         "; SRJPEdata ", substr(sha, 1, 7)),
       caption = paste(
         "Dots: each season's monthly max (filled: seasons with telemetry releases in the survival fit).",
         "Lines: 16th, 50th and 84th percentiles across seasons for each gauge and month.",
         "Dotted grey: discharge range of the telemetry releases in the survival fit; outside it survival is extrapolated.",
         sep = "\n")) +
  theme_bw() +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))

ggsave(here("results", "figures", "12_discharge_scenarios.png"), p12, width = 13, height = 5, dpi = 300)
message("Wrote results/figures/12_discharge_scenarios.png")
