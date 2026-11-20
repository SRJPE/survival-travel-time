# Survival and travel time files for runSRJPE

Made by `v1/` in the SRJPE survival-travel-time repo (`scripts/export/19_write_handoff.R`). `MANIFEST.csv`
lists every file, the step that made it, its checksum, the fits it came from and the
SRJPEdata commit.

Survival and travel time are from the RST to the Sacramento station (Delta entry), at
any fork length and monthly max discharge, for the 7 RSTs in `survival_sites.csv`.
Discharge outside the fitted range is clamped to its ends and flagged (`clamped`).

| File | What |
|---|---|
| `survival_pars.csv.gz` | survival parameter draws (8,000), one row per draw |
| `survival_sites.csv` | RSTs: survival curve, distances, gauge, discharge scaling, fitted range |
| `travel_time_pars.csv.gz` | travel-time parameter draws, same `draw` index as survival |
| `travel_time_sites.csv` | `survival_sites.csv` plus reach lengths `km1`-`km3` |
| `lorenzen_pars.csv.gz` | Lorenzen exponent `c` per draw; only for `small_fish = "lorenzen"` |
| `survival_predict_functions.R` | `predict_survival()`, `predict_travel_time()`, `draw_re()` |
| `discharge/` | dry / average / wet seasons and monthly max discharge by season |
| `fit_checksums.csv` | md5 of the fits these files came from (fits are in SharePoint) |

## Use

```r
library(dplyr); library(tidyr); library(readr)
source("survival_predict_functions.R")
pars     <- read_csv("survival_pars.csv.gz")
sites    <- read_csv("survival_sites.csv")
tt_pars  <- read_csv("travel_time_pars.csv.gz")
tt_sites <- read_csv("travel_time_sites.csv")

# release-group effects for survival (e1-e3) and travel time (t1-t3)
re <- draw_re(pars, sites, months = c(11, 12, 1:5), mode = "month", tt_pars = tt_pars)

# wk: one row per site x week with month, fl_mm and cfs
out <- wk |>
  inner_join(re, by = c("site", "month"), relationship = "many-to-many") |>
  predict_survival(pars, sites) |>          # S
  predict_travel_time(tt_pars, tt_sites)    # tt (median days), tt_mean, tt_sdlog
```

## Fish smaller than the tagged fish (70 mm)

Default is `small_fish = "extrapolate"`: the telemetry size slope extended as a straight
line. To use the Lorenzen rule (hazard scaled by (70 / L)^c from the CWT analysis):

```r
lz <- read_csv("lorenzen_pars.csv.gz")
predict_survival(x, pars, sites, lz = lz, small_fish = "lorenzen")
```

`"clamp"` (70 mm survival for all smaller fish) is an upper bound. At 70 mm and above
the three rules are identical.

Travel time also extrapolates by default. `predict_travel_time(x, tt_pars, tt_sites,
small_fish = "clamp")` gives small fish the 70 mm travel time instead, as a bound.
There is no Lorenzen-type rule for travel time.

## Regenerating these files

From the fits (in SharePoint), without refitting: export steps 09 to 19 in `v1/`.
Step 10 rebuilds an intermediate file of about 11 GB, so it needs a machine with plenty
of memory. Refitting the telemetry model takes about
8 hours (refit steps 01 to 03). See `v1/WORKFLOW.md` in survival-travel-time.
