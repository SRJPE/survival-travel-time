# Survival model, V1 (`survival-travel-time/v1/`)

Juvenile Chinook survival and travel time on the Sacramento River, from acoustic
telemetry (Flora Cordoleani's CJS model, corrected) with a CWT-based size rule for
fish smaller than the tagged fish. This repo makes the files handed to the
integrated model (runSRJPE). Bill Pine, November 2026.

Open `v1/v1.Rproj` before running anything; every script uses `here()`, which resolves to `v1/`
(the `.here` file marks it). Flora's original code and data in the rest of the repo are untouched.

## What is handed off (written by step 19)

`handoff/<date>/runSRJPE/` goes into the runSRJPE repo. A JPE run needs only these files.

| File | Made by | What |
|---|---|---|
| `survival_pars.csv.gz` | 12 | survival parameter draws (8,000), one row per draw |
| `survival_sites.csv` | 12 | the 7 RSTs: survival curve, distances, gauge, discharge scaling, fitted range |
| `travel_time_pars.csv.gz` | 12b | travel-time parameter draws, same draw index as survival |
| `travel_time_sites.csv` | 12b | `survival_sites.csv` plus reach lengths in km |
| `lorenzen_pars.csv.gz` | 13 | Lorenzen exponent c for each draw (fish under 70 mm, optional rule) |
| `survival_predict_functions.R` | R/ | `predict_survival()`, `predict_travel_time()`, `draw_re()` |
| `discharge/` | 15 | dry / average / wet seasons and monthly max discharge by season |
| `fit_checksums.csv` | 09 | md5 of the fits the files came from |
| `README.md`, `MANIFEST.csv` | 19 | how to use them; sizes, checksums, SRJPEdata commit, versions |

`handoff/<date>/sharepoint/` goes to a SharePoint folder for the same date (not in git): the telemetry and CWT
gamma fits, needed only to re-export without refitting, plus an archive copy of
`handoff/runSRJPE/`.

A sent handoff is never changed. Each run of 19 writes a new dated folder and stops
if that date already exists. Code changes after a handoff go through git; a new set of
files for SharePoint gets a new date.

Survival and travel time are from the RST to the Sacramento station (Delta entry), at
any fork length and monthly max discharge. Discharge outside the fitted range is clamped
and flagged. Below 70 mm (the smallest tagged fish), `predict_survival()` uses
`small_fish = "extrapolate"` by default (telemetry slope, straight line); `"lorenzen"`
(c from the CWT fit) and `"clamp"` (70 mm survival) are options. See
`scripts/export/16_check_size_rule.R`. The stm is no longer handed off.

## Two tiers

**Export** (from the saved fits, no refitting). This is what normally runs.

| Step | Script | Was (old project) | Output |
|---|---|---|---|
| 09 | `export/09_check_fits.R` | new | checks the fits against `data/fit_checksums.csv` (writes it the first time) |
| 10 | `export/10_predict_stm.R` | 05 | `results/stm_MaxFlowFL_KLdelta.Rdata` (about 11 GB) and figures |
| 11 | `export/11_add_tt_spread.R` | 05d | adds travel-time spread (Pro_sd) to the stm |
| 12 | `export/12_export_survival_pars.R` | 05c | `results/survival_export/survival_pars.csv.gz`, `survival_sites.csv` |
| 12b | `export/12b_export_travel_time_pars.R` | new | `results/survival_export/travel_time_pars.csv.gz`, `travel_time_sites.csv` |
| 13 | `export/13_export_lorenzen_pars.R` | 07b | `results/survival_export/lorenzen_pars.csv.gz` |
| 15 | `export/15_discharge_scenarios.R` | 12 | `results/tables/12_*.csv` |
| 16 | `export/16_check_size_rule.R` | 07c | checks and figure, small-fish rule |
| 17 | `export/17_check_travel_time.R` | new | travel time by fork length, both small-fish rules: tables and figure |
| 19 | `export/19_write_handoff.R` | new | `handoff/<date>/runSRJPE/`, `handoff/<date>/sharepoint/`; stops if that date exists |

If the stm already exists (copied by setup), skip 10 and 11. Restart R before 10, 11,
12 and 12b (they load the stm). 12 and 12b each end with a check that the exported
draws reproduce the stm's forecasts; they stop if not.

**Refit** (hours; only when the telemetry or CWT data change).

| Step | Script | Was | Time | Output |
|---|---|---|---|---|
| 00 | `refit/00_save_taggedfish.R` | new | seconds | `data/TaggedFish.csv.gz` (once; needs internet) |
| 01 | `refit/01_prep_data.R` (+ 01b) | 01, 01b | minutes | `data/Sac_data.csv`, `FeaBut_data.csv`, detection histories, `hatchery_relrec*.csv` |
| 02 | `refit/02_get_data.R` | 02 | | sourced by 03 to 07 and 10 to 12 |
| 03 | `refit/03_fit_telemetry.R` | 03g | about 8 h | `results/fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata` |
| 04 | `refit/04_check_telemetry.R` | 03h | minutes | `results/tables/diagnostics_03g.csv` |
| 05 | `refit/05_fit_cwt_weekly.R` | 06w | about 1 h | `results/fit_CWT_weekly.Rdata` (+ pooled fit, reused if present) |
| 06 | `refit/06_fit_cwt_oneway.R` | 06o | about 30 min | `results/fit_CWT_oneway_Hmed.Rdata` |
| 07 | `refit/07_fit_cwt_gamma.R` | 07 | about 30 min | `results/fit_CWT_gamma_G1.Rdata` |

01 ends by comparing the rebuilt model data with the 8 Oct 2026 data behind the
current fit (`data/reference/`). `release_trib` is expected to differ (typo fixed);
anything else differing means the data changed and the fit should be rerun.

## Inputs (`data/`, in git)

- `SacRiverAllYears/`, `FeatherRiverAllYears/`, `ButteCreekAllYears/`: ERDDAP detections
  pulled 8 Oct 2026, gzipped. `helper_fxt.R` can still pull live, but 01 reads these.
- `TaggedFish.csv.gz`: ERDDAP tagged-fish table (step 00).
- SRJPEdata, pinned: `remotes::install_github("SRJPE/SRJPEdata@d3584db2a7230a85d3fdb06b5df1b9f794c3a197")`.
  Discharge (`forecast_covariates`), CWT releases and Knights Landing recoveries come from it.
- `kl_weekly_pcap.csv`: Knights Landing weekly capture probability (logit mean, SD) by
  KL year and run week, from the BT-SPAS-X abundance fits (1996-2025). The script that
  made it is part of the BT-SPAS-X work, not this repo. Needed only for refit steps 05 to 07.
- `pCap_one_site_skew_re_knights landing.rds`: older BT-SPAS-X KL efficiency fit, prior
  for the pooled CWT model in step 05.
- `FPTflow.csv.gz`: Freeport discharge, read by 01 but replaced in 02 by the Knights
  Landing series (`delta_flow_source = "knights"`).
- `reference/`: model data behind the current fit, for 01's check.

## Input provenance

| File | Source |
|---|---|
| `data/SacRiverAllYears/`, `FeatherRiverAllYears/`, `ButteCreekAllYears/` | ERDDAP JSATS detections, archived 8 Oct 2026 (gzipped) |
| `data/TaggedFish.csv.gz` | ERDDAP FED_JSATS_taggedfish, saved by `refit/00` (checked against the fitted fish) |
| `data/Sac_data.csv`, `FeaBut_data.csv`, `DetectionHistory*.csv`, `hatchery_relrec*.csv` | made by `refit/01` and `01b` from the above plus SRJPEdata |
| `data/reference/*_20261008.csv` | the model data behind the current fit, for `refit/01`'s comparison |
| `data/FPTflow.csv.gz`, `data/cdec-water-year-type-jun-2025.xlsx` | Freeport flow; CDEC water-year types |
| `data/kl_weekly_pcap.csv` | Knights Landing weekly logit capture probability (posterior mean and SD) by KL year and run week, from the BT-SPAS-X abundance fits (25 KL years, 1996-2025). Copied from the simple joint simulations' input folder; the script that made it belongs to the BT-SPAS-X work, not this repo |
| `data/pCap_one_site_skew_re_knights landing.rds` | Knights Landing pooled capture-probability posterior from the BT-SPAS-X work; prior for `refit/05` |

Known quirk: `refit/01` gives about 23 warnings "multiple rows match for GEN=Releasepoint:
first taken" from base `reshape()` when building detection histories (some fish have more
than one release-point row). The release point only marks release, so nothing is lost, and
the data behind the current fit were built the same way (01's reference check matches).

`archive/` holds the scripts from the original analysis folder (history, not runnable);
see `archive/INDEX.csv`.

## Not in git (`results/`, `handoff/*/sharepoint/`)

Fits, the stm, CmdStan output, tables and figures. The fits go to the SharePoint folder
with each handoff; `data/fit_checksums.csv` (in git) records which fits the current
handoff came from, and step 09 checks against it.

## Built from

The old working folder (`Documents/Ecometric/.../20260924_flora_josh_mortality`) by
`setup/build_from_old_project.R`, which lists every edit made to the old scripts.
History, side analyses (flow vs year, joint models, travel-time checks) and simulations
stay in the old folder.
