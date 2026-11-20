# Changes from Flora's scripts (Phase 2, September 2026)

Every edit in the code is marked with a `FIX` comment. Stan files are unchanged.

## Fixes applied

### A1. Travel times were attached to the wrong fish (02_GetData.R)
Travel times (`TTfR`, `TTind`, `ObsTT`) were built in FishID order, while every
other model input (`CH`, `Sz`, `CovX`, `rgind`, `ReachKM_ind`) is in `d_Sac_sort`
order (year, release location). Stan reads row i of `TTind` as fish i of
`d_Sac_sort`, so 7,857 of 9,161 mainstem fish carried a travel time from a fish
in a different study. Feather/Butte had the same problem (5,358 fish).
- Travel times now built by joining onto `d_Sac_sort` / `d_FeaBut_sort` in model
  row order, with `stopifnot()` guards on row count and order.
- `TravTime_Sac.csv` / `TravTime_FeaBut.csv` now bind columns in the same order
  (the old versions had `RelGp` misaligned).
- **Affects:** all travel-time parameters in CovIndCont, CovIndCont_syre and the
  combined model (including `T_bSz_shared`, which soft-shares the CJS travel-time
  size slope into CWT). Survival estimates are not affected. CovIndWY and NoCov
  have no travel-time component.

### A2. Flow standardized two different ways (02, 05, 08, 09)
CJS `Maxflow.z` was standardized in 01 over every month in the Knights Landing
table (mean 13,782, sd 8,282 cfs). The CWT covariate and the flow prediction
grids used mean/sd across CJS fish (13,445 / 6,203). Same flow, different z
(30,000 cfs: 1.96 vs 2.67).
- 02 defines one reference scale (`flow_mu_sac`, `flow_sd_sac`, plus Butte and
  Feather) and stops if it doesn't reproduce the fitted z-scores.
- 05, 08, 09 use the reference scale. 05's flow-plot x-axis now matches the grid
  the forecasts were computed on.

### A2b. CWT flow came from a different gauge (01_PrepData.R)
`hatchery_release$monthly_max_flow` is pre-joined in SRJPEdata at the Red Bluff
gauge (USGS 11377100). The CJS covariate is Knights Landing (USGS 11390500).
They differ by up to 78,100 cfs over 190 matching year-months, so the combined
model's shared flow slope linked two different covariates.
- CWT groups now join the same `forecast_covariates` Knights Landing series by
  release year and month. Red Bluff kept as `monthly_max_flow_rbdd`.
- 01 prints how many CWT groups have no Knights Landing flow (dropped downstream).
- Gauge choice is a science decision; both tag types must use the same one.

### A3. 08 used a different fork-length scale than the 03 fit
03 standardizes FL on mainstem + tributary fish; 08 used mainstem only. Fixed.

### A4. CWT pCap prior no longer falls back silently (06, 09)
If the Knights Landing pCap `.rds` was missing, 06/09 used a vague normal(0, 1.5)
prior and ran anyway, leaving CWT survival unidentified (Phase 1 Stage 2). Both
now stop with a message.

### A5. CWT groups released when Knights Landing wasn't reading tags (01, new 01b)
The Knights Landing CWT recovery data start in December 1999 and have almost no
reads in 2013 and none in 2022, although the trap was fishing. Groups released
then (128 of ~520: 1994-1999, 2013, 2022) had zero recoveries by construction,
which the binomial likelihood treated as "none survived". This inflated
`sd_cwtSre`, pulled pCap down (posterior ~0.04% vs prior ~0.21%) and, since
2013 and 2022 were dry years, likely inflated the CWT flow slope.
- New `01b_filter_cwt_reading_years.R`: keeps a group only if at least one CWT
  (any tag, any group) was read at Knights Landing within 60 days of its release.
  Writes the filtered `data/hatchery_relrec.csv`, keeps the unfiltered
  `data/hatchery_relrec_all.csv`, and writes a per-group and per-year table to
  `results/tables/cwt_groups_kl_reading_filter*.csv`.
- `01` now writes `hatchery_relrec_all.csv` and sources `01b`. `01b` also runs
  standalone (no ERDDAP pull).
- 2023 is kept (Knights Landing read 167 CWTs that year) but 19 of 26 groups had
  zero recoveries; worth checking whether processing was partial.

### Housekeeping
- 05, 08 read `fit_CovIndCont_MaxFlow_FL.Rdata` (the name 03 saves), not the
  dated `_0226` copy.
- helper_fxt.R: Program MARK setup no longer hard-codes Flora's user folder and
  only runs if MARK is installed. RMark is optional.

## Rerun order
1. `01_PrepData.R`, or just `01b_filter_cwt_reading_years.R` if 01 already ran (check the "CWT release groups with no Knights Landing flow" message)
2. `scripts/checks/check_A1_A2.R`: A1 `identical_order` should be TRUE for both;
   A2 `MaxflowSac.z` mean/sd will still show ~-0.04 / 0.75 (expected: fish are a
   subset of months), but the two scale rows are no longer mixed anywhere.
3. `03_Run_CJS_Model.R` (run 15), then `04`, `05`
4. `06_Run_CWT_Model.R` (CWT flow covariate changed), then `07`, `08`
5. `09_Run_Combined_CJS_CWT_Model.R`

## Not yet done (from the Phase 2 code review)
Code bugs:
- `01`: `release_trib` case_when has a broken string, so Butte_blwsanborn and
  Butte_weir2 fish are labeled Feather (affects `indsitet`, 05 plots; `trib_ind`
  is fine).
- `01`: CWT filter spells `"BATTLE CREEK NFK WILDAT"`; the data say `WILDCAT`.
- `CovIndCont_syre.stan`: tributary fish index site effects 1-2 (Sac, Battle)
  instead of 5-6.

Model changes (from Phase 1 simulation results):
- Combined prediction (generated quantities): replace flat 50/50 blend with
  precision weighting in R (Phase 1 Script 22); drop the unidentifiable Sac-Delta
  reach from `mean(S_bCJS)`; define the endpoint as release to Knights Landing
  (product of reach survivals) instead of the ~400 km fish-count-weighted
  average; propagate release-group RE variance into forecasts.
- Soft sharing: report CWT-specific slopes alongside shared; run a no-sharing
  comparison fit; prior sensitivity on `sigma_*_model`; watch
  `sigma_T_bSz_model` convergence. Small sigma is not evidence sharing is safe.
- pCap: one scalar pCap for 1988-2024 CWT groups; check trial years vs CWT
  years; consider group-level or flow-dependent pCap.
- Travel time: cumulative times from one fish enter as independent observations;
  CWT travel times truncated at 60 days without a truncated likelihood.
- Restrict size prediction grids to observed sizes (03/05 use 10-150 mm).
- Non-center random effects; add priors to CovIndCont.stan (RE_sd is improper);
  use posterior::rhat / ess_bulk / ess_tail with Rhat > 1.01 or ESS < 400.
- Freeport (Delta reach) flow comes from a local CDEC file, not SRJPEdata.
