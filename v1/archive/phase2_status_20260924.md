# srJPE Phase 2 status, updated 30 September 2026

## Where things live
- Project home: `C:\Users\flori\Documents\Ecometric\20260820_flora_mortality\20260924_flora_josh_mortality` (scripts/, data/, results/, CHANGES.md, README.md, .here). Its project file is still named `20260820_flora_mortality.Rproj` (same name as the old root's); renaming it to `20260924_flora_josh_mortality.Rproj` was suggested. Always check `here::here()` ends in `20260924_flora_josh_mortality`.
- Laptop working copy (while on the road): `C:\Users\owner\Documents\20260924_flora_josh_mortality`. Run models from local disk, not from the Google Drive "Other computers" copy of the Dell (`G:/Other computers/dell_work_station/...`). Copy results back to the Dell afterwards.
- The parent folder `20260820_flora_mortality` is the old root: loose root scripts and `results_mystery_file` (38 MB gzip R file, not yet identified) should go to `old/`.
- Latest script bundle: `20260925_sims_and_script_updates.zip` (unzip into the project home): scripts/sims/ plus current 01, 01b, 06, CWT_nc.stan.

## Purpose of the simulation work
All simulation work exists to improve Flora's original production code (02, 05, 06, 07, 09; CovIndCont.stan, CovIndCont_syre.stan, CovIndWY.stan, CWT.stan, Combined_CJS_CWT_FL_SharedMainstem.stan). Order: finish the simulation stages first, then build new versions of her models from what they show. 05's `stm` object is the current interface to the integrated model team; what we hand them (including whether flow stays in) is decided after the sims.

## Fixes applied to Flora's scripts (details in CHANGES.md)
- A1: travel times were attached to the wrong fish (TTfR built in FishID order; model rows in year/release order). 7,857 of 9,161 mainstem fish had another study's travel time; Feather/Butte also misaligned. Fixed in 02 with guards. Affects all travel-time estimates; survival unaffected.
- A2: flow standardized two ways (monthly-table mean 13,782 / sd 8,282 cfs vs per-fish 13,445 / 6,203). One reference scale defined in 02; used in 05, 08, 09.
- A2b: CWT flow came from Red Bluff (USGS 11377100, pre-joined in SRJPEdata hatchery_release) while CJS used Knights Landing (USGS 11390500); up to 78,100 cfs apart. 01 now joins CWT groups to the Knights Landing series.
- A3: 08 used a different fork-length standardization than the 03 fit.
- A4: 06/09 stop if the Knights Landing pCap prior .rds is missing.
- A5: CWT groups released when Knights Landing wasn't reading tags dropped by `01b_filter_cwt_reading_years.R` (any CWT read at KL within 60 days of release). Real run: kept 394 of 521, dropped 127 (1994, 1996-1999, 2013, 2022). 58 zero-recovery groups remain in reading years (2023: 19 of 26; ask about partial processing).
- A6: `CWT_nc.stan` (non-centered group random effects) and 06 set to 4 chains x 4000 iterations (2000 warmup). Run 30 Sep; see results below.
- A7 (30 Sep): 06 fits with cmdstanr instead of rstan. On the laptop (R 4.6), rstan's C++ compile fails and its error handler crashes with "Error in sink(type = "output") : invalid connection", hiding the compiler message. cmdstanr works (stage 1 sims used it). The cmdstanr fit is converted to a stanfit so the rest of 06 is unchanged. rstan::read_stan_csv() cannot parse newer CmdStan CSV headers (`save_warmup = false`), so 06 copies the CSVs with header true/false rewritten as 1/0 (`*_rstan.csv`) and reads those. 06 now also saves `fit_CWT.Rdata` (the original save line was lost in the edit; 07 needs it) and `fit_CWT_cmdstan.rds`; draws are kept in `results/cmdstan_cwt/`. 03 and 09 (which fit models) will need the same swap on the laptop; 05, 07 and 08 only read fits and run under rstan as-is.
- A8 (30 Sep): 05 per-fish predictions omitted fish size. 03 fits `S_bSz*Sz[i]` (survival) and `T_bSz*Sz[i]` (travel time), but 05's pred_surv, pred_surv_per100, SurvRelSac, TT_reach, TT_RelSac, pred_survT, pred_survT_per100 and TT_reachT left both out, so all were computed at mean fork length. Fixed by adding the size terms (mainstem Sz, trib SzT), with `Sz <- (FL - musz)/sdsz` and `SzT <- (FL_T - musz)/sdsz` added after `Xsz` (same fork-length branch and scale as the 03 fit; 03 does not export Sz, it builds it before fitting). Forecast arrays already included size and are unchanged.
- A9 (30 Sep): 05 memory. Replaced `as_draws_array(fit)`, `as_draws_df(fit)` and a full `rstan::extract(fit)` (each several GB because the 03 fit stores per-fish transformed parameters; draws and draws_df were never used) with `rstan::extract(fit, pars = pars_needed)` for the 19 parameters 05 uses. `fit` stays loaded for extract_log_lik().

## Results so far
### 03 CJS-only (CovIndCont run 15, flow + FL, corrected travel-time alignment), finished 25 Sep
- Sampler: 1 divergence in 1,500 draws, 0 treedepth hits.
- Covariate effects converge well (Rhat <= 1.002, ESS > 900): S_bCov 0.70 (0.54 to 0.85) per SD of KL flow on the reference scale; S_bSz 0.10 (0.06 to 0.15); TT_bCov -0.15 (-0.20 to -0.11); T_bSz -0.025 (-0.041 to -0.010); trib flow effects S_bCovT 2.1 (Butte), 2.3 (Feather); TT_bCovT 0.80 (Butte), -0.06 (Feather).
- Problems concentrated in the terminal Sac-Delta reach, as Phase 1 Stage 1 predicted: TTRE_sd[4] Rhat 1.12, ESS 21 (SD near 0, classic funnel); RE_sd[4] Rhat 1.03, ESS 82; S_bReach[4] ESS 93; P_b[,4] poorly identified. S_bReach[3] Rhat 1.014. Many S_RE below ESS 400 with only 1,500 draws.
- Read: main effects usable; reach-4 parameters not interpretable. Production fix: non-center S_RE / TT_RE / P_b, more iterations, keep reach 4 out of predictions.
- Run time under 7 hours (started about 20:10 on 24 Sep, done before 03:12 on 25 Sep).

### 06 CWT-only, Knights Landing flow
| Parameter | Unfiltered | Filtered (A5), centered | Filtered (A5), non-centered (30 Sep) |
|---|---|---|---|
| S_bCov | 0.30 (0.18 to 0.42) | 0.09 (0.00 to 0.19) | 0.09 (0.002 to 0.19) |
| S_bSz | 0.36 (0.23 to 0.50) | 0.06 (-0.04 to 0.15) | 0.06 (-0.04 to 0.16) |
| S_bReach | | | -0.09 (-0.39 to 0.16) |
| sd_cwtSre | 1.17 | 0.82 | 0.82 (0.67 to 0.97) |
| logit_pCap_Sim | -7.78 (0.042%) | -7.51 (0.055%) | -7.49 (-7.79 to -7.07), about 0.056% |
| T_bReach | | | 1.09 (1.06 to 1.13) |
| T_bSz | -0.23 | -0.22 | -0.22 (-0.26 to -0.17) |
| TT_bCov | -0.02 | -0.02 | -0.02 (-0.06 to 0.01) |
| sd_cwtTTre, Pro_sd | | | 0.26, 0.42 |
- Prior on logit pCap: mean -6.14 (about 0.21%), sd 1.39.
- Covariates standardized on the CWT groups' own mean and SD, so slopes are not on the CJS scale.

#### Non-centered run (A6/A7), 30 Sep, laptop
- 394 release groups (2000 to 2024; 372 Battle, 22 RBDD), so the A5 filter was applied.
- 4 chains x 2000 post-warmup draws = 8,000 draws; about 18 minutes (mean chain 938 s).
- Sampler: 0 divergences, 0 treedepth hits. Step size 0.017 to 0.026, treedepth 7 to 8, partly from adapt_delta 0.99; 0.9 would likely halve run time.
- Convergence fixed: worst Rhat 1.004 (T_bReach), lowest ESS 967 (S_bCov), vs Rhat about 1.07 and ESS about 55 in the centered run. All 394 group survivals and travel times Rhat <= 1.001, ESS > 2,400.
- Estimates essentially unchanged from the centered run: the centered fit found the right answer but could not certify it. Same pattern as stage 1 for the CJS model.
- Group-level predictions: median survival release to Knights Landing 0.12 (range 0.008 to 0.62); median travel time 8.6 days (4.6 to 22.7).

#### Flags from the non-centered run
1. The weak CWT flow effect is now a reliable finding: S_bCov 0.09 with the interval just above zero, against CJS 0.70. Scales differ, but a gap this size is not a scaling artifact. The discharge scenarios and the no-sharing fit in 09 need to explain it.
2. pCap separation rests on an assumption. The binomial sees only survival x pCap (same last-interval problem as CJS reach 4). Posterior pCap is much tighter than the prior (SD 0.19 vs 1.39) and sits about one prior SD below its mean. Likely source of the separation: the two release distances (Battle 300 km, RBDD 240 km). With per-km survival assumed equal from both sites and a shared pCap, the recovery ratio between sites identifies per-km survival. That rests on 22 RBDD groups and the equal-per-km-survival assumption. Stage 2 should test this directly (pCap recovery with and without the prior, with and without the RBDD contrast).
3. 2023 mean group survival 0.025, the lowest of any year. 2023 has 19 of 26 zero-recovery groups, suspected partial CWT processing, not mortality. Exclude 2023 or confirm with the KL crew before these survivals are used.
4. Remaining pCap gap (about 4x below the prior mean) is consistent with only a fraction of captured ad-clipped fish having CWTs processed; ask the KL crew.

### CJS vs CWT tension to watch in 09
- Flow on survival: CJS 0.70 vs CWT 0.09 (different standardizations: CJS on the monthly-table reference scale, CWT-only 06 on CWT groups' own mean/sd; the combined model puts both on one scale). Confirmed by the converged non-centered 06. A flow-dependent RST efficiency would bias the CWT flow slope down; this is the first discharge simulation scenario.
- Travel time on size: CJS -0.025 per SD vs CWT -0.22 per SD (again different size scales).
- Phase 1: soft sharing pulls CWT about 30 to 50% toward CJS and sigma may not flag realistic divergence, so the no-sharing comparison fit is essential.
- Quantified in 08 (30 Sep, below): flow disagrees strongly, size on survival agrees, size on travel time disagrees.

### 05 CJS predictions (30 Sep, laptop, with A8 and A9)
- Ran from the 25 Sep 03 fit; wrote `stm_MaxFlowFL.Rdata` and figures.
- Release-to-Sacramento survival by study shows a strong wet/dry split: wet-year groups mostly 0.25 to 0.46, dry-year groups mostly below 0.1, several near zero (MillCk 2015, CNFH_FMR 2021, Rbdd 2021). This is S_bCov 0.70 visible in the data. For telemetry fish, year conditions move survival by about an order of magnitude; if flow is dropped (integrated team), that variation goes into the release-group REs and widens every forecast interval.
- Check: Wild_stock_Chinook_Rbdd_2024 is a wet year at about 0.02 survival (real, or a data issue with that study).
- Trib Sac-Delta survival (reach 4) intervals about 0.05 to 0.9 in most years: prior plus hierarchy, as stage 1 predicted. The few tight years (SB 2017, FR 2019, all three 2023) are years where survival x Delta detection is high, which forces both up; not true separation. Do not use the Sac-Delta panels downstream.
- LOO: elpd_loo -17,821 (SE 122), p_loo 125, looic 35,641. Pareto k count not yet recorded.
- Plotting notes for new versions: wet/dry year lists hard-coded and inconsistent (2014 dry in trib list only); red/blue should become Okabe-Ito.

### 07 CWT predictions (30 Sep, laptop)
- Ran as-is from the non-centered fit_CWT.Rdata.
- **Release site and fish size are almost completely confounded:** RBDD groups 43 to 60 mm, Battle groups 61 to 87 mm, essentially no overlap. Each site's curve is drawn across the full size range but has data on only half. The RBDD-above-Battle survival gap is built in by the distance exponent (same per-km survival), not estimated. A lower recovery rate for small RBDD fish can be explained by the size slope, the distance exponent or pCap, so S_bSz, S_bReach and pCap are harder to interpret than their tight intervals suggest (see 06 flag 2).
- No wet/dry pattern in CWT survival by year: dry 2004 (0.48) and 2010 (0.38) among the highest; wet 2000 (0.08) and 2023 (0.025) among the lowest. Most CWT years (2000 to 2012) predate telemetry; in overlapping years, wet 2017, 2019, 2024 are among the higher recent values.
- Travel time: fits Battle roughly; misses RBDD. Observed RBDD times spread 0 to 55 days with a long tail (fry likely rearing on the way, not directed migration); model peaks about 13 days. The 60-day truncation cuts the tail. RBDD predicted slower than Battle despite being 60 km closer, purely from T_bSz (-0.22) acting on small fish (same confounding). The predicted-vs-observed density figure compares group means with individual fish (Pro_sd 0.42 not included), so the model curve should look narrower; the RBDD misfit is real regardless.
- Reading notes: "by year" bands average interval endpoints across groups (descriptive, not a credible interval); size and flow curves set the group RE to zero (typical group, not a new-group forecast; sd_cwtSre 0.82).
- Implication: CWT informs survival and travel time mainly for 60 to 87 mm Battle fish; below 60 mm it is RBDD fry with different behavior. The Phase 1 fry vs smolt gap is also a site gap.

### 08 CWT vs CJS comparison (30 Sep, laptop)
- Ran as-is. Scales handled correctly: each model's slope converted from cfs and mm on the scale it was fit on (A2, A3 for CJS; CWT's own for CWT).
- Endpoints differ (CWT release to Knights Landing, 240 or 300 km averaged across sites; CJS release to Sacramento, reaches 1 to 3), so compare ratios, not levels. REs set to zero, so bands are typical-group, not new-group, and overlap is not a formal test. Each curve extends beyond its data (rug marks).

**Flow (FL fixed at 86 mm): strong disagreement**
| Flow (cfs) | CJS survival | CWT survival |
|---|---|---|
| 4,480 | 0.041 | 0.136 |
| about 7,600 (lowest CWT group) | about 0.066 | about 0.143 |
| about 24,000 (highest CJS group) | about 0.38 | about 0.186 |
| 29,200 | 0.51 | 0.20 |
| Ratio over full range | 12.6x | 1.5x |
- Not an extrapolation artifact: within the shared range (about 7,600 to 24,000 cfs) CJS rises about 5.8x, CWT about 1.3x. Curves cross near 14,500 cfs.
- Travel time: CJS 6.4 to 4.1 days across the range; CWT 5.2 to 4.9.
- Caution on CJS: the flow slope rests on about 25 release groups over 12 years, few distinct flow values, and flow is tied up with year. The 05 wet/dry figure supports it, but it rests on fewer independent points than 9,000 fish suggest.

**Size (flow fixed at 11,500 cfs)**
- Survival vs size agrees: 43 to 138 mm, CWT 0.089 to 0.28, CJS 0.046 to 0.26. Both roughly triple; curves nearly parallel.
- Travel time vs size disagrees: CJS 6.5 to 4.8 days; CWT 23 days at 43 mm to 0.85 days at 138 mm (over 300 km/day, impossible). Exponential size term extrapolated past the largest CWT group (87 mm). At 100 mm, CWT 3.1 days vs CJS 5.4. Phase 1 travel-time extrapolation warning confirmed on real data; driven by the RBDD fry vs Battle confounding.

**Implications for 09 (combined model shares three slopes softly)**
| Shared slope | Verdict |
|---|---|
| Size on survival | Defensible; the data agree |
| Flow on survival | Problematic; a shared value fits neither data set (pulls CWT up, CJS down) and sigma may not flag it (Phase 1) |
| Size on travel time | Should not be shared; CWT slope confounded with site and extrapolates to impossible values |
- The no-sharing comparison fit is essential; drop size-on-travel-time sharing.
- Candidate explanations for the flow gap (discharge scenarios): (1) Knights Landing trap efficiency drops at high flow; one pCap for all groups would flatten the CWT flow slope (first scenario, leading suspect). (2) Flow measure matches migration differently (monthly max at release for CWT vs max during migration for CJS). (3) Different fish (Battle hatchery fall-run about 70 mm vs tagged mix about 90 mm). (4) CJS slope partly a year effect with only about 12 years.
- Possible later addition: survival per 100 km (no distance exponent; CJS averaged over reaches 1 to 3) for a like-for-like level comparison.

## Simulation recovery, stage 1: CJS on the real design (run 25 to 29 Sep)
### Setup
- Files: scripts/sims/sim_functions.R, 10_sim_recovery_cjs.R, 11_summarise_sim_recovery.R, CovIndCont_nc.stan, README_sims.md.
- Design: real fish, release groups, years, FL, flows, distances from 02; 20% of fish per release group (min 20). Truth = 03 posterior means (sdPb set to 1). Both models (production CovIndCont.stan and non-centered CovIndCont_nc.stan with explicit priors) fit the same simulated data per replicate; 3 chains x 1000 iterations (1,500 post-warmup draws). One .rds per fit in results/sims/cjs_recovery/run1/.
- Completed: 20 non-centered replicates, 5 production replicates. Outputs: sampler_summary.csv, summary_by_parameter.csv, coverage and error plots.

### Speed and sampler behavior
| | Non-centered | Production (Flora's centered) |
|---|---|---|
| Fits | 20 | 5 |
| Minutes per fit | 27 | 609 (about 10 hours) |
| Total compute | 8.9 h | 50.7 h |
| Fits with divergences | 1 of 20 (1 divergence total) | 5 of 5 (mean 52) |
| Fits with treedepth hits | 0 | 2 |
- Same model, data and settings; only the parameterization differs. Production is about 23x slower and never clean. This reproduces Flora's "days of computer time" on a 20% design.
- Sporadic "lognormal_lpdf location parameter is inf" rejection messages during warmup are harmless (exp-then-log of travel time overflowing on extreme warmup proposals); non-centered fits still had 0 divergences.

### Recovery, non-centered (20 replicates)
- Main effects unbiased: S_bCov 0.674 (truth 0.70), TT_bCov -0.152 (-0.15), S_bSz 0.112 (0.10), T_bSz -0.022 (-0.025). Trib intercepts and slopes, travel-time intercepts, Pro_sd and Pro_sdT essentially unbiased.
- Coverage near nominal throughout; values of 0.85 to 0.90 are within binomial noise for 20 replicates (0.85 = 17 of 20). Lowest: RE_sd[3] 0.80.
- Production, even where it runs: S_bCov biased low 22% (0.55 vs 0.70) vs 4% non-centered. Production coverage comes in steps of 0.2 with 5 replicates; use its bias, Rhat and ESS, not coverage.

### Key finding: Sac-Delta (reach 4) is structurally weakly identified
- Delta is the final detection station, so Sac-Delta survival (S_bReach[4]) and Delta detection (P_b[,4]) form the CJS last interval: the likelihood sees only their product. They are separated only by the year-level hierarchy on P_b and the priors. More fish will not fix this; detections below the Delta array would.
- Non-centered: reach-4 parameters converge but with very wide intervals (S_bReach[4] logit CI width 2.6, bias +0.23; muPb[4] CI width 4.5).
- Production: the same parameters fail in every replicate (S_bReach[4] median Rhat 1.08, ESS 30; RE_sd[4] Rhat 1.12, ESS 20; muPb[4] bias +1.36; sdPb[3], sdPb[4] +42% and +48%, Rhat about 1.06). The trib Sac-Delta term fails the same way (RE_sdT[2] Rhat 1.11, ESS 30).
- So the real-data reach-4 problems in 03 have two causes: the centered funnel (fixable by non-centering) and the last-interval confound (structural).
- Flora's mainstem forecasts in 05 use reaches 1 to 3 only (correct). But 05's trib predictions use S_bReach[,Nreaches], and the combined model's predictions use mean(S_bCJS) over all four reaches.

### Other findings
- TTRE_sd[4] is prior-dominated: truth 0.09, both models estimate about 0.38 (+316%). An SD near zero with few groups cannot be pinned down.
- Random-effect SDs biased high 10 to 22% in the non-centered model (RE_sd[3] +22%, RE_sdT[2] +23%): partly the right skew of posterior means for hierarchical SDs, partly fewer fish per group on the 20% design. Relevant to forecasts, which draw a new release group from RE_sd, so intervals may be slightly too wide. Check with posterior medians.
- S_bReach[1] mixes slowly (40% of fits Rhat > 1.01, median ESS 270 in 1,500 draws; a few outlying replicates). RE_sd[1] is large (2.0), so the reach-1 intercept is poorly determined. Needs more iterations, not a model change.

### Caveats
- Data were simulated from the fitted model, so this tests estimability, not robustness to misspecification.
- 20% design: fish-driven intervals should shrink about 2x on the full data; group-level SDs and the reach-4 confound will not improve (same number of release groups).
- Truth came from 03, whose reach-4 values were themselves unreliable. Fine for a recovery test; those truths are not estimates of the real river.
- Tests the Stan models on the real design, not 02's own assembly of detection histories.

### Implications for the new model versions (after all sim stages)
1. Non-center all random effects (S_RE, TT_RE, S_REt, TT_RET, P_b, CWT group effects). Settled; confirmed on real CWT data 30 Sep.
2. Treat Sac-Delta survival as non-identifiable: report survival to the Sacramento station, or add a downstream array. Fix 05 trib predictions and the combined model's mean(S_bCJS).
3. Drop or constrain the reach-4 travel-time random effect.
4. Check RE_sd bias with posterior medians before building forecast intervals from it.
5. Run longer chains (S_bReach[1] ESS).
6. Compute travel time on the log scale (eta + log(km/100), log_sum_exp for cumulative reaches) instead of exp then log.
7. Fit with cmdstanr (works on both machines; rstan fails to compile on the laptop).

## Simulation stages still to run
- Stage 2: CWT recovery (CWT_nc.stan) with the pCap prior. Add: pCap recovery with and without the prior, and with and without the Battle vs RBDD distance contrast (see 06 flag 2), and with site and size separated vs confounded (see 07).
- Stage 3: combined model; soft sharing vs no sharing; 50/50 vs precision-weighted blend; single covariate scales. Add: simulate CJS and CWT from one shared process with flow-dependent KL trap efficiency and a CWT tag offset, and compare the joint model (below) against Flora's soft-sharing model and no sharing; tests the delta_CWT vs p_KL trade-off.
- Discharge scenarios: flow-dependent RST trap efficiency (first), flow vs year/release-group RE confounding, wrong gauge or flow summary, different flow ranges for CJS vs CWT fish, forecast extrapolation.
- Later (real data, after sims): no-flow version of 03 vs the flow model marginalized over telemetry-year flows, since the integrated model team is considering dropping flow.

## Design candidate (30 Sep): joint CJS + CWT model (leading option for the combined-model redesign)
Proposed by bp: instead of estimating CWT survival separately and linking it by soft-shared slopes, use the CWT and telemetry data jointly to estimate one movement and survival process. Template for part of it: Perry et al. 2018, CJFAS 75:1886-1901 (joint Bayesian model of travel time, routing and survival; covariates assigned at reach entry, arrival times imputed for undetected fish; noncentered release-group REs; detection probability a function of discharge).

**Structure**
- CWT fish move through the same reaches under the same survival and travel-time process as telemetry fish, observed through a different observation model. Several telemetry groups share release sites with CWT (ColemanFall, RBDD); KL lies inside CJS reach 3 (Flora's Rmult_Kni, ReachKM_Kni).
- CWT survival to KL: S_g = product over reaches 1, 2 and the partial reach 3 of S_j^(km_j/100), with S_j = logit^-1(S_bReach[j] + flow + size + delta_CWT + group RE). Uses the telemetry model's reach parameters plus a tag-type offset delta_CWT. No separate CWT intercept, no soft sharing, no 50/50 blend.
- CWT observation model (daily recoveries): y_{g,d} ~ Poisson(N_g x S_g x A_g(d) x p_KL(d)).
  - A_g(d): probability of arriving at KL on day d, from the shared travel-time model (Flora's lognormal on summed reach means, discretized to days).
  - p_KL(d): daily KL trap capture probability; can depend on flow, is 0 on days the trap was not fishing, and multiplies by the fraction of ad-clipped fish processed for CWTs that year.
  - Replaces both 06 pieces (binomial on total recoveries, separate travel-time lognormal).

**What it fixes**
1. Flow-dependent trap efficiency (leading suspect for the flat CWT flow slope): modeled in p_KL(d) instead of absorbed into survival.
2. Biased CWT recovery times: fish are recovered more often on high-efficiency days, so observed CWT times are not a random sample of travel times. Also handles the 60-day truncation and trap outages.
3. Last-interval confound: with survival to KL predicted from telemetry-estimated reaches, CWT recoveries mainly inform trap efficiency and between-year variation, not the absolute survival level.

**What CWT would add**
- 394 groups over 22 years (2000 to 2024) vs about 25 telemetry groups over 12 years with flow tied up with year; CWT covers flows and sizes the telemetry does not. CWT informs relative effects (flow, size, year variation); telemetry anchors the absolute level. That split suits a JPE.
- Not resolved by itself: delta_CWT and p_KL trade off (as survival and pCap do now). Either assume exchangeability (delta_CWT = 0; CWT then mostly estimates p_KL) or supply independent trap efficiency and estimate delta_CWT. Independent source: the KL efficiency trials behind the BT-SPAS-X pCap prior, ideally entered as a function of flow, not one constant.

**Assumptions and risks**
1. Exchangeability: acoustic fish about 90 mm with tag burden; CWT fish 43 to 87 mm, lightly tagged hatchery fall-run. delta_CWT and the size term cover part of this, but only across the 61 to 87 mm Battle overlap.
2. RBDD fry do not behave like migrants (07 travel-time misfit). First version: Battle groups only, or a rearing delay / migrant-rearer mixture for small fish.
3. CWT processing fraction by year (2023 zero-recovery problem) needs KL crew records; otherwise it goes into a year effect on p_KL, competing with year variation in survival.
4. Compute: Perry imputed arrival times per fish in JAGS (5 days). Stan cannot impute discrete latent states. CWT groups (tens of thousands of fish each) must be group-level daily probabilities anyway (cheap). Fully time-varying covariates for telemetry fish are a later step; first version keeps Flora's group-level flow.
5. Detection varying with flow: Perry models detection as a function of discharge; our P_b varies only by year. Flow-dependent detection at telemetry stations could affect the CJS flow slope too.
- Context: Perry found survival increased with inflow only in reaches transitioning from tidal to riverine flow, while travel time decreased with flow everywhere. Our reaches are upstream of theirs, so not directly transferable. Recollection to verify: Michel and colleagues reported flow-survival thresholds in the upper Sacramento; could matter for the functional form with about 25 telemetry groups.

**Order**
1. Cheap real-data check first: does KL trap efficiency fall with flow in the efficiency trials? If yes, the CJS vs CWT flow disagreement likely has a mundane explanation and the joint model is clearly the way to go.
2. Stage 3 sims (above) to test identifiability before building on real data.
3. Build as a new model version after the sims (not a change to Flora's code).

## Next (real data)
1. ~~Run 06 with CWT_nc.stan~~ Done 30 Sep (converged; see above).
2. ~~05 (CJS predictions) from the 03 fit~~ Done 30 Sep (A8, A9).
3. ~~07, 08~~ Done 30 Sep. Next: 09 (combined). Consider non-centered combined model and a no-sharing fit before/alongside 09; drop size-on-travel-time sharing. Use cmdstanr on the laptop (09 is centered, adapt_delta 0.995, max_treedepth 15: expect it to be slow).

## Decision (24 Sep): keep working with simulated data
bp: if real-data runs fail, keep working with simulated data, through the production code path. Compute: use thinned designs, cmdstanr if needed, both machines.

## Still open from the code review
- Combined prediction in Combined_CJS_CWT_FL_SharedMainstem.stan: flat 50/50 blend (precision weighting in R, Script 22), mean(S_bCJS) includes unidentifiable reach 4 (confirmed by stage 1), endpoint about 400 km instead of Knights Landing, random-effect variance not propagated.
- Soft sharing: report CWT-specific slopes, no-sharing comparison fit, sigma prior sensitivity, watch sigma_T_bSz_model convergence.
- 05 stm output: per-fish predictions omit the size term; reach-1 exponent inconsistent (Rmult[1,1] vs RmultSac); scaling constants must match 02's reference scale; trib forecasts rely on reach 4.
- 01 release_trib string bug; "WILDAT" vs "WILDCAT"; CovIndCont_syre tributary site-effect indexing; cumulative travel-time pseudo-replication; 60-day CWT travel-time truncation; size prediction grid 10-150 mm; priors for CovIndCont.stan (RE_sd improper); Freeport flow from local CDEC file.
- CWT 2023 zero-recovery groups (partial processing?) before using CWT survivals. Needed for the joint model: yearly CWT processing fraction and KL trap operating days from the KL crew.
- KL trap efficiency vs flow from the efficiency trials (BT-SPAS-X): first check for the joint-model design.
- Stage 2 design: include CWT site-size confounding (simulate with size and site separated vs confounded) alongside the pCap tests.
- Record Pareto k counts from the 05 LOO.
- SurvRelSacSz in 05 uses Rmult[1,1] for reach 1 where everything else uses RmultSac (saved in stm, not plotted).

## Documents
- Plain-language Phase 1 report (docx + pdf, 24 Sep): Figure 5 (stage21_traveltime_curves_plot.png) still missing.
- Figure and table catalog (docx, 24 Sep).
