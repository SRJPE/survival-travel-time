# Simulation recovery on the real design (Phase 2)

Fits the production CJS model (scripts/sims/CovIndCont_sim.stan, identical to
scripts/CovIndCont.stan except per-fish quantities are not saved) and the non-centered
revision (scripts/sims/CovIndCont_nc.stan) to data simulated from known values
on the REAL design: every release group, year, fork length, flow and reach
distance from 02_GetData.R, with a fraction of fish kept per group.

Truth = posterior means of the corrected 03 fit (25 Sep 2026); see
sim_functions.R. Both models see the same simulated data in each replicate.

## One-time setup (cmdstanr; rstan failed to compile on the laptop's R 4.6)
```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan(cores = 2)   # downloads and builds CmdStan, 10 to 20 min
```
To use rstan instead, set `backend <- "rstan"` in 10_sim_recovery_cjs.R.

## Run
1. Open the project (.Rproj) so here::here() points at the project folder.
2. Smoke test: open 10_sim_recovery_cjs.R, leave `smoke_test <- TRUE`, source it.
   Takes a few minutes; writes results/sims/cjs_recovery/smoke/.
3. Real run: set `smoke_test <- FALSE`, source again. 20 replicates x 2 models,
   20% of fish, 3 chains x 1000 iterations. Prints time per fit and hours left.
4. Any time: source 11_summarise_sim_recovery.R (set run <- "smoke" or "run1").
   Writes summary_by_parameter.csv, sampler_summary.csv, coverage_plot.png,
   error_plot.png into the run folder.

Stopping is safe: each finished fit is saved as its own .rds and skipped on
restart. Keep the computer from sleeping.

## Tested
Pipeline tested end to end on a synthetic design with both rstan 2.32 and
cmdstanr 0.9 / CmdStan 2.36 (smoke settings): both models compile, sample and save. In that tiny test the
production model had 6 divergences and max Rhat 1.38; the non-centered model
had 0 divergences, max Rhat 1.04, and ran about 3 times faster. Not evidence
yet (one tiny replicate), which is what the real run is for.

## Scope
Tests the Stan models on the real design. It does not test 02's own assembly
of detection histories (the simulator writes model-ready data directly); that
is a later step.
