# 10_sim_recovery_cjs.R
#
# Simulation-recovery test of the CJS survival + travel-time model on the REAL
# study design. Each replicate:
#   1. takes the real fish, release groups, years, fork lengths, flows and reach
#      distances from 02_GetData.R, keeping a fraction of fish in every group
#      (thin_frac) so fits run in a few hours instead of overnight,
#   2. simulates capture histories and travel times from known "true" values
#      (the posterior means of the corrected 03 fit, see sim_functions.R),
#   3. fits the production model (CovIndCont.stan) and the non-centered
#      revision (CovIndCont_nc.stan) to the SAME simulated data,
#   4. saves estimates, truth, coverage and sampler diagnostics to one small
#      .rds per model per replicate in results/sims/cjs_recovery/.
#
# Safe to stop and restart: finished replicates are skipped. Models are
# interleaved (production rep 1, nc rep 1, production rep 2, ...) so both
# build up evenly if the run is cut short.
#
# Run the smoke test first (smoke_test <- TRUE): a tiny design and short chains
# that finish in minutes and prove the pipeline works on this machine. Then set
# smoke_test <- FALSE and start the real run.
#
# Summarise any time (even mid-run) with 11_summarise_sim_recovery.R.
#
# Needs: the project's data/ folder (01 outputs), tidyverse, SRJPEdata, and cmdstanr with
# CmdStan installed (or rstan, if backend <- "rstan").

library(here)

source(here("scripts", "02_GetData.R")) # clears the workspace, builds the real design
source(here("scripts", "sims", "sim_functions.R"))

# ---- the one setting to change (set AFTER 02, which clears the workspace) ----
smoke_test <- TRUE # TRUE = quick pipeline check; FALSE = the real run
backend <- "cmdstanr" # "cmdstanr" (recommended on Windows) or "rstan"

cfg <- if (smoke_test) {
  list(thin_frac = 0.03, min_per_group = 5, iter = 300, chains = 2, n_reps = 1)
} else {
  list(thin_frac = 0.20, min_per_group = 20, iter = 1000, chains = 3, n_reps = 20)
}
cfg$cores <- min(cfg$chains, max(1, parallel::detectCores() - 1))
cfg$base_seed <- 20260925
cfg$truth <- truth_03
cfg$backend <- backend

# CovIndCont_sim.stan = CovIndCont.stan with per-fish quantities kept out of the
# saved output (same parameters, priors and likelihood; see its header)
models <- c(production = "CovIndCont_sim.stan", noncentered = "CovIndCont_nc.stan")
model_paths <- c(production = here("scripts", "sims", "CovIndCont_sim.stan"),
                 noncentered = here("scripts", "sims", "CovIndCont_nc.stan"))
stopifnot(all(file.exists(model_paths)))

out_dir <- here("results", "sims", "cjs_recovery", if (smoke_test) "smoke" else "run1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

design_full <- make_design()
message("Real design: ", design_full$main$N, " mainstem fish, ", design_full$trib$N,
        " Butte/Feather fish, ", design_full$Nrg, " + ", design_full$NrgT,
        " release groups, ", design_full$Nyrs, " years")

compiled <- compile_models(model_paths, cfg$backend)

jobs <- tidyr::expand_grid(rep_id = seq_len(cfg$n_reps), model_name = names(models))

run_start <- Sys.time()
pwalk(jobs, \(rep_id, model_name) {
  run_replicate(rep_id, model_name, compiled[[model_name]], design_full, cfg, out_dir)
  done <- read_sim_results(out_dir)$sampler
  avg <- mean(done$minutes)
  left <- nrow(jobs) - nrow(done)
  message(sprintf("  %d of %d fits done; average %.0f min per fit; about %.1f h left",
                  nrow(done), nrow(jobs), avg, left * avg / 60))
})
message("All done in ", round(as.numeric(difftime(Sys.time(), run_start, units = "hours")), 1), " h")
