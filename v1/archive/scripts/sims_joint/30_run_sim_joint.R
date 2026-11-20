# Run the joint CJS + CWT simulation (draft, 30 Sep 2026)
#
# For each scenario x replicate: simulate one data set, then fit all four models to it.
#   A_cwt_only_constp : CWT only, one constant KL pCap with prior (like 06)
#   B_joint_constp    : joint CJS + CWT + efficiency trials, pCap not flow-dependent
#   C_joint_flowp     : joint, pCap depends on daily flow (the proposed model)
#   D_cjs_only        : telemetry only (like 03, mainstem)
#   E_joint_trialnoise_artifact, F_joint_trialnoise_real : joint with skew-normal
#     efficiency-trial noise (BT-SPAS-X KL structure), S3 and S4 only
#   G_tag0_noise_artifact, G_tag0_noise_real : as E and F with the acoustic tag effect
#     fixed at 0 (exchangeability, as in Perry et al. and Hance et al.), S3 and S4 only.
#     Truth in the sims has a tag effect of -0.3, so G also shows the cost of assuming none.
#   E2_tagprior_artifact : as E with a literature-informed prior on the tag effect, normal(0, 0.15)
# Scenario S3_tag01: S3 with a smaller true tag effect (-0.1); models E, E2, G1.
# Scenarios: S1 KL efficiency constant in flow; S2 KL efficiency falls with flow;
#   S3, S4 KL efficiency from the BT-SPAS-X KL fit (rises with flow, skewed trial noise),
#   noise a trial artifact (S3) or real weekly variation (S4).
# Models per scenario: scenario_models in sim_joint_functions.R
#
# Restart-safe: one .rds per fit; existing fits are skipped.

library(here)
source(here("scripts", "sims_joint", "sim_joint_functions.R"))

# Settings ---------------------------------------------------------------------------
cfg <- list(
  backend = "cmdstanr",          # "rstan" also works (used for testing on Linux)
  n_reps = 20,
  # set scen_run before sourcing (the 30_S1 / 30_S2 wrappers do) to run one scenario only,
  # so two R sessions can run the scenarios in parallel without touching the same files
  scenarios = if (exists("scen_run")) scen_run else names(scenarios),
  # set models_run before sourcing to fit only some models (e.g. the 30_S3_G_only wrapper)
  models_run = if (exists("models_run")) models_run else NULL,
  # set wait_for_data <- TRUE before sourcing for a second job on a scenario whose data another
  # job is writing (the *_G_only wrappers)
  wait_for_data = exists("wait_for_data") && isTRUE(wait_for_data),
  chains = 4,
  iter_warmup = 750,
  iter_sampling = 750,
  adapt_delta = 0.9,
  out_dir = here("results", "sims", "joint_cjs_cwt", "run2")   # run1 = draft test fits (rstan, 2 x 400/400)
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

model <- if (cfg$backend == "cmdstanr") {
  cmdstanr::cmdstan_model(here("scripts", "sims_joint", "joint_cjs_cwt.stan"))
} else {
  rstan::stan_model(here("scripts", "sims_joint", "joint_cjs_cwt.stan"))
}

# One scenario x replicate -----------------------------------------------------------
run_one_rep <- function(scen, r) {
  data_file <- file.path(cfg$out_dir, sprintf("data_%s_rep%02d.rds", scen, r))
  # A second job on the same scenario (e.g. G only) waits for the first job to write the
  # data file instead of simulating it too (same seed, but two jobs writing one file can collide)
  if (cfg$wait_for_data) {
    while (!file.exists(data_file) ||
           difftime(Sys.time(), file.mtime(data_file), units = "secs") < 60) Sys.sleep(60)
  }
  rep <- if (file.exists(data_file)) {
    readRDS(data_file)
  } else {
    set.seed(1000 * match(scen, names(scenarios)) + r)
    tr <- modifyList(default_truth(), scenarios[[scen]])
    x <- simulate_replicate(default_design(), tr)
    saveRDS(x, data_file)
    x
  }
  truth <- true_values(rep$truth)

  mods <- scenario_models[[scen]]
  if (!is.null(cfg$models_run)) mods <- intersect(mods, cfg$models_run)
  walk(mods, \(mod) {
    fit_file <- file.path(cfg$out_dir, sprintf("fit_%s_rep%02d_%s.rds", scen, r, mod))
    if (file.exists(fit_file)) return(invisible())
    f <- fit_model(make_stan_data(rep, model_flags[[mod]]), model,
                   backend = cfg$backend, chains = cfg$chains,
                   iter_warmup = cfg$iter_warmup, iter_sampling = cfg$iter_sampling,
                   adapt_delta = cfg$adapt_delta, seed = r)
    out <- f$summary |>
      left_join(truth, by = "parameter") |>
      mutate(scenario = scen, rep = r, model = mod,
             n_div = f$n_div, n_treedepth = f$n_treedepth, minutes = f$minutes)
    saveRDS(out, fit_file)
    message(sprintf("%s rep %02d %-18s %5.1f min, %d divergences, max Rhat %.3f",
                    scen, r, mod, f$minutes, f$n_div, max(out$rhat, na.rm = TRUE)))
  })
}

# Run: replicates in the outer loop so every scenario and model fills in evenly -----
expand_grid(r = seq_len(cfg$n_reps), scen = cfg$scenarios) |>
  pwalk(\(r, scen) run_one_rep(scen, r))
