# Run the joint CJS + CWT simulation (draft, 30 Sep 2026)
#
# For each scenario x replicate: simulate one data set, then fit all four models to it.
#   A_cwt_only_constp : CWT only, one constant KL pCap with prior (like 06)
#   B_joint_constp    : joint CJS + CWT + efficiency trials, pCap not flow-dependent
#   C_joint_flowp     : joint, pCap depends on daily flow (the proposed model)
#   D_cjs_only        : telemetry only (like 03, mainstem)
# Scenarios: S1 KL efficiency constant in flow; S2 KL efficiency falls with flow.
#
# Restart-safe: one .rds per fit; existing fits are skipped.

library(here)
source(here("scripts", "sims_joint", "sim_joint_functions.R"))

# Settings ---------------------------------------------------------------------------
cfg <- list(
  backend = "cmdstanr",          # "rstan" also works (used for testing on Linux)
  n_reps = 20,
  scenarios = "S1_const_pKL",
  models = names(model_flags),
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

  walk(cfg$models, \(mod) {
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
