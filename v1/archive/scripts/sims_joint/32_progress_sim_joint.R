# Progress of the joint simulation: fits done per scenario and model, and minutes per fit
library(here)
library(dplyr)
library(purrr)
run_dir <- here("results", "sims", "joint_cjs_cwt", "run2")
f <- list.files(run_dir, pattern = "^fit_.*\\.rds$", full.names = TRUE)
done <- map(f, \(x) readRDS(x) |> distinct(scenario, rep, model, minutes)) |> list_rbind()
source(here("scripts", "sims_joint", "sim_joint_functions.R"))
total <- 20 * sum(lengths(scenario_models))   # n_reps x models in every scenario
cat(sprintf("%d of %d fits done (%.0f%%)\n", nrow(done), total, 100 * nrow(done) / total))
done |>
  group_by(scenario, model) |>
  summarise(fits = n(), last_rep = max(rep), mean_minutes = round(mean(minutes), 1), .groups = "drop") |>
  print(n = Inf)
