# Scenario S3_tag01: model G1 (tag effect fixed at 0) only, alongside 30_S3tag01.R.
# Waits for each replicate's data file written by the other job.
#   rstudioapi::jobRunScript(here::here("scripts", "sims_joint", "30_S3tag01_G_only.R"), workingDir = here::here())
scen_run <- "S3_tag01"
models_run <- "G_tag0_noise_artifact"
wait_for_data <- TRUE
source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))
