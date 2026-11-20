# Fit only the G models (tag effect fixed at 0) on scenario S3, as an RStudio background job,
# alongside the job that is still fitting A, C, E, F on S3:
#   rstudioapi::jobRunScript(here::here("scripts", "sims_joint", "30_S3_G_only.R"), workingDir = here::here())
# Uses replicate data files written by the other S3 job; waits for a replicate's data file if it
# does not exist yet. Restart-safe like the other wrappers (existing fits are skipped).
scen_run <- "S3_KLup_trialnoise"
models_run <- c("G_tag0_noise_artifact", "G_tag0_noise_real")
wait_for_data <- TRUE
source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))
