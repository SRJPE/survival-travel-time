# Fit only model E2 (tag-effect prior normal(0, 0.15)) on scenario S3, as an RStudio background job.
# S3 data files already exist (run2), so this job only adds the E2 fits.
#   rstudioapi::jobRunScript(here::here("scripts", "sims_joint", "30_S3_E2_only.R"), workingDir = here::here())
scen_run <- "S3_KLup_trialnoise"
models_run <- "E2_tagprior_artifact"
source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))
