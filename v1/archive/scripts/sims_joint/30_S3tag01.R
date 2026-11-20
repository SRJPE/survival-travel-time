# Scenario S3_tag01 (S3 with true tag effect -0.1): models E and E2, as an RStudio background job.
# This job simulates and writes the replicate data files; 30_S3tag01_G_only.R fits G1 on them.
#   rstudioapi::jobRunScript(here::here("scripts", "sims_joint", "30_S3tag01.R"), workingDir = here::here())
scen_run <- "S3_tag01"
models_run <- c("E_joint_trialnoise_artifact", "E2_tagprior_artifact")
source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))
