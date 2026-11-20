# Run one scenario of the joint CJS + CWT simulation (S3_KLup_trialnoise) as an RStudio background job:
#   rstudioapi::jobRunScript(here::here("scripts", "sims_joint", "30_S3_KLup_trialnoise.R"), workingDir = here::here())
# Compile the model in the console first, and start the second job a minute after the first.
scen_run <- "S3_KLup_trialnoise"
source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))
