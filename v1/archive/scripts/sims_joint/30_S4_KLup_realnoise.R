# Run one scenario of the joint CJS + CWT simulation (S4_KLup_realnoise) as an RStudio background job:
#   rstudioapi::jobRunScript(here::here("scripts", "sims_joint", "30_S4_KLup_realnoise.R"), workingDir = here::here())
# Compile the model in the console first, and start the second job a minute after the first.
scen_run <- "S4_KLup_realnoise"
source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))
