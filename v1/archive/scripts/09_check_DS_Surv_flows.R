# 09_check_DS_Surv_flows.R
# Are the fixed flows in For_Input.xlsx (DS_Surv tab: 15,945 mainstem, 896 Butte,
# 21,028 Feather) the stm's mean flows (the centre of its flow scaling)?
# Prints every single-number object in the team's stm and ours, flags any within
# 1% of a DS_Surv value, and lists objects whose names suggest flow scaling.
library(tidyverse)
library(here)

ds_surv <- c(mainstem = 15945, Butte = 896, Feather = 21028)
stm_files <- c(team = here("data", "stm.Rdata"),
               ours = here("results", "stm_MaxFlowFL_KLdelta.Rdata"))

scalars <- function(f) {
  e <- new.env(); load(f, envir = e)
  objs <- mget(ls(e), envir = e)
  # objects saved as a list (e.g. stm <- list(...)) are unpacked one level
  objs <- c(objs[!map_lgl(objs, is.list)], flatten(objs[map_lgl(objs, is.list)]))
  keep(objs, \(x) is.numeric(x) && length(x) == 1) |>
    imap(\(x, nm) tibble(object = nm, value = as.numeric(x))) |> list_rbind()
}
res <- imap(stm_files, \(f, who) scalars(f) |> mutate(stm = who)) |> list_rbind() |>
  mutate(matches = map_chr(value, \(v) {
    hit <- names(ds_surv)[abs(v - ds_surv) / ds_surv < 0.01]
    if (length(hit) == 0) "" else paste(hit, collapse = ",")
  }))
cat("\nSingle-number objects that match a DS_Surv flow (within 1%):\n")
print(as.data.frame(filter(res, matches != "")))
cat("\nSingle-number objects whose names look like flow scaling (mu, sd, mean, flow):\n")
print(as.data.frame(filter(res, str_detect(tolower(object), "^mu|^sd|mean|flow"))))
