# check_A1_A2.R
# Two quick checks on Flora's pipeline before any refitting.
# Run from the project root (same place you run 03_Run_CJS_Model.R).
# 02_GetData.R clears the workspace, so it is sourced first.

source("scripts/02_GetData.R")
library(tidyverse)

# =============================================================================
# A1. Are travel-time rows aligned with the fish rows the model uses?
# =============================================================================
# Stan treats row i of TTind as fish i of d_Sac_sort (same order as CH, Sz,
# CovX, rgind, ReachKM_ind). TTind is built from release_times, sorted by FishID.

check_alignment <- function(model_ids, tt_ids, label) {
  n_match <- sum(model_ids == tt_ids)
  tibble(
    data            = label,
    n_model_rows    = length(model_ids),
    n_tt_rows       = length(tt_ids),
    same_fish_set   = setequal(model_ids, tt_ids),
    identical_order = identical(model_ids, tt_ids),
    rows_misaligned = length(model_ids) - n_match,
    pct_misaligned  = round(100 * (1 - n_match / length(model_ids)), 1)
  )
}

A1 <- bind_rows(
  check_alignment(d_Sac_sort$FishID,    release_times$FishID,  "Sacramento mainstem"),
  check_alignment(d_FeaBut_sort$FishID, release_timesT$FishID, "Feather/Butte")
)
print(A1)

# First few misaligned rows: which fish's travel time lands on which fish
A1_examples <- tibble(
  row          = seq_len(Nind),
  model_fish   = d_Sac_sort$FishID,
  tt_from_fish = release_times$FishID,
  model_study  = d_Sac_sort$StudyID,
  tt_study     = d_Sac_sort$StudyID[match(release_times$FishID, d_Sac_sort$FishID)]
) %>%
  filter(model_fish != tt_from_fish)

print(head(A1_examples, 10))
cat("Misaligned rows where the travel time comes from a DIFFERENT study:",
    sum(A1_examples$model_study != A1_examples$tt_study, na.rm = TRUE), "\n\n")

# =============================================================================
# A2. Is the CJS flow covariate on the same scale as the CWT flow covariate?
# =============================================================================
# CJS CovX = Maxflow.z, standardized in 01_PrepData.R over every month in
# forecast_covariates (Knights Landing monthly max flow).
# CWT CovX_cwt (09) and the flow prediction grids (05, 08, 09) use the mean/sd
# of MaxflowSac across CJS fish. If the two agree, Maxflow.z has mean ~0, sd ~1
# across fish.

# Recover the mean and sd actually used to build Maxflow.z:
# z = (x - m) / s  ->  x = m + s * z, so the lm intercept = m and slope = s
fit_scale <- lm(MaxflowSac ~ MaxflowSac.z)

A2 <- tibble(
  scale              = c("Used for CJS CovX (01, monthly table)",
                         "Used for CWT CovX_cwt and grids (09/05/08, per fish)"),
  mean_cfs           = c(coef(fit_scale)[[1]], mean(MaxflowSac, na.rm = TRUE)),
  sd_cfs             = c(coef(fit_scale)[[2]], sd(MaxflowSac, na.rm = TRUE))
)
print(A2)

cat("\nMaxflowSac.z across CJS fish: mean =", round(mean(MaxflowSac.z, na.rm = TRUE), 3),
    " sd =", round(sd(MaxflowSac.z, na.rm = TRUE), 3), "(both ~0 and ~1 if scales agree)\n")

# What the mismatch does: the same flow gets two different z values
example_flows <- c(10000, 20000, 30000)
A2_example <- tibble(
  flow_cfs     = example_flows,
  z_as_CJS_sees_it = (example_flows - A2$mean_cfs[1]) / A2$sd_cfs[1],
  z_as_CWT_sees_it = (example_flows - A2$mean_cfs[2]) / A2$sd_cfs[2]
) %>%
  mutate(z_difference = z_as_CWT_sees_it - z_as_CJS_sees_it)
print(A2_example)
