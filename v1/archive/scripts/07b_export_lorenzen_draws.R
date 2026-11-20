# ============================================================
# 07b_export_lorenzen_draws.R
#
# Exports the Lorenzen size exponent c from the CWT gamma fit (G1, run by
# 07_Run_CWT_gamma.R on the laptop) so predict_surv_fl() can set survival
# for fish smaller than the tagged fish (70 mm):
#   S(L) = S(70)^((70 / L)^c)
#
# Why only c: in CWT_oneway_gamma.stan the cumulative hazard is
#   M_r * (70 / L)^c * exp(delta + h) * T^gamma.
# Scaling telemetry survival at 70 mm to length L needs only the ratio
# (70 / L)^c; M_r cancels. gamma (about 0.02) makes the time term about 1.
# delta and the group effects h describe CWT fish specifically (hatchery,
# tag, capture), so they are exported for reference but not applied.
#
# Draw pairing: the gamma fit is one-way (telemetry fixed at its median
# hazard), so its draws are not paired with the telemetry draws. Each of
# the 8,000 survival draws gets one G1 draw at random (fixed seed).
#
# Inputs : results/fit_CWT_gamma_G1.Rdata         (copy from the thumb drive:
#                                                   cwt_gamma/results/)
#          results/survival_export/survival_pars.csv.gz   (05c)
# Output : results/survival_export/lorenzen_pars.csv.gz
#          one row per survival draw: draw, g1_draw, c, gamma, delta, M_r, L_ref
# Seconds. No model fitting.
# ============================================================

library(tidyverse)
library(here)

g1_file   <- here("results", "fit_CWT_gamma_G1.Rdata")
pars_file <- here("results", "survival_export", "survival_pars.csv.gz")
out_file  <- here("results", "survival_export", "lorenzen_pars.csv.gz")
L_ref <- 70
seed  <- 20261118

stopifnot("Copy fit_CWT_gamma_G1.Rdata from the thumb drive (cwt_gamma/results/) into results/" =
            file.exists(g1_file),
          "Run 05c first" = file.exists(pars_file))

# ---- G1 draws ---------------------------------------------------------------------
e <- new.env(); load(g1_file, envir = e)                       # object: fit
keep <- c("c", "gamma_raw", "delta", "M_r")
g1 <- rstan::extract(e$fit, pars = keep) |>
  map(as.vector) |>
  as_tibble() |>
  rename(gamma = gamma_raw) |>
  mutate(g1_draw = row_number(), .before = 1)

sp <- rstan::get_sampler_params(e$fit, inc_warmup = FALSE)
cat("G1 sampler: draws", nrow(g1), "| divergences", sum(map_dbl(sp, \(x) sum(x[, "divergent__"]))), "\n")
rstan::summary(e$fit, pars = c("c", "gamma_raw", "delta"),
               probs = c(0.025, 0.5, 0.975))$summary[, c("2.5%", "50%", "97.5%", "n_eff", "Rhat")] |>
  round(3) |> print()
rm(e); invisible(gc())

# Guard: these should match the 9 Oct result (c 0.46, gamma 0.017, delta 0.13)
stopifnot("c median far from the 9 Oct G1 fit (0.46): is this the right file?" =
            abs(median(g1$c) - 0.46) < 0.1)

# ---- Pair with the survival draws ---------------------------------------------------
n_surv <- nrow(read_csv(pars_file, col_select = "draw", show_col_types = FALSE))
set.seed(seed)
lz <- tibble(draw = seq_len(n_surv),
             g1_draw = sample(nrow(g1), n_surv, replace = n_surv > nrow(g1))) |>
  left_join(g1, by = "g1_draw") |>
  mutate(L_ref = L_ref)

write_csv(lz, out_file)

# ---- What the rule implies ---------------------------------------------------------
# multiplier on the route hazard relative to a 70 mm fish
mult <- expand_grid(fl_mm = c(35, 40, 50, 60, 70), draw = lz$draw) |>
  left_join(select(lz, draw, c), by = "draw") |>
  mutate(m = (L_ref / fl_mm)^c) |>
  summarise(median = median(m), q10 = quantile(m, 0.1), q90 = quantile(m, 0.9), .by = fl_mm)
cat("\nHazard multiplier vs a 70 mm fish, (70 / L)^c (median, 80% interval):\n")
print(mult |> mutate(across(-fl_mm, \(x) round(x, 2))))
cat("Example: if S(70) = 0.30, a 40 mm fish gets 0.30^", round(mult$median[mult$fl_mm == 40], 2),
    " = ", round(0.30^mult$median[mult$fl_mm == 40], 3), "\n", sep = "")
cat("\nWrote", nrow(lz), "rows to", out_file, "\n")
