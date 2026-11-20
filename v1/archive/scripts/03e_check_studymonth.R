# ============================================================
# 03e_check_studymonth.R
# Checks the study x month fit (03d): convergence, where the chain-1
# divergences are, and the key estimates next to the study-group fit.
# ============================================================
library(rstan); library(purrr); library(dplyr); library(tibble); library(here)

load(here("results", "fit_CovIndCont_studymonth_ncRE_MaxFlow_FL.Rdata")); fit_sm <- fit
load(here("results", "fit_CovIndCont_centered_alignedTT_long_MaxFlow_FL.Rdata")); fit_st <- fit
rm(fit); gc()

pars <- c("S_bReach","T_bReach","S_bCov","TT_bCov","S_bSz","T_bSz","S_bTrib","T_bTrib",
          "S_bCovT","TT_bCovT","RE_sd","RE_sdT","TTRE_sd","TTRE_sdT","Pro_sd","Pro_sdT","muPb","sdPb")

# 1. Convergence and estimates, study x month vs study groups
summ <- function(f, label) {
  summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary |>
    as.data.frame() |> rownames_to_column("parameter") |>
    transmute(parameter, fit = label, median = `50%`, q2.5 = `2.5%`, q97.5 = `97.5%`, rhat = Rhat, ess = n_eff)
}
tab <- bind_rows(summ(fit_sm, "study_month"), summ(fit_st, "study")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  tidyr::pivot_wider(names_from = fit, values_from = c(median, q2.5, q97.5, rhat, ess))
print(as.data.frame(tab))

# 2. Per-chain medians of the random-effect SDs: is chain 1 somewhere different?
sd_pars <- c("RE_sd","RE_sdT","TTRE_sd","TTRE_sdT","sdPb")
dr <- as.array(fit_sm, pars = sd_pars)
apply(dr, c(2, 3), median) |> t() |> round(3) |> print()

# 3. Divergent draws (chain 1): which SDs are small there?
sp <- get_sampler_params(fit_sm, inc_warmup = FALSE)
div <- sp[[1]][, "divergent__"] == 1
cat("chain 1 divergences:", sum(div), "\n")
rbind(divergent = apply(dr[div, 1, , drop = FALSE], 3, median),
      chain1_all = apply(dr[, 1, ], 2, median),
      other_chains = apply(dr[, -1, ], 3, median)) |> round(3) |> t() |> print()
