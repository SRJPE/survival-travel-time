# ---- srjpe_survival: export step 11, add travel-time spread to the stm
# ---- made from scripts/05d_add_tt_spread_to_stm.R by setup/build_from_old_project.R (exact edits listed there)

# ============================================================
# 05d_add_tt_spread_to_stm.R
# Adds the travel-time pieces the arrival-timing step needs to the stm
# built by 05, without rerunning 05. Additions only; nothing in the stm
# is changed or removed.
#
#   Pro_sd, Pro_sdT : fish-level lognormal SD of cumulative travel time
#                     (mainstem, tributary), one value per posterior draw.
#                     Individual TT ~ lognormal(log(TTForecast*), Pro_sd).
#   Xvec_raw        : mainstem discharge grid in cfs (Wilkins Slough)
#   XvecT_raw       : Butte and Feather discharge grids in cfs
#   Xsz_raw         : fork-length grid in mm
#
# Draw alignment: rstan::extract() orders draws by the permutation stored
# in the stanfit, so extracting Pro_sd here gives the same draw order 05
# used. The script checks this against stm$T_bReach and TT_bCov and stops
# if they differ.
#
# Inputs : results/fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata (03g)
#          results/stm_MaxFlowFL_KLdelta.Rdata (05)
# Output : results/stm_MaxFlowFL_KLdelta.Rdata, overwritten with the additions
#          (a copy of the original is kept as stm_MaxFlowFL_KLdelta_pre05d.Rdata)
# Then   : 12, 12b, 13.
# Memory : loads the fit, extracts, frees it, then loads the stm (about 10 GB).
#          Restart R first.
# ============================================================

library(rstan)
library(here)
library(purrr)

fit_file <- here("results", "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata")
stm_file <- here("results", "stm_MaxFlowFL_KLdelta.Rdata")
bak_file <- here("results", "stm_MaxFlowFL_KLdelta_pre05d.Rdata")

# ---- 1. Draws from the fit, in 05's order -------------------------------------
e <- new.env(); load(fit_file, envir = e)
post <- rstan::extract(e$fit, pars = c("Pro_sd", "Pro_sdT", "T_bReach", "TT_bCov"))
rm(e); invisible(gc())

# ---- 2. stm --------------------------------------------------------------------
if (!file.exists(bak_file)) stopifnot(file.copy(stm_file, bak_file))
load(stm_file)                                   # object: stm

stopifnot("draw order differs from 05" =
            identical(post$T_bReach, stm$T_bReach) && identical(post$TT_bCov, stm$TT_bCov))

# ---- 3. Raw grids --------------------------------------------------------------
add <- list(
  Pro_sd    = post$Pro_sd,
  Pro_sdT   = post$Pro_sdT,
  Xvec_raw  = stm$Xvec * stm$sdx + stm$mux,
  XvecT_raw = cbind(butte   = stm$XvecT[[1]] * stm$sdxB + stm$muxB,
                    feather = stm$XvecT[[2]] * stm$sdxF + stm$muxF),
  Xsz_raw   = stm$Xsz * stm$sdsz + stm$musz
)

# Grid ends must equal the observed range of tagged-fish discharge (05 builds them that way)
source(here("scripts", "refit", "02_get_data.R"), local = (d <- new.env()))
stopifnot(isTRUE(all.equal(range(add$Xvec_raw), range(d$MaxflowSac, na.rm = TRUE))),
          isTRUE(all.equal(range(add$XvecT_raw[, "butte"]),   range(d$MaxflowB[d$trib_ind == 1]))),
          isTRUE(all.equal(range(add$XvecT_raw[, "feather"]), range(d$MaxflowF[d$trib_ind == 2]))),
          isTRUE(all.equal(range(add$Xsz_raw), c(10, 150))))
rm(d)

stm <- c(stm, add)
save(stm, file = stm_file)

# ---- 4. Report -----------------------------------------------------------------
q <- \(x) round(quantile(x, c(0.025, 0.5, 0.975)), 3)
cat("Added:", paste(names(add), collapse = ", "), "\n")
cat("Pro_sd  (mainstem):  ", q(add$Pro_sd), "\n")
cat("Pro_sdT (tributary): ", q(add$Pro_sdT), "\n")
cat("Fish-level 80% range around a 10-day median (mainstem, posterior median Pro_sd):",
    round(10 * exp(qnorm(c(0.1, 0.9)) * median(add$Pro_sd)), 1), "days\n")
cat("Mainstem grid (cfs):", round(range(add$Xvec_raw)), "\n")
cat("Original stm kept as", bak_file, "\n")
