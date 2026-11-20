# ============================================================
# 06L_lorenzen_check.R
#
# Lorenzen check linking telemetry and CWT survival. No model fitting:
# reads the saved fits.
#
# Idea: mortality per day scales with fish length,
#   M(L) = M_r * (L_r / L)^c        (Lorenzen; c = 1 is the classic value)
# so survival over a migration of T days is exp(-M(L) * T).
#   1. Telemetry (03f fit) gives each mainstem fish's survival from release
#      to the Sacramento station (SurvRelSac) and travel time (TT_RelSac).
#      With c fixed, M_r = sum(-log S_i) / sum(T_i * (L_r / L_i)^c) per
#      posterior draw (total hazard over size-adjusted days at risk).
#   2. Each CWT group gets Lorenzen survival to KL from its own length and
#      its travel time to KL (06w weekly fit; travel time was the same in
#      the pooled and weekly fits).
#   3. Compare (a) with the CWT fits' survival (pooled and weekly), and
#      (b) observed CWT recoveries with those expected from Lorenzen
#      survival x BT-SPAS-X weekly capture probability (prior only, not
#      updated by the CWT data). O/E near 1 = telemetry + Lorenzen and
#      BT-SPAS-X agree with the recoveries.
#
# c = 0.5, 1, 1.5 as sensitivity. L_r = 70 mm (typical CWT fish).
# Assumes mortality accrues per day (not per km), no growth during
# migration, no flow effect on the CWT side, and CWT fish exchangeable
# with telemetry fish apart from length (no tag effect).
#
# Inputs (results/): fit_CovIndCont_studymonth_ncRE_MaxFlow_FL.Rdata (03f),
# fit_CWT_weekly.Rdata and fit_CWT_pooled_samegroups.Rdata (06w).
# Fork lengths of telemetry fish come from 02.
# Outputs: results/tables/lorenzen_*.csv, results/figures/lorenzen_*.png
# ============================================================

rm(list = ls())

library(tidyverse)
library(here)
library(rstan)

source("scripts/02_GetData.R")   # for FL (mainstem telemetry fish, model row order)

c_values <- c(0.5, 1, 1.5)
L_r <- 70
set.seed(1)

# Telemetry: per-fish survival and travel time to Sacramento, per draw -----------------
e <- new.env(); load(here("results", "fit_CovIndCont_studymonth_ncRE_MaxFlow_FL.Rdata"), envir = e)
S_tel <- rstan::extract(e$fit, "SurvRelSac")[[1]]    # draws x fish
T_tel <- rstan::extract(e$fit, "TT_RelSac")[[1]]
rm(e); gc()
stopifnot(ncol(S_tel) == length(FL))
ok <- is.finite(FL)
S_tel <- S_tel[, ok]; T_tel <- T_tel[, ok]; L_tel <- FL[ok]
hazard_tel <- rowSums(-log(S_tel))                     # total hazard per draw
nd <- nrow(S_tel)
message("Telemetry fish: ", length(L_tel), "; length ", paste(range(L_tel), collapse = " to "), " mm")

# CWT fits ----------------------------------------------------------------------------
load(here("results", "fit_CWT_weekly.Rdata"))          # fit_cwt_weekly, cwt_data, cwt_groups, segs, wk_tbl, ...
e <- new.env(); load(here("results", "fit_CWT_pooled_samegroups.Rdata"), envir = e)
fit_cwt_pooled <- e$fit_cwt_pooled; rm(e)

pick <- \(x) x[round(seq(1, nrow(x), length.out = nd)), , drop = FALSE]   # match telemetry draw count
T_cwt <- pick(rstan::extract(fit_cwt_weekly, "pTT")[[1]])                # median travel time to KL, draws x groups
S_weekly <- pick(rstan::extract(fit_cwt_weekly, "surv_cwt")[[1]])
S_pooled <- pick(rstan::extract(fit_cwt_pooled, "surv_cwt")[[1]])
L_cwt <- cwt_groups$avg_length

# year averaging (groups x years)
Y <- model.matrix(~ factor(year) - 1, cwt_groups)
colnames(Y) <- sort(unique(cwt_groups$year))
Yavg <- sweep(Y, 2, colSums(Y), "/")

# BT-SPAS-X capture probability for each group, prior only --------------------------
# weekly mean of inv_logit(normal(lp_mean, lp_sd)) and, for reference, the median
wk_tbl <- wk_tbl |>
  mutate(p_mean = map2_dbl(lp_mean, lp_sd, \(m, s) mean(plogis(rnorm(20000, m, s)))),
         p_med = plogis(lp_mean))
mu_g <- log(apply(T_cwt, 2, median))
sd_tt <- median(rstan::extract(fit_cwt_weekly, "Pro_sd")[[1]])
p_bt <- segs |>
  left_join(wk_tbl |> select(wk, p_mean, p_med), by = "wk") |>
  mutate(f = plnorm(hi, mu_g[grp], sd_tt) - plnorm(lo, mu_g[grp], sd_tt)) |>
  group_by(grp) |>
  summarise(p_mean = sum(f * p_mean), p_med = sum(f * p_med)) |>
  arrange(grp)
stopifnot(identical(p_bt$grp, cwt_groups$grp))

# Lorenzen for each c -----------------------------------------------------------------
q3 <- \(x) c(median = median(x), q2.5 = unname(quantile(x, 0.025)), q97.5 = unname(quantile(x, 0.975)))

lorenzen <- function(c) {
  M_r <- hazard_tel / drop(T_tel %*% ((L_r / L_tel)^c))      # per day at L_r, per draw
  S_L <- exp(-sweep(T_cwt, 2, (L_r / L_cwt)^c, `*`) * M_r)   # draws x CWT groups (row d uses M_r[d])
  E_mean <- sweep(S_L, 2, cwt_groups$Nrel * p_bt$p_mean, `*`) %*% Y   # expected recoveries by year
  E_med  <- sweep(S_L, 2, cwt_groups$Nrel * p_bt$p_med, `*`) %*% Y
  list(M_r = M_r, S_L = S_L, E_mean = E_mean, E_med = E_med)
}
L <- set_names(map(c_values, lorenzen), paste0("c", c_values))

# 1. Mortality scale ---------------------------------------------------------------------
Mr_tbl <- imap(L, \(x, nm) tibble(c = nm, !!!as.list(q3(x$M_r)),
                                   surv_per_day_at_median = exp(-median(x$M_r)))) |>
  list_rbind() |>
  mutate(across(where(is.numeric), \(x) signif(x, 3)))
print(as.data.frame(Mr_tbl))
write.csv(Mr_tbl, here("results", "tables", "lorenzen_Mr.csv"), row.names = FALSE)

# 2. Mean group survival to KL: Lorenzen vs the CWT fits ------------------------------
surv_draws <- c(map(L, "S_L"), list(cwt_weekly = S_weekly, cwt_pooled = S_pooled))
overall_tbl <- imap(surv_draws, \(S, nm) tibble(source = nm, !!!as.list(q3(rowMeans(S))))) |>
  list_rbind() |>
  mutate(across(where(is.numeric), \(x) signif(x, 3)))
print(as.data.frame(overall_tbl))
write.csv(overall_tbl, here("results", "tables", "lorenzen_overall_survival.csv"), row.names = FALSE)

surv_year <- imap(surv_draws, \(S, nm) {
  m <- S %*% Yavg
  tibble(source = nm, year = as.integer(colnames(m)),
         median = apply(m, 2, median), q2.5 = apply(m, 2, quantile, 0.025),
         q97.5 = apply(m, 2, quantile, 0.975))
}) |> list_rbind()
write.csv(surv_year, here("results", "tables", "lorenzen_survival_by_year.csv"), row.names = FALSE)

# 3. Observed vs expected recoveries by year -------------------------------------------
obs_year <- cwt_groups |> group_by(year) |>
  summarise(groups = n(), released = sum(Nrel), observed = sum(cwt_recaptures))
oe <- imap(L, \(x, nm) tibble(
  c = nm, year = as.integer(colnames(x$E_mean)),
  exp_median = apply(x$E_mean, 2, median),
  exp_q2.5 = apply(x$E_mean, 2, quantile, 0.025),
  exp_q97.5 = apply(x$E_mean, 2, quantile, 0.975),
  exp_btmedian_pcap = apply(x$E_med, 2, median))) |>
  list_rbind() |>
  left_join(obs_year, by = "year") |>
  mutate(OE = observed / exp_median, OE_lo = observed / exp_q97.5, OE_hi = observed / exp_q2.5)
oe_tot <- oe |> group_by(c) |>
  summarise(observed = sum(observed), expected = sum(exp_median),
            expected_btmedian_pcap = sum(exp_btmedian_pcap), OE = observed / expected)
print(as.data.frame(oe_tot), digits = 3)
oe |> filter(c == "c1") |>
  select(year, groups, observed, exp_median, exp_q2.5, exp_q97.5, OE) |>
  mutate(across(where(is.double), \(x) signif(x, 3))) |>
  as.data.frame() |> print()
write.csv(oe, here("results", "tables", "lorenzen_obs_vs_expected_by_year.csv"), row.names = FALSE)

# Figures (Okabe-Ito colours plus shapes) -------------------------------------------------
src_lab <- c(c0.5 = "Lorenzen c = 0.5", c1 = "Lorenzen c = 1", c1.5 = "Lorenzen c = 1.5",
             cwt_weekly = "CWT fit, weekly pCap", cwt_pooled = "CWT fit, pooled pCap")
src_col <- c(c0.5 = "#56B4E9", c1 = "#0072B2", c1.5 = "#009E73",
             cwt_weekly = "#E69F00", cwt_pooled = "#D55E00")
src_shp <- c(c0.5 = 1, c1 = 16, c1.5 = 2, cwt_weekly = 15, cwt_pooled = 17)

p1 <- ggplot(surv_year, aes(year, median, ymin = q2.5, ymax = q97.5,
                            colour = source, shape = source)) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) +
  scale_colour_manual(values = src_col, labels = src_lab) +
  scale_shape_manual(values = src_shp, labels = src_lab) +
  labs(x = NULL, y = "Mean group survival to Knights Landing", colour = NULL, shape = NULL,
       title = "CWT survival: Lorenzen from telemetry vs the CWT fits",
       subtitle = paste0("Lorenzen: telemetry mortality per day scaled to CWT length (L_r = ", L_r,
                         " mm) over CWT travel time")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "lorenzen_survival_by_year.png"), p1,
       width = 10, height = 5.5, dpi = 200, bg = "white")

p2 <- ggplot(oe, aes(year, OE, ymin = OE_lo, ymax = OE_hi, colour = c, shape = c)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.3) +
  scale_y_log10() +
  scale_colour_manual(values = src_col[1:3], labels = src_lab[1:3]) +
  scale_shape_manual(values = src_shp[1:3], labels = src_lab[1:3]) +
  labs(x = NULL, y = "Observed / expected CWT recoveries (log scale)", colour = NULL, shape = NULL,
       title = "CWT recoveries vs Lorenzen survival x BT-SPAS-X capture",
       subtitle = "1 = agreement; BT-SPAS-X weekly mean capture probability, not updated by CWT data") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "lorenzen_obs_vs_expected.png"), p2,
       width = 10, height = 5.5, dpi = 200, bg = "white")
p1; p2
