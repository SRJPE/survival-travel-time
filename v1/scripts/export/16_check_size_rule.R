# ---- srjpe_survival: export step 16, checks on the small-fish survival rule
# ---- made from scripts/07c_check_lorenzen_size_rule.R by setup/build_from_old_project.R (exact edits listed there)

# ============================================================
# 07c_check_lorenzen_size_rule.R
#
# Checks on the small-fish survival rule (predict_surv_fl). No model fitting.
#   1. At 70 mm and above, all three rules give exactly predict_surv().
#   2. Survival vs fork length, three rules side by side (the team figure):
#      extrapolate (current), clamp at 70 mm, Lorenzen. One site per curve,
#      mean discharge, no release-group effect.
#   3. How the gamma fit (G1) does on the small CWT groups: observed vs
#      expected recoveries by size class. This is the evidence behind c below
#      70 mm; under 60 mm it is a handful of RBDD groups (2003, 2007, 2008).
#
# Inputs : results/survival_export/survival_pars.csv.gz, survival_sites.csv (05c)
#          results/survival_export/lorenzen_pars.csv.gz (13)
#          results/fit_CWT_gamma_G1.Rdata (thumb drive), results/fit_CWT_oneway_Hmed.Rdata (06o)
#          scripts/survival_predict_functions.R
# Outputs: results/tables/07c_lorenzen/07c_identity_check.csv
#          results/tables/07c_lorenzen/07c_survival_by_rule.csv
#          results/tables/07c_lorenzen/07c_cwt_obs_vs_exp_by_size.csv
#          results/tables/07c_lorenzen/07c_cwt_small_groups.csv
#          results/figures/07c_lorenzen/07c_survival_by_rule.png
# A minute or two (uses 2,000 of the 8,000 draws for the figure).
# ============================================================

library(tidyverse)
library(here)
source(here("R", "survival_predict_functions.R"))

ex      <- here("results", "survival_export")
fig_dir <- here("results", "figures", "07c_lorenzen")
tab_dir <- here("results", "tables", "07c_lorenzen")
walk(c(fig_dir, tab_dir), \(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))

L_ref      <- 70
fl_grid    <- seq(30, 130, by = 5)
fl_table   <- c(35, 40, 50, 60, 70, 90)
sites_show <- c("ubc", "mill creek", "okie dam", "lfc")   # Battle, Mill, Butte, Feather
site_lab   <- c(ubc = "Battle Ck (mainstem curve)", `mill creek` = "Mill Ck (mainstem curve)",
                `okie dam` = "Butte Ck", lfc = "Feather River (LFC)")
rule_lab   <- c(extrapolate = "Extrapolate telemetry slope (current)",
                clamp       = "Clamp at 70 mm",
                lorenzen    = "Lorenzen, c from CWT")
rule_col   <- c(extrapolate = "#999999", clamp = "#E69F00", lorenzen = "#0072B2")   # Okabe-Ito
rule_lty   <- c(extrapolate = "dashed", clamp = "dotdash", lorenzen = "solid")
probs      <- c(0.1, 0.5, 0.9)

pars  <- read_csv(file.path(ex, "survival_pars.csv.gz"), show_col_types = FALSE)
sites <- read_csv(file.path(ex, "survival_sites.csv"), show_col_types = FALSE)
lz    <- read_csv(file.path(ex, "lorenzen_pars.csv.gz"), show_col_types = FALSE)
stopifnot(all(sites_show %in% sites$site), nrow(lz) == nrow(pars))

# ---- Grid: mean discharge for each gauge, no release-group effect --------------------
use_draws <- pars$draw[seq(1, nrow(pars), length.out = 2000)]
grid <- expand_grid(draw = use_draws, site = sites_show, fl_mm = fl_grid) |>
  left_join(select(sites, site, cfs = mux), by = "site") |>
  mutate(e1 = 0, e2 = 0, e3 = 0)

res <- names(rule_lab) |>
  set_names() |>
  map(\(r) predict_surv_fl(grid, pars, sites, lz, small_fish = r, L_ref = L_ref) |> mutate(rule = r)) |>
  list_rbind()

# ---- 1. Identity at and above 70 mm ----------------------------------------------------
ident <- res |>
  filter(fl_mm >= L_ref) |>
  select(draw, site, fl_mm, rule, S) |>
  pivot_wider(names_from = rule, values_from = S) |>
  summarise(max_diff_clamp = max(abs(clamp - extrapolate)),
            max_diff_lorenzen = max(abs(lorenzen - extrapolate)), .by = site)
write_csv(ident, file.path(tab_dir, "07c_identity_check.csv"))
print(ident)
stopifnot("rules differ at or above 70 mm" = all(ident$max_diff_clamp < 1e-12, ident$max_diff_lorenzen < 1e-12))
cat("Check 1 passed: all rules equal predict_surv() at 70 mm and above.\n")

# ---- 2. Survival by rule ---------------------------------------------------------------
summ <- res |>
  summarise(q10 = quantile(S, probs[1]), med = quantile(S, probs[2]), q90 = quantile(S, probs[3]),
            .by = c(site, fl_mm, rule))
write_csv(summ |> filter(fl_mm %in% fl_table) |> mutate(across(q10:q90, \(x) round(x, 4))) |>
            arrange(site, fl_mm, rule),
          file.path(tab_dir, "07c_survival_by_rule.csv"))

cat("\nSurvival to the Sacramento station, posterior median, mean discharge, no release-group effect:\n")
summ |>
  filter(fl_mm %in% fl_table) |>
  select(site, fl_mm, rule, med) |>
  mutate(med = round(med, 3)) |>
  pivot_wider(names_from = rule, values_from = med) |>
  arrange(site, fl_mm) |>
  print(n = Inf)

p <- summ |>
  mutate(site = factor(site, levels = sites_show, labels = site_lab[sites_show])) |>
  ggplot(aes(fl_mm, med, colour = rule, fill = rule, linetype = rule)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = L_ref, linetype = "dotted", colour = "grey40") +
  facet_wrap(~ site, scales = "free_y") +
  scale_colour_manual(values = rule_col, labels = rule_lab) +
  scale_fill_manual(values = rule_col, labels = rule_lab) +
  scale_linetype_manual(values = rule_lty, labels = rule_lab) +
  labs(x = "Fork length at the RST (mm)", y = "Survival, RST to Sacramento station",
       colour = NULL, fill = NULL, linetype = NULL,
       title = "Survival of fish smaller than the tagged fish: three rules",
       subtitle = "Mean discharge, typical release group; median and 80% interval",
       caption = paste0("Dotted: 70 mm, about the smallest tagged fish. Above it the three rules are identical.\n",
                        "Lorenzen: S(L) = S(70)^((70/L)^c), c from the CWT gamma fit (G1).")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(file.path(fig_dir, "07c_survival_by_rule.png"), p, width = 10, height = 7, dpi = 300)

# ---- 3. G1 fit on the small CWT groups ---------------------------------------------------
g1_file <- here("results", "fit_CWT_gamma_G1.Rdata")
ow_file <- here("results", "fit_CWT_oneway_Hmed.Rdata")
if (file.exists(g1_file) && file.exists(ow_file)) {
  e <- new.env(); load(g1_file, envir = e); pKL <- rstan::extract(e$fit, "pKL")[[1]]; rm(e)
  e <- new.env(); load(ow_file, envir = e); grps <- e$cwt_groups; rm(e); invisible(gc())
  stopifnot("G1 and one-way group counts differ" = ncol(pKL) == nrow(grps))

  exp_draws <- sweep(pKL, 2, grps$Nrel, `*`)            # expected recoveries, draws x groups
  grps <- grps |>
    mutate(col = row_number(),
           size_class = cut(avg_length, c(0, 60, 70, 80, Inf), right = FALSE,
                            labels = c("< 60 mm", "60 to 70 mm", "70 to 80 mm", "80 mm +")))

  oe <- grps |>
    summarise(groups = n(), years = paste(sort(unique(year)), collapse = " "),
              released = sum(Nrel), observed = sum(cwt_recaptures),
              cols = list(col), .by = size_class) |>
    mutate(exp = map(cols, \(cc) {
      v <- quantile(rowSums(exp_draws[, cc, drop = FALSE]), probs, names = FALSE)
      tibble(exp_q10 = v[1], exp_med = v[2], exp_q90 = v[3])
    })) |>
    unnest(exp) |>
    select(-cols) |>
    mutate(OE = observed / exp_med) |>
    arrange(size_class)
  write_csv(oe, file.path(tab_dir, "07c_cwt_obs_vs_exp_by_size.csv"))
  cat("\nG1 gamma fit: observed vs expected CWT recoveries by size class:\n")
  print(oe |> select(-years) |> mutate(across(where(is.double), \(x) signif(x, 3))))

  small <- grps |>
    filter(avg_length < 60) |>
    mutate(exp_med = map_dbl(col, \(cc) median(exp_draws[, cc]))) |>
    transmute(year, relloc_area, avg_length, Nrel, observed = cwt_recaptures,
              exp_med = signif(exp_med, 3))
  write_csv(small, file.path(tab_dir, "07c_cwt_small_groups.csv"))
  cat("\nCWT groups under 60 mm:\n"); print(as.data.frame(small))
} else {
  message("Skipping check 3: need ", g1_file, " and ", ow_file)
}

cat("\nFigure:", file.path(fig_dir, "07c_survival_by_rule.png"), "\nTables:", tab_dir, "\n")
