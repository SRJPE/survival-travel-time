# ============================================================
# scripts/export/17_check_travel_time.R      (new, 19 Nov 2026)
#
# Travel time from each RST to the Sacramento station by fork length, from the
# exported draws (12b) and predict_travel_time(). For Alex (WY2027 arrival
# timing) and for checking how far the size effect is extrapolated below 70 mm.
#
# Fish under 70 mm (the smallest tagged fish) are shown two ways:
#   extrapolate (default)  telemetry size slope extended in a straight line (log scale)
#   clamp                  70 mm travel time (a bound)
# Discharge: the site's mean and mean +/- 1 SD of the tagged-fish discharge
# (monthly max), clamped to the fitted range as predict_travel_time() does.
# Release-group effect = 0 (a typical release group). Posterior: 2,000 of the
# 8,000 draws (set n_draws below).
#
# Inputs : results/survival_export/travel_time_pars.csv.gz, travel_time_sites.csv (12b)
#          R/survival_predict_functions.R
# Writes : results/tables/17_travel_time/
#            17_tt_by_fl.csv          median and 80% interval, every site, length, flow, rule
#            17_tt_summary.csv        mean flow, selected lengths, both rules, ratio to 70 mm
#            17_tt_fish_spread.csv    10th-90th percentile of individual fish, mean flow
#          results/figures/17_travel_time/17_tt_by_fl.png
# ============================================================

library(tidyverse)
library(here)
source(here("R", "survival_predict_functions.R"))

n_draws <- 2000
L_ref   <- 70
fl_grid <- seq(30, 150, by = 5)
fl_show <- c(35, 40, 50, 60, 70, 90, 110)

ex       <- here("results", "survival_export")
tab_dir  <- here("results", "tables", "17_travel_time")
fig_dir  <- here("results", "figures", "17_travel_time")
walk(c(tab_dir, fig_dir), \(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))

tt_pars  <- read_csv(file.path(ex, "travel_time_pars.csv.gz"), show_col_types = FALSE)
tt_sites <- read_csv(file.path(ex, "travel_time_sites.csv"), show_col_types = FALSE)

set.seed(20261119)
draws <- sort(sample(tt_pars$draw, min(n_draws, nrow(tt_pars))))

site_order <- c("ucc", "ubc", "mill creek", "deer creek", "okie dam", "hallwood", "lfc")
site_lab   <- c(ucc = "Clear Creek (UCC)", ubc = "Battle Creek (UBC)", `mill creek` = "Mill Creek",
                `deer creek` = "Deer Creek", `okie dam` = "Butte Creek (Okie Dam)",
                hallwood = "Feather (Hallwood)", lfc = "Feather (LFC)")

# Discharge levels on the site's gauge (mean, +/- 1 SD), back to input cfs
flows <- tt_sites |>
  select(site, mux, sdx, cfs_mult) |>
  cross_join(tibble(flow = factor(c("low", "mean", "high"), levels = c("low", "mean", "high")),
                    k = c(-1, 0, 1))) |>
  transmute(site, flow, cfs = (mux + k * sdx) / cfs_mult)

grid <- flows |>
  cross_join(tibble(fl_mm = fl_grid)) |>
  cross_join(tibble(draw = draws)) |>
  mutate(t1 = 0, t2 = 0, t3 = 0)

pred <- c("extrapolate", "clamp") |>
  set_names() |>
  map(\(r) predict_travel_time(grid, tt_pars, tt_sites, small_fish = r, L_ref = L_ref)) |>
  list_rbind(names_to = "rule")

# Check: the two rules agree at and above 70 mm
chk <- pred |>
  filter(fl_mm >= L_ref) |>
  select(rule, site, flow, fl_mm, draw, tt) |>
  pivot_wider(names_from = rule, values_from = tt) |>
  summarise(max_rel_diff = max(abs(clamp / extrapolate - 1)))
stopifnot("clamp and extrapolate differ at 70 mm and above" = chk$max_rel_diff < 1e-12)
cat("Check passed: rules agree at", L_ref, "mm and above.\n")

by_fl <- pred |>
  summarise(tt_med = median(tt), tt_q10 = quantile(tt, 0.1), tt_q90 = quantile(tt, 0.9),
            sdlog = median(tt_sdlog), cfs_used = first(cfs_used), clamped = first(clamped),
            .by = c(rule, site, flow, fl_mm))
write_csv(by_fl, file.path(tab_dir, "17_tt_by_fl.csv"))

flow_used <- by_fl |>
  distinct(site, flow, cfs_used, clamped) |>
  pivot_wider(names_from = flow, values_from = c(cfs_used, clamped))
cat("\nDischarge used (site gauge, cfs; TRUE = clamped to the fitted range):\n")
print(as.data.frame(flow_used), row.names = FALSE)

summ <- by_fl |>
  filter(flow == "mean", fl_mm %in% fl_show) |>
  left_join(filter(by_fl, flow == "mean", fl_mm == L_ref) |> select(rule, site, tt70 = tt_med),
            by = c("rule", "site")) |>
  mutate(ratio_to_70 = tt_med / tt70) |>
  select(site, fl_mm, rule, tt_med, tt_q10, tt_q90, ratio_to_70) |>
  arrange(factor(site, site_order), fl_mm, rule)
write_csv(summ, file.path(tab_dir, "17_tt_summary.csv"))

cat("\nTravel time to Sacramento (days), mean discharge, median [80% interval]:\n")
summ |>
  mutate(txt = sprintf("%.1f [%.1f-%.1f]", tt_med, tt_q10, tt_q90)) |>
  select(site, fl_mm, rule, txt) |>
  pivot_wider(names_from = rule, values_from = txt) |>
  left_join(filter(summ, rule == "extrapolate") |> select(site, fl_mm, x_vs_70 = ratio_to_70),
            by = c("site", "fl_mm")) |>
  mutate(x_vs_70 = round(x_vs_70, 2)) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Fish-to-fish spread at the posterior median: 10th-90th percentile of individual fish
spread <- by_fl |>
  filter(flow == "mean", fl_mm %in% fl_show, rule == "extrapolate") |>
  mutate(fish_p10 = tt_med * exp(qnorm(0.1) * sdlog),
         fish_p90 = tt_med * exp(qnorm(0.9) * sdlog)) |>
  select(site, fl_mm, tt_med, sdlog, fish_p10, fish_p90)
write_csv(spread, file.path(tab_dir, "17_tt_fish_spread.csv"))
cat("\nIndividual fish, mean discharge (days; 10th-90th percentile around the median):\n")
print(mutate(spread, across(where(is.numeric), \(v) round(v, 1))) |> as.data.frame(), row.names = FALSE)

# ---- Figure ------------------------------------------------------------------------
okabe <- c(low = "#E69F00", mean = "#000000", high = "#0072B2")
fd <- by_fl |> mutate(site = factor(site, site_order, site_lab[site_order]))
ext <- filter(fd, rule == "extrapolate")
clp <- filter(fd, rule == "clamp", fl_mm <= L_ref)

p <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = L_ref, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_ribbon(data = filter(ext, flow == "mean"),
              aes(fl_mm, ymin = tt_q10, ymax = tt_q90), fill = "grey60", alpha = 0.4) +
  geom_line(data = ext, aes(fl_mm, tt_med, colour = flow), linewidth = 0.8) +
  geom_line(data = clp, aes(fl_mm, tt_med, colour = flow), linewidth = 0.6, linetype = "22") +
  geom_vline(xintercept = L_ref, linewidth = 0.3, colour = "grey40") +
  facet_wrap(~site, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = okabe, name = "Discharge",
                      labels = c(low = "mean - 1 SD", mean = "mean", high = "mean + 1 SD")) +
  expand_limits(y = 0) +
  labs(x = "Fork length (mm)", y = "Travel time to Sacramento (days)",
       title = "Travel time from each RST to the Sacramento station",
       subtitle = paste0("Solid: size effect extrapolated below ", L_ref,
                         " mm (shaded; no tagged fish). Dashed: ", L_ref, " mm travel time (clamp).\n",
                         "Band: 80% posterior interval at mean discharge. Typical release group (no release-group effect).")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "17_tt_by_fl.png"), p, width = 11, height = 6.5, dpi = 200)
cat("\nFigure:", file.path(fig_dir, "17_tt_by_fl.png"), "\nTables:", tab_dir, "\n")
