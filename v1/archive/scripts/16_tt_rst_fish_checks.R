# ============================================================
# 16_tt_rst_fish_checks.R
# Two data-only checks (no model, no stm) on how to treat RST fish in the
# arrival step (15: option 1 = reach-1 rate and group SD, 2 = reach-2 rate and
# group SD, 3 = reach-1 rate with reach-2 group SD).
#
#  1. Do fish starting at an RST travel like migrants?
#     Observed travel rate (days per 100 km) in reach 1 (release to Woodson)
#     vs reach 2 (Woodson to Butte) and reach 3 (Butte to Sacramento), by
#     release group, for wild fish caught and released at an RST (Mill, Deer,
#     RBDD wild) vs hatchery releases.
#       Wild RST fish: reach 1 about as fast as reach 2 -> supports option 2 or 3.
#       Wild RST fish slow in reach 1 like hatchery fish -> the slow upper river
#       is the channel, not release behavior -> option 1 or 3.
#     Also the between-group SD of log rate per reach (raw analog of TTRE_sd;
#     it still contains size and discharge differences, so it runs high).
#
#  2. Is the fish-level spread (Pro_sd) inflated by reach 1?
#     Pooled within-group SD of log travel time for the same fish over:
#     release to Woodson, release to Sacramento, Woodson to Sacramento,
#     Butte to Sacramento. Pro_sd is one SD on log cumulative travel time from
#     release, shared by every station. If Woodson to Sacramento is clearly
#     tighter than release to Sacramento, the spread for RST fish should be
#     smaller than Pro_sd.
#
# Inputs : scripts/02_GetData.R (sourced; clears the workspace)
#          results/tables/diagnostics_03g.csv (03h; model values for reference, optional)
# Outputs: results/tables/16_tt_rst_fish_checks/16_rates_by_group.csv
#          results/tables/16_tt_rst_fish_checks/16_reach1_vs_reach2_by_class.csv
#          results/tables/16_tt_rst_fish_checks/16_between_group_sd.csv
#          results/tables/16_tt_rst_fish_checks/16_within_group_spread.csv
#          results/figures/16_tt_rst_fish_checks/16_rates_by_reach.png
#          results/figures/16_tt_rst_fish_checks/16_within_group_spread.png
# Seconds. Restart R first.
# ============================================================

source(here::here("scripts", "02_GetData.R"))   # d_Sac_sort, TTfR, ReachKM_ind, rgind

library(tidyverse)
library(here)

fig_dir <- here("results", "figures", "16_tt_rst_fish_checks")
tab_dir <- here("results", "tables", "16_tt_rst_fish_checks")
walk(c(fig_dir, tab_dir), \(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

min_n <- 3                                   # fish per group for group-level summaries
theme_set(theme_bw(base_size = 12))

# Model values from 03h, for reference lines and columns (optional)
diag_file <- here("results", "tables", "diagnostics_03g.csv")
diag <- if (file.exists(diag_file)) read_csv(diag_file, show_col_types = FALSE) else NULL
model_val <- \(p) if (is.null(diag) || !p %in% diag$parameter) NA_real_ else diag$median[diag$parameter == p]
pro_sd <- model_val("Pro_sd")

# ---- Fish table ------------------------------------------------------------------
# Origin: wild if fish_type or StudyID says so; check the printed table below.
fish <- tibble(
  rg = rgind, study = d_Sac_sort$StudyID, site = d_Sac_sort$rel_site,
  fish_type = d_Sac_sort$fish_type, year = d_Sac_sort$year,
  km1 = ReachKM_ind[[1]], km2 = ReachKM_ind[[2]], km3 = ReachKM_ind[[3]],
  t_wood = TTfR$TTfR1, t_butte = TTfR$TTfR2, t_sac = TTfR$TTfR3       # days from release (02, A1 order)
) |>
  mutate(origin = if_else(str_detect(str_to_lower(paste(fish_type, study)), "wild"), "wild", "hatchery"),
         class  = paste(site, origin))

cat("\nOrigin assignment (check it):\n")
fish |> count(site, origin, fish_type, name = "fish") |> as.data.frame() |> print()

# ---- 1. Travel rate by reach ------------------------------------------------------
reach_lv <- c("1 Release-Woodson", "2 Woodson-Butte", "3 Butte-Sacramento")
rates <- fish |>
  transmute(rg, study, site, origin, class, year,
            `1 Release-Woodson`  = t_wood / km1 * 100,
            `2 Woodson-Butte`    = (t_butte - t_wood) / km2 * 100,
            `3 Butte-Sacramento` = (t_sac - t_butte) / km3 * 100) |>
  pivot_longer(all_of(reach_lv), names_to = "reach", values_to = "rate") |>
  filter(is.finite(rate), rate > 0) |>
  mutate(reach = factor(reach, levels = reach_lv))

grp_rates <- rates |>
  summarise(n = n(), med_rate = median(rate), .by = c(rg, study, site, origin, class, year, reach)) |>
  filter(n >= min_n)
write_csv(grp_rates |> mutate(med_rate = round(med_rate, 3)), file.path(tab_dir, "16_rates_by_group.csv"))

# Paired within each group: reach-1 rate relative to reach 2
paired <- grp_rates |>
  select(rg, class, reach, med_rate) |>
  pivot_wider(names_from = reach, values_from = med_rate) |>
  mutate(log_r1_r2 = log(`1 Release-Woodson` / `2 Woodson-Butte`)) |>
  filter(is.finite(log_r1_r2))

r12 <- bind_rows(paired, mutate(paired, class = "All")) |>
  summarise(groups = n(),
            median_log_ratio = median(log_r1_r2),
            min_log_ratio = min(log_r1_r2), max_log_ratio = max(log_r1_r2),
            .by = class) |>
  mutate(reach1_over_reach2 = exp(median_log_ratio),
         model_reach1_over_reach2 = exp(model_val("T_bReach[1]") - model_val("T_bReach[2]")),
         across(where(is.double), \(x) round(x, 2))) |>
  arrange(class != "All", class)
write_csv(r12, file.path(tab_dir, "16_reach1_vs_reach2_by_class.csv"))

between_sd <- grp_rates |>
  summarise(groups = n(), sd_log_rate = sd(log(med_rate)), .by = reach) |>
  mutate(model_TTRE_sd = map_dbl(as.integer(reach), \(j) model_val(paste0("TTRE_sd[", j, "]"))),
         across(where(is.double), \(x) round(x, 2))) |>
  arrange(reach)
write_csv(between_sd, file.path(tab_dir, "16_between_group_sd.csv"))

p1 <- grp_rates |>
  ggplot(aes(reach, med_rate, group = rg)) +
  geom_line(colour = "grey60") +
  geom_point(aes(size = n, shape = origin, fill = origin), colour = "black") +
  scale_shape_manual(values = c(hatchery = 21, wild = 24)) +
  scale_fill_manual(values = c(hatchery = "#E69F00", wild = "#0072B2")) +
  scale_y_log10() +
  facet_wrap(~ site) +
  labs(x = NULL, y = "Group median travel rate (days per 100 km, log)", size = "Fish",
       shape = "Origin", fill = "Origin",
       caption = "Lines join the same release group across reaches. Wild fish at Mill, Deer and RBDD were caught and released at the RST.") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.caption = element_text(hjust = 0))
ggsave(file.path(fig_dir, "16_rates_by_reach.png"), p1, width = 11, height = 7.5, dpi = 300)

# ---- 2. Fish-level spread within groups ------------------------------------------
seg_lv <- c("Release to Woodson", "Release to Sacramento", "Woodson to Sacramento", "Butte to Sacramento")
spread <- fish |>
  filter(is.finite(t_wood), is.finite(t_butte), is.finite(t_sac),
         t_wood > 0, t_butte > t_wood, t_sac > t_butte) |>       # same fish in every segment
  transmute(rg, class,
            `Release to Woodson`    = log(t_wood),
            `Release to Sacramento` = log(t_sac),
            `Woodson to Sacramento` = log(t_sac - t_wood),
            `Butte to Sacramento`   = log(t_sac - t_butte)) |>
  pivot_longer(all_of(seg_lv), names_to = "segment", values_to = "log_tt")

pooled_sd <- function(d) {       # pooled within-group SD of log travel time
  g <- d |> summarise(n = n(), v = var(log_tt), .by = rg) |> filter(n >= min_n)
  tibble(groups = nrow(g), fish = sum(g$n),
         sd_within = sqrt(sum((g$n - 1) * g$v) / sum(g$n - 1)))
}

spread_tab <- bind_rows(mutate(spread, class = "All"), spread) |>
  nest(.by = c(class, segment)) |>
  mutate(s = map(data, pooled_sd)) |>
  select(-data) |>
  unnest(s) |>
  mutate(segment = factor(segment, levels = seg_lv),
         model_Pro_sd = pro_sd,
         across(where(is.double), \(x) round(x, 3))) |>
  arrange(class != "All", class, segment)
write_csv(spread_tab, file.path(tab_dir, "16_within_group_spread.csv"))

p2 <- spread_tab |>
  filter(groups > 0) |>
  mutate(is_all = class == "All") |>
  ggplot(aes(segment, sd_within, group = class)) +
  { if (is.finite(pro_sd)) geom_hline(yintercept = pro_sd, linetype = "dashed", colour = "#D55E00") } +
  geom_line(aes(linewidth = is_all), colour = "grey50") +
  geom_point(aes(shape = class, size = is_all), fill = "#0072B2", colour = "black") +
  scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.4), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.2), guide = "none") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 3, 4, 8)[seq_len(n_distinct(spread_tab$class))]) +
  labs(x = NULL, y = "Pooled within-group SD of log travel time", shape = NULL,
       caption = paste0("Same fish (detected at Woodson, Butte and Sacramento) in every segment. ",
                        "Dashed: model Pro_sd (", pro_sd, "). Large points: all groups.")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), plot.caption = element_text(hjust = 0))
ggsave(file.path(fig_dir, "16_within_group_spread.png"), p2, width = 10, height = 6, dpi = 300)

# ---- Console ---------------------------------------------------------------------
cat("\n=== 1a. Reach-1 rate relative to reach 2, within groups (>1 = reach 1 slower) ===\n")
print(as.data.frame(r12))
cat("\n=== 1b. Between-group SD of log travel rate (raw; includes size and discharge) ===\n")
print(as.data.frame(between_sd))
cat("\n=== 2. Pooled within-group SD of log travel time, same fish ===\n")
print(as.data.frame(spread_tab))
cat("\nFigures:", fig_dir, "\nTables: ", tab_dir, "\n")
