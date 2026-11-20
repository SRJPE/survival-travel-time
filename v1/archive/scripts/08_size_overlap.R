# ============================================================
# 08_size_overlap.R
#
# How different in size are CWT and acoustic telemetry fish, and does the
# overlap change among years? Data only, no model fitting (about 1 minute).
#
# Inputs (project data/ folder, as used by 02 and 06o):
#   data/Sac_data.csv       telemetry, mainstem releases (fish_length, year)
#   data/FeaBut_data.csv    telemetry, Butte and Feather releases
#   data/hatchery_relrec.csv CWT releases and KL recoveries
# CWT release groups are built exactly as in 06o (Battle Creek and RBDD
# releases, group 873 dropped), but for ALL years, including those without
# a BT-SPAS-X fit at KL.
#
# CWT sizes: hatchery_relrec has one mean fork length per release group
# (avg_length), not individual fish lengths at release. Fish-level CWT
# lengths exist only for recovered fish, measured at KL (after some growth).
#
# Outputs (results/figures, results/tables):
#   fig_size_by_year.png         telemetry fish (boxes) and CWT group means (points) by year
#   fig_size_distributions.png   pooled size distributions, all years and 2013-2024
#   fig_size_overlap_by_year.png share of CWT fish released in groups as large as telemetry fish
#   size_by_year.csv             the numbers behind the figures
# ============================================================

library(tidyverse)
library(here)

dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
oi <- c(tel_main = "#0072B2", tel_trib = "#56B4E9", cwt_battle = "#E69F00",
        cwt_rbdd = "#D55E00", cwt_kl = "#999999")

# ---- telemetry fish ----------------------------------------------------------------------
tel <- bind_rows(
  read.csv(here("data", "Sac_data.csv")) |> transmute(source = "tel_main", year, FL = fish_length),
  read.csv(here("data", "FeaBut_data.csv")) |> transmute(source = "tel_trib", year, FL = fish_length)) |>
  filter(is.finite(FL))
tel_q <- quantile(tel$FL, c(0.1, 0.5, 0.9))
message("Telemetry fish: ", nrow(tel), "; pooled 10th, 50th, 90th percentile FL: ",
        paste(tel_q, collapse = ", "), " mm")

# ---- CWT release groups and recovered fish (as 06o) ----------------------------------------
release_distance_km <- c(Battle = 300.443064082, RBDD = 240.3631)
rr <- read.csv(here("data", "hatchery_relrec.csv"), stringsAsFactors = FALSE) |>
  filter(!is.na(relloc_area), relloc_area %in% names(release_distance_km), release_group_id != 873)
cwt <- rr |>
  group_by(release_group_id) |>
  summarise(relloc_area = first(relloc_area), year = first(year),
            Nrel = first(group_total_marked_N), FL = first(avg_length), .groups = "drop") |>
  filter(is.finite(Nrel), Nrel > 0, is.finite(FL)) |>
  mutate(source = if_else(relloc_area == "Battle", "cwt_battle", "cwt_rbdd"))
cwt_kl <- rr |>
  filter(!is.na(tag_code), is.finite(forklength)) |>
  distinct(release_group_id, tag_code, date, forklength, .keep_all = TRUE) |>
  transmute(source = "cwt_kl", year, FL = forklength)
message("CWT groups: ", nrow(cwt), " (", min(cwt$year), " to ", max(cwt$year), "); recovered fish with a length: ", nrow(cwt_kl))

# ---- by-year summary ------------------------------------------------------------------------
wq <- function(x, w, p) { o <- order(x); x <- x[o]; w <- cumsum(w[o]) / sum(w); x[which(w >= p)[1]] }
by_year <- full_join(
  tel |> group_by(year) |>
    summarise(tel_fish = n(), tel_q10 = quantile(FL, 0.1), tel_median = median(FL), tel_q90 = quantile(FL, 0.9),
              .groups = "drop"),
  cwt |> group_by(year) |>
    summarise(cwt_groups = n(), cwt_fish = sum(Nrel),
              cwt_min_group = min(FL), cwt_median_fish = wq(FL, Nrel, 0.5), cwt_max_group = max(FL),
              # share of CWT fish released in groups whose mean FL reaches the telemetry 10th percentile
              cwt_share_ge_tel_q10 = sum(Nrel[FL >= tel_q[[1]]]) / sum(Nrel),
              cwt_share_ge_70 = sum(Nrel[FL >= 70]) / sum(Nrel), .groups = "drop"),
  by = "year") |>
  left_join(cwt_kl |> group_by(year) |> summarise(cwt_kl_fish = n(), cwt_kl_median = median(FL), .groups = "drop"),
            by = "year") |>
  arrange(year) |>
  mutate(across(where(is.double), \(x) round(x, 2)))
write.csv(by_year, here("results", "tables", "size_by_year.csv"), row.names = FALSE)
print(as.data.frame(by_year))

# ---- figure 1: size by year -------------------------------------------------------------------
lab <- c(tel_main = "Telemetry, mainstem releases (fish)", tel_trib = "Telemetry, Butte and Feather releases (fish)",
         cwt_battle = "CWT, Battle Creek releases (group mean)", cwt_rbdd = "CWT, RBDD releases (group mean)")
tel_box <- tel |>
  group_by(source, year) |>
  summarise(q10 = quantile(FL, 0.1), q25 = quantile(FL, 0.25), q50 = median(FL),
            q75 = quantile(FL, 0.75), q90 = quantile(FL, 0.9), .groups = "drop") |>
  mutate(x = year + if_else(source == "tel_main", -0.18, 0.18))
p1 <- ggplot() +
  geom_hline(yintercept = 70, linetype = "dashed", colour = "grey40") +
  geom_boxplot(data = tel_box, aes(x = x, ymin = q10, lower = q25, middle = q50, upper = q75, ymax = q90,
                                   group = interaction(source, year), fill = source),
               stat = "identity", width = 0.32, colour = "grey20", linewidth = 0.3) +
  geom_point(data = cwt, aes(year, FL, colour = source, shape = source, size = Nrel),
             position = position_jitter(width = 0.15, height = 0, seed = 1), alpha = 0.7) +
  scale_fill_manual(values = oi, labels = lab, name = NULL) +
  scale_colour_manual(values = oi, labels = lab, name = NULL) +
  scale_shape_manual(values = c(cwt_battle = 16, cwt_rbdd = 17), labels = lab, name = NULL) +
  scale_size_area(max_size = 3, guide = "none") +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(x = NULL, y = "Fork length (mm)",
       title = "Fish size by year: acoustic telemetry fish and CWT release groups",
       subtitle = "Boxes = telemetry fish (10th, 25th, 50th, 75th, 90th percentiles)\nPoints = CWT group mean length (size = number released); dashed line = 70 mm") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "fig_size_by_year.png"), p1, width = 10, height = 6, dpi = 200, bg = "white")

# ---- figure 2: pooled distributions -----------------------------------------------------------
lab2 <- c(tel = "Telemetry fish at release", cwt_rel = "CWT fish at release (group means, weighted by number released)",
          cwt_kl = "CWT fish recovered at KL (measured at KL)")
dist <- bind_rows(
  tel |> transmute(type = "tel", year, FL, w = 1),
  cwt |> transmute(type = "cwt_rel", year, FL, w = Nrel),
  cwt_kl |> transmute(type = "cwt_kl", year, FL, w = 1)) |>
  group_by(type) |> mutate(w = w / sum(w)) |> ungroup()
dist2 <- bind_rows(dist |> mutate(period = "All years"),
                   dist |> filter(year >= 2013) |> group_by(type) |> mutate(w = w / sum(w)) |> ungroup() |>
                     mutate(period = "2013 to 2024 (years with telemetry)"))  |>
  mutate(period = factor(period, levels = c("All years", "2013 to 2024 (years with telemetry)")))
p2 <- ggplot(dist2, aes(FL, weight = w, fill = type, colour = type)) +
  geom_density(alpha = 0.35, adjust = 1.2, linewidth = 0.5) +
  geom_vline(xintercept = 70, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ period, ncol = 1) +
  scale_fill_manual(values = c(tel = "#0072B2", cwt_rel = "#E69F00", cwt_kl = "#999999"), labels = lab2, name = NULL) +
  scale_colour_manual(values = c(tel = "#0072B2", cwt_rel = "#E69F00", cwt_kl = "#999999"), labels = lab2, name = NULL) +
  coord_cartesian(xlim = c(35, 150)) +
  labs(x = "Fork length (mm)", y = "Density",
       title = "How much do CWT and telemetry fish sizes overlap?",
       subtitle = "CWT release sizes are group means, so the CWT release curve is narrower\nthan the real spread of fish") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1))
ggsave(here("results", "figures", "fig_size_distributions.png"), p2, width = 8, height = 7, dpi = 200, bg = "white")

# ---- figure 3: overlap by year ------------------------------------------------------------------
ov <- by_year |>
  filter(is.finite(cwt_share_ge_70)) |>
  select(year, cwt_share_ge_70, cwt_share_ge_tel_q10) |>
  pivot_longer(-year, names_to = "threshold", values_to = "share") |>
  mutate(threshold = recode(threshold,
                            cwt_share_ge_70 = "Group mean 70 mm or more",
                            cwt_share_ge_tel_q10 = paste0("Group mean ", tel_q[[1]], " mm or more (telemetry 10th percentile)")))
p3 <- ggplot(ov, aes(year, share, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("#E69F00", "#0072B2"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(x = NULL, y = "Share of CWT fish released",
       title = "Share of CWT fish released in groups as large as typical telemetry fish",
       subtitle = "By release year; based on CWT group mean fork length") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(here("results", "figures", "fig_size_overlap_by_year.png"), p3, width = 9, height = 4.5, dpi = 200, bg = "white")
message("Wrote results/tables/size_by_year.csv and results/figures/fig_size_*.png")
