# ============================================================
# 17_tt_trib_delay.R
# A tributary delay for fish starting at a tributary RST, estimated from the
# wild fish caught and released at the Mill and Deer RSTs, and a fourth
# option for RST to Delta entry travel time built on it.
#
# Background (16): within release groups, reach-1 and reach-2 travel rates are
# about equal for hatchery fish and for wild fish released at the mainstem RBDD
# RST. Wild fish released at the Mill and Deer RSTs are 3 to 9 times slower in
# reach 1, over only 18 to 26 km: a fixed delay in or near the creek, not a
# slower river. The 03g model fits that as a rate x distance, which the
# forecast then scales to 93 to 125 km at Upper Clear and Battle.
#
# Part 1, delay: for each wild Mill / Deer fish detected at Woodson,
#   delay = (release to Woodson time) - (reach-1 km x its group's reach-2 rate)
# The group's reach-2 rate is the group median of (Woodson to Butte time /
# 88 km); groups with fewer than min_r2 such fish use the median over all
# wild Mill / Deer fish. Delays are summarized per fish and per group.
#
# Part 2, option 4 (proposed structure):
#   travel time = delay + reach-2 rate over RST to Woodson + reaches 2 and 3,
#   with the reach-2 group SD (as option 2), and fish-level spread from 16
#   (Woodson to Sacramento, about 0.42) instead of Pro_sd (0.54).
#   The delay is a group-level quantity: each posterior draw takes one group
#   median delay at random (negative values set to 0). Tributary RSTs (Upper
#   Clear, Battle, Mill, Deer) get the delay; mainstem RSTs (Redding, Red Bluff)
#   do not. Options 1 to 3 are as in 15, for comparison.
#
# Caution: the delay rests on a handful of small groups (counts are printed).
# Treat it as a flagged assumption for the team, not an estimate to rely on.
#
# Inputs : scripts/02_GetData.R, results/handoff/stm.Rdata (after 05d and 11),
#          results/tables/16_tt_rst_fish_checks/16_within_group_spread.csv (from 16)
# Outputs: results/tables/17_tt_trib_delay/17_delay_by_fish.csv
#          results/tables/17_tt_trib_delay/17_delay_by_group.csv
#          results/tables/17_tt_trib_delay/17_tt_by_option.csv
#          results/figures/17_tt_trib_delay/17_delay_by_group.png
#          results/figures/17_tt_trib_delay/17_tt_by_option.png
# Restart R first.
# ============================================================

library(tidyverse)
library(here)

stm_file    <- here("results", "handoff", "stm.Rdata")
spread_file <- here("results", "tables", "16_tt_rst_fish_checks", "16_within_group_spread.csv")
fig_dir <- here("results", "figures", "17_tt_trib_delay")
tab_dir <- here("results", "tables", "17_tt_trib_delay")
walk(c(fig_dir, tab_dir), \(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

trib_sites     <- c("Mill", "Deer")      # wild fish released at these RSTs
min_r2         <- 2                      # fish needed for a group's own reach-2 rate
min_fish_delay <- 2                      # fish needed for a group median delay
fl_clamp <- 70; z_flow <- 0
fl_grid  <- seq(30, 150, by = 5); fl_table <- c(70, 90, 110)
probs    <- c(0.1, 0.5, 0.9)
seed     <- 1234
theme_set(theme_bw(base_size = 12))

rst_names <- c("Upper Clear Ck", "Battle Ck", "Mill Ck", "Deer Ck", "Redding", "Red Bluff")
rst_trib  <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)    # tributary RST gets the delay
opt_lab <- c(`1` = "1: reach-1 rate and group SD (current stm)",
             `2` = "2: reach-2 rate and group SD",
             `3` = "3: reach-1 rate, reach-2 group SD",
             `4` = "4: tributary delay + reach-2 rate and group SD, fish SD from data")
opt_col <- c(`1` = "#000000", `2` = "#0072B2", `3` = "#E69F00", `4` = "#CC79A7")   # Okabe-Ito
opt_lty <- c(`1` = "solid", `2` = "longdash", `3` = "dotdash", `4` = "twodash")

# ---- Data from 02 ---------------------------------------------------------------
d <- new.env()
sys.source(here("scripts", "02_GetData.R"), envir = d)
fish <- tibble(rg = d$rgind, study = d$d_Sac_sort$StudyID, site = d$d_Sac_sort$rel_site,
               fish_type = d$d_Sac_sort$fish_type, year = d$d_Sac_sort$year,
               km1 = d$ReachKM_ind[[1]], km2 = d$ReachKM_ind[[2]],
               t_wood = d$TTfR$TTfR1, t_butte = d$TTfR$TTfR2) |>
  mutate(origin = if_else(str_detect(str_to_lower(paste(fish_type, study)), "wild"), "wild", "hatchery"))
km_rst <- d$dist_rstwoodson
km_23  <- d$ReachKM[2:3]
rm(d); invisible(gc())

# ---- Part 1: delay ----------------------------------------------------------------
trib <- fish |>
  filter(site %in% trib_sites, origin == "wild") |>
  mutate(rate2 = if_else(is.finite(t_butte) & is.finite(t_wood) & t_butte > t_wood,
                         (t_butte - t_wood) / km2 * 100, NA_real_))

pooled_r2 <- median(trib$rate2, na.rm = TRUE)
grp_r2 <- trib |>
  summarise(n_r2 = sum(is.finite(rate2)), grp_rate2 = median(rate2, na.rm = TRUE), .by = rg)

delay_fish <- trib |>
  filter(is.finite(t_wood), t_wood > 0) |>
  left_join(grp_r2, by = "rg") |>
  mutate(r2_used   = if_else(n_r2 >= min_r2, grp_rate2, pooled_r2),
         r2_source = if_else(n_r2 >= min_r2, "own group", "pooled Mill/Deer"),
         migrate   = km1 * r2_used / 100,
         delay     = t_wood - migrate)
write_csv(delay_fish |> mutate(across(where(is.double), \(x) round(x, 3))),
          file.path(tab_dir, "17_delay_by_fish.csv"))

delay_grp <- delay_fish |>
  summarise(n = n(), km1 = first(km1), r2_used = first(r2_used), r2_source = first(r2_source),
            median_t_wood = median(t_wood), median_migrate = median(migrate),
            median_delay = median(delay), .by = c(rg, study, site, year)) |>
  mutate(used = n >= min_fish_delay) |>
  arrange(site, year)
write_csv(delay_grp |> mutate(across(where(is.double), \(x) round(x, 2))),
          file.path(tab_dir, "17_delay_by_group.csv"))

grp_delays <- pmax(delay_grp$median_delay[delay_grp$used], 0)
if (length(grp_delays) < 3) warning("Fewer than 3 groups for the delay; treat option 4 as illustrative only")

p1 <- ggplot(delay_fish, aes(factor(year), delay)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_jitter(width = 0.15, height = 0, shape = 21, fill = "#56B4E9", colour = "black", alpha = 0.7) +
  geom_point(data = delay_grp, aes(y = median_delay, shape = used), size = 4, fill = "#0072B2") +
  scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 4),
                     labels = c(`TRUE` = "group median (used)", `FALSE` = paste0("group median (< ", min_fish_delay, " fish)"))) +
  facet_wrap(~ site, scales = "free_x") +
  labs(x = "Release year", y = "Delay before reaching Woodson at the reach-2 rate (days)", shape = NULL,
       caption = "Delay = release to Woodson time minus reach-1 distance at the group's Woodson-Butte rate. Wild fish released at the Mill and Deer RSTs.") +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(file.path(fig_dir, "17_delay_by_group.png"), p1, width = 9, height = 5.5, dpi = 300)

# ---- Fish-level SD from 16 ---------------------------------------------------------
sd_fish <- if (file.exists(spread_file)) {
  read_csv(spread_file, show_col_types = FALSE) |>
    filter(class == "All", segment == "Woodson to Sacramento") |>
    pull(sd_within)
} else NA_real_
if (!is.finite(sd_fish)) stop("Run 16 first: need the Woodson to Sacramento within-group SD")

# ---- Part 2: options 1 to 4 --------------------------------------------------------
e <- new.env(); load(stm_file, envir = e); s <- e$stm; rm(e); invisible(gc())
need <- c("T_bReach", "TTRE_sd", "T_bSz", "TT_bCov", "musz", "sdsz", "Pro_sd")
stopifnot("stm is missing objects (run 05d, then 11)" = all(need %in% names(s)))

nd <- nrow(s$T_bReach); T_bSz <- as.vector(s$T_bSz); TT_bCov <- as.vector(s$TT_bCov)
set.seed(seed)
u      <- matrix(rnorm(nd * 3), nd, 3)        # new-group draws, shared across options
u_fish <- rnorm(nd)                           # fish-level draw, shared across options
delay_draw <- sample(grp_delays, nd, replace = TRUE)

z_of <- \(fl) (pmax(fl, fl_clamp) - s$musz) / s$sdsz

tt_group <- function(r, z_sz, rate1, sd1, add_delay) {
  lin <- \(j_rate, j_sd, uj) s$T_bReach[, j_rate] + TT_bCov * z_flow + T_bSz * z_sz + s$TTRE_sd[, j_sd] * uj
  (if (add_delay && rst_trib[r]) delay_draw else 0) +
    exp(lin(rate1, sd1, u[, 1])) * km_rst[r] / 100 +
    exp(lin(2, 2, u[, 2])) * km_23[1] / 100 +
    exp(lin(3, 3, u[, 3])) * km_23[2] / 100
}

qt3 <- function(x) {
  v <- quantile(x, probs, names = FALSE, na.rm = TRUE)
  tibble(lo = v[1], med = v[2], hi = v[3])
}

options_tbl <- tribble(~option, ~rate1, ~sd1, ~add_delay, ~fish_sd_src,
                       "1",     1,      1,    FALSE,      "Pro_sd",
                       "2",     2,      2,    FALSE,      "Pro_sd",
                       "3",     1,      2,    FALSE,      "Pro_sd",
                       "4",     2,      2,    TRUE,       "data")

res <- expand_grid(rst = seq_along(rst_names), options_tbl, fl = sort(union(fl_grid, fl_table))) |>
  mutate(draws = pmap(list(rst, fl, rate1, sd1, add_delay), \(r, f, r1, s1, a) tt_group(r, z_of(f), r1, s1, a)),
         grp   = map(draws, qt3),
         fish  = map2(draws, fish_sd_src, \(x, src) qt3(x * exp((if (src == "data") sd_fish else s$Pro_sd) * u_fish)))) |>
  select(-draws, -rate1, -sd1, -add_delay) |>
  unnest_wider(grp, names_sep = "_") |>
  unnest_wider(fish, names_sep = "_") |>
  mutate(rst = factor(rst_names[rst], levels = rst_names)) |>
  mutate(ratio_to_opt1 = grp_med / grp_med[option == "1"], .by = c(rst, fl))

tab <- res |>
  filter(fl %in% fl_table) |>
  transmute(rst, fork_length_mm = fl, option,
            group_med = round(grp_med, 1), group_lo80 = round(grp_lo, 1), group_hi80 = round(grp_hi, 1),
            fish_sd = fish_sd_src, fish_lo80 = round(fish_lo, 1), fish_hi80 = round(fish_hi, 1),
            ratio_to_opt1 = round(ratio_to_opt1, 2)) |>
  arrange(rst, fork_length_mm, option)
write_csv(tab, file.path(tab_dir, "17_tt_by_option.csv"))

p2 <- res |>
  filter(fl %in% fl_grid) |>
  ggplot(aes(fl, grp_med, colour = option, fill = option, linetype = option)) +
  geom_ribbon(aes(ymin = grp_lo, ymax = grp_hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = fl_clamp, linetype = "dotted") +
  facet_wrap(~ rst, scales = "free_y") +
  scale_colour_manual(values = opt_col, labels = opt_lab) +
  scale_fill_manual(values = opt_col, labels = opt_lab) +
  scale_linetype_manual(values = opt_lty, labels = opt_lab) +
  labs(x = "Fork length (mm)", y = "RST to Delta entry travel time (days)",
       colour = NULL, fill = NULL, linetype = NULL,
       title = "RST to Delta entry travel time: current structure vs a tributary delay",
       subtitle = "New release group, mean discharge; median and 80% interval of the group-level travel time",
       caption = paste0("Option 4 delay: a group median delay drawn from ", length(grp_delays),
                        " Mill/Deer wild groups, added at tributary RSTs only (not Redding, Red Bluff).\n",
                        "Dotted: fork-length clamp at ", fl_clamp, " mm.")) +
  guides(colour = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
ggsave(file.path(fig_dir, "17_tt_by_option.png"), p2, width = 11, height = 8, dpi = 300)

# ---- Console -------------------------------------------------------------------------
cat("\n=== Part 1. Delay by group (wild fish released at Mill and Deer RSTs) ===\n")
cat("Pooled Mill/Deer reach-2 rate:", round(pooled_r2, 2), "days per 100 km\n")
print(as.data.frame(delay_grp |> mutate(across(where(is.double), \(x) round(x, 2)))))
cat("\nGroup delays used (days, negatives set to 0):", round(grp_delays, 2), "\n")
cat("Fish-level SD for option 4 (16, Woodson to Sacramento):", sd_fish, "; Pro_sd:", round(median(s$Pro_sd), 3), "\n")
cat("\n=== Part 2. RST to Delta entry travel time (days) by option ===\n")
print(as.data.frame(tab))
cat("\nFigures:", fig_dir, "\nTables: ", tab_dir, "\n")
