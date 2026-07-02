# CWT survival and travel time model prediction figures

rm(list = ls())

library(tidyverse)
library(here)
library(lubridate)
library(readxl)
library(rstan)

source("scripts/02_GetData.R")

load(here("results", "fit_CWT.Rdata"))

surv_summary <- read.csv(
  here("results", "tables", "cwt_survival_summary.csv"),
  check.names = FALSE
)
tt_summary <- read.csv(
  here("results", "tables", "cwt_travel_time_summary.csv"),
  check.names = FALSE
)

use_size_effect <- cwt_data$UseSizeEffect

release_distance_km <- c(
  Battle = 300.443064082,
  RBDD = 240.3631
)

size_mu <- mean(cwt_groups$avg_length, na.rm = TRUE)
size_sd <- sd(cwt_groups$avg_length, na.rm = TRUE)
flow_mu <- mean(cwt_groups$Maxflowsac, na.rm = TRUE)
flow_sd <- sd(cwt_groups$Maxflowsac, na.rm = TRUE)

# CWT posterior distribution plots -------------------------------------------------------------------
cwt_post <- rstan::extract(
  fit_cwt,
  pars = c("surv_cwt", "pTT", "pKL", "S_bReach", "T_bReach",
           "S_bCov", "TT_bCov", "S_bSz", "T_bSz")
)

plot_group_info <- cwt_groups %>%
  mutate(
    mid_release_date = as.Date(mid_release_date),
    release_group_id = as.factor(release_group_id)
  ) %>%
  select(release_group_id, mid_release_date, year, month, release_location_name,
         relloc_area, Nrel, cwt_recaptures, avg_length, Maxflowsac)

surv_dist_df <- as.data.frame(cwt_post$surv_cwt) %>%
  setNames(as.character(plot_group_info$release_group_id)) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "release_group_id", values_to = "survival") %>%
  left_join(plot_group_info, by = "release_group_id")

tt_dist_df <- as.data.frame(cwt_post$pTT) %>%
  setNames(as.character(plot_group_info$release_group_id)) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "release_group_id", values_to = "travel_time") %>%
  left_join(plot_group_info, by = "release_group_id")

tt_obs_dist_df <- cwt_tt_obs %>%
  mutate(release_group_id = as.factor(release_group_id)) %>%
  left_join(plot_group_info %>% select(release_group_id, relloc_area),
            by = "release_group_id")

(cwt_survival_plot <- ggplot(surv_dist_df,
                             aes(x = survival, fill = relloc_area, color = relloc_area)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    theme_bw() +
    labs(
      x = "Predicted survival",
      y = "Posterior density",
      fill = "Release area",
      color = "Release area",
      title = "CWT Stan Predicted Survival Distribution"
    ))

(cwt_travel_time_plot <- ggplot() +
    geom_density(data = tt_dist_df,
                 aes(x = travel_time, fill = relloc_area, color = relloc_area),
                 alpha = 0.35, linewidth = 0.8) +
    geom_density(data = tt_obs_dist_df,
                 aes(x = cwt_tt, color = relloc_area),
                 linetype = "dashed", linewidth = 0.8,
                 na.rm = TRUE) +
    theme_bw() +
    labs(
      x = "Travel time to Knights Landing (days)",
      y = "Density",
      fill = "Release area",
      color = "Release area",
      title = "CWT Stan Predicted Travel Time Distribution",
      subtitle = "Dashed lines are observed CWT travel-time distributions"
    ))

ggsave(plot = cwt_survival_plot,
       filename = here("results/figures", "cwt_pred_survival.png"),
       width = 11, height = 6)

ggsave(plot = cwt_travel_time_plot,
       filename = here("results/figures", "cwt_pred_travel_time.png"),
       width = 11, height = 6)

# Annual prediction plots ------------------------------------------------------------
summarise_prediction_year <- function(dat) {
  dat %>%
    group_by(year) %>%
    summarise(
      n_release_groups = n(),
      mean = mean(mean, na.rm = TRUE),
      lwr = mean(`2.5%`, na.rm = TRUE),
      med = mean(`50%`, na.rm = TRUE),
      upr = mean(`97.5%`, na.rm = TRUE),
      .groups = "drop"
    )
}

cwt_surv_by_year <- summarise_prediction_year(surv_summary) %>%
  mutate(metric = "Survival")

cwt_tt_by_year <- summarise_prediction_year(tt_summary) %>%
  mutate(metric = "Travel time")

cwt_surv_by_year <- cwt_surv_by_year %>%
  left_join(cdec_wyt, by = "year") %>%
  mutate(dry_wet = factor(dry_wet, levels = c("Dry", "Wet")))

cwt_tt_by_year <- cwt_tt_by_year %>%
  left_join(cdec_wyt, by = "year") %>%
  mutate(dry_wet = factor(dry_wet, levels = c("Dry", "Wet")))

cwt_predictions_by_year <- bind_rows(cwt_surv_by_year, cwt_tt_by_year) %>%
  select(metric, year, WYT, dry_wet, n_release_groups, mean, lwr, med, upr) %>%
  arrange(metric, year)

write.csv(cwt_predictions_by_year,
          here("results/tables", "cwt_predictions_by_year.csv"),
          row.names = FALSE)

all_prediction_years <- sort(unique(cwt_predictions_by_year$year))

(cwt_survival_by_year_plot <- ggplot(cwt_surv_by_year,
                                     aes(x = year, y = mean, ymin = lwr, ymax = upr)) +
    geom_ribbon(alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.8, color = "grey35") +
    geom_point(aes(color = dry_wet), size = 2.2) +
    scale_x_continuous(breaks = all_prediction_years) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c("Dry" = "red", "Wet" = "blue"), na.translate = FALSE) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    labs(
      x = "Year",
      y = "Predicted survival",
      color = "Water year",
      title = "CWT Predicted Survival by Year"
    ))

(cwt_tt_by_year_plot <- ggplot(cwt_tt_by_year,
                               aes(x = year, y = mean, ymin = lwr, ymax = upr)) +
    geom_ribbon(alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.8, color = "grey35") +
    geom_point(aes(color = dry_wet), size = 2.2) +
    scale_x_continuous(breaks = all_prediction_years) +
    scale_color_manual(values = c("Dry" = "red", "Wet" = "blue"), na.translate = FALSE) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    labs(
      x = "Year",
      y = "Predicted travel time to Knights Landing (days)",
      color = "Water year",
      title = "CWT Predicted Travel Time by Year"
    ))

ggsave(plot = cwt_survival_by_year_plot,
       filename = here("results/figures", "cwt_pred_survival_by_year.png"),
       width = 11, height = 6)

ggsave(plot = cwt_tt_by_year_plot,
       filename = here("results/figures", "cwt_pred_travel_time_by_year.png"),
       width = 11, height = 6)

# Covariate prediction plots ---------------------------------------------------------
# Population-level predictions use fixed effects only, with CWT random effects set to 0.
summarise_prediction <- function(size_mm, maxflow, area) {
  reach_km <- unname(release_distance_km[area])
  rmult <- reach_km / 100
  size_z <- (size_mm - size_mu) / size_sd
  flow_z <- (maxflow - flow_mu) / flow_sd
  
  surv_draw <- plogis(cwt_post$S_bReach +
                        cwt_post$S_bCov * flow_z +
                        use_size_effect * cwt_post$S_bSz * size_z)^rmult
  
  tt_draw <- exp(cwt_post$T_bReach +
                   use_size_effect * cwt_post$T_bSz * size_z +
                   cwt_post$TT_bCov * flow_z) * reach_km / 100
  
  tibble(
    survival_mean = mean(surv_draw),
    survival_lwr = quantile(surv_draw, 0.025),
    survival_upr = quantile(surv_draw, 0.975),
    travel_time_mean = mean(tt_draw),
    travel_time_lwr = quantile(tt_draw, 0.025),
    travel_time_upr = quantile(tt_draw, 0.975)
  )
}

length_pred_grid <- expand_grid(
  avg_length = seq(min(cwt_groups$avg_length, na.rm = TRUE),
                   max(cwt_groups$avg_length, na.rm = TRUE),
                   length.out = 100),
  Maxflowsac = median(cwt_groups$Maxflowsac, na.rm = TRUE),
  relloc_area = names(release_distance_km)
) %>%
  mutate(pred = pmap(list(avg_length, Maxflowsac, relloc_area), summarise_prediction)) %>%
  unnest(pred)

flow_pred_grid <- expand_grid(
  avg_length = median(cwt_groups$avg_length, na.rm = TRUE),
  Maxflowsac = seq(min(cwt_groups$Maxflowsac, na.rm = TRUE),
                   max(cwt_groups$Maxflowsac, na.rm = TRUE),
                   length.out = 100),
  relloc_area = names(release_distance_km)
) %>%
  mutate(pred = pmap(list(avg_length, Maxflowsac, relloc_area), summarise_prediction)) %>%
  unnest(pred)

(cwt_survival_by_length_plot <- ggplot(length_pred_grid,
                                       aes(x = avg_length, y = survival_mean,
                                           color = relloc_area, fill = relloc_area)) +
    geom_ribbon(aes(ymin = survival_lwr, ymax = survival_upr), alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_rug(data = cwt_groups, aes(x = avg_length, color = relloc_area), inherit.aes = FALSE,
             sides = "b", alpha = 0.25) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    theme_bw() +
    labs(
      x = "Average fork length (mm)",
      y = "Predicted survival",
      color = "Release area",
      fill = "Release area",
      title = "CWT Predicted Survival by Fish Length",
      subtitle = paste0("Red Bluff max flow held at median CWT value: ",
                        round(median(cwt_groups$Maxflowsac, na.rm = TRUE), 0))
    ))

(cwt_tt_by_length_plot <- ggplot(length_pred_grid,
                                 aes(x = avg_length, y = travel_time_mean,
                                     color = relloc_area, fill = relloc_area)) +
    geom_ribbon(aes(ymin = travel_time_lwr, ymax = travel_time_upr), alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_rug(data = cwt_groups, aes(x = avg_length, color = relloc_area), inherit.aes = FALSE,
             sides = "b", alpha = 0.25) +
    theme_bw() +
    labs(
      x = "Average fork length (mm)",
      y = "Predicted travel time to Knights Landing (days)",
      color = "Release area",
      fill = "Release area",
      title = "CWT Predicted Travel Time by Fish Length",
      subtitle = paste0("Red Bluff max flow held at median CWT value: ",
                        round(median(cwt_groups$Maxflowsac, na.rm = TRUE), 0))
    ))

(cwt_survival_by_flow_plot <- ggplot(flow_pred_grid,
                                     aes(x = Maxflowsac, y = survival_mean,
                                         color = relloc_area, fill = relloc_area)) +
    geom_ribbon(aes(ymin = survival_lwr, ymax = survival_upr), alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_rug(data = cwt_groups, aes(x = Maxflowsac, color = relloc_area), inherit.aes = FALSE,
             sides = "b", alpha = 0.25) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    theme_bw() +
    labs(
      x = "Monthly max flow at Red Bluff",
      y = "Predicted survival",
      color = "Release area",
      fill = "Release area",
      title = "CWT Predicted Survival by Red Bluff Max Flow",
      subtitle = paste0("Fish length held at median CWT value: ",
                        round(median(cwt_groups$avg_length, na.rm = TRUE), 0), " mm")
    ))

(cwt_tt_by_flow_plot <- ggplot(flow_pred_grid,
                               aes(x = Maxflowsac, y = travel_time_mean,
                                   color = relloc_area, fill = relloc_area)) +
    geom_ribbon(aes(ymin = travel_time_lwr, ymax = travel_time_upr), alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_rug(data = cwt_groups, aes(x = Maxflowsac, color = relloc_area), inherit.aes = FALSE,
             sides = "b", alpha = 0.25) +
    theme_bw() +
    labs(
      x = "Monthly max flow at Red Bluff",
      y = "Predicted travel time to Knights Landing (days)",
      color = "Release area",
      fill = "Release area",
      title = "CWT Predicted Travel Time by Red Bluff Max Flow",
      subtitle = paste0("Fish length held at median CWT value: ",
                        round(median(cwt_groups$avg_length, na.rm = TRUE), 0), " mm")
    ))

ggsave(plot = cwt_survival_by_length_plot,
       filename = here("results/figures", "cwt_pred_survival_by_length.png"),
       width = 11, height = 6)

ggsave(plot = cwt_tt_by_length_plot,
       filename = here("results/figures", "cwt_pred_travel_time_by_length.png"),
       width = 11, height = 6)

ggsave(plot = cwt_survival_by_flow_plot,
       filename = here("results/figures", "cwt_pred_survival_by_maxflow.png"),
       width = 11, height = 6)

ggsave(plot = cwt_tt_by_flow_plot,
       filename = here("results/figures", "cwt_pred_travel_time_by_maxflow.png"),
       width = 11, height = 6)
