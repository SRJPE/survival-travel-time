# Compare fork-length based survival and travel-time predictions from:
#   1) scripts/CWT.stan
#   2) scripts/CovIndCont.stan, treated here as the CJS model
#
# Population-level predictions are computed with random effects set to 0.
# CWT release locations are combined into one curve. CJS/CovInd uses only the
# Sacramento prediction curve.
# The fork-length range is min(CWT average fork length) to 100 mm.
# Max-flow predictions hold fork length at the median observed fork length.

rm(list = ls())

library(rstan)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(here)

# Change this block when comparing different model outputs.
compare_config <- list(
  cwt_fit_file = here("results", "fit_CWT.Rdata"),
  cjs_fit_file = here("results", "fit_CovIndCont_MaxFlow_FL.Rdata"),
  output_tag = "CWT_CJS"
)

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

table_file <- function(prefix) {
  here("results", "tables", paste0(prefix, "_", compare_config$output_tag, ".csv"))
}

fig_file <- function(prefix) {
  here("results", "figures", paste0(prefix, "_", compare_config$output_tag, ".png"))
}

inv_logit <- function(x) 1 / (1 + exp(-x))

summarise_draws <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    q2.5 = unname(quantile(x, 0.025, na.rm = TRUE)),
    q50 = unname(quantile(x, 0.5, na.rm = TRUE)),
    q97.5 = unname(quantile(x, 0.975, na.rm = TRUE))
  )
}

load_rdata <- function(path) {
  env <- new.env(parent = globalenv())
  objects <- load(path, envir = env)
  list(env = env, objects = objects)
}

# Load CWT model objects: fit_cwt, cwt_groups, cwt_data.
cwt_loaded <- load_rdata(compare_config$cwt_fit_file)
fit_cwt <- get("fit_cwt", envir = cwt_loaded$env)
cwt_groups <- get("cwt_groups", envir = cwt_loaded$env)

cwt_post <- rstan::extract(
  fit_cwt,
  pars = c("S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz")
)

cwt_size_mu <- mean(cwt_groups$avg_length, na.rm = TRUE)
cwt_size_sd <- sd(cwt_groups$avg_length, na.rm = TRUE)
cwt_flow_mu <- mean(cwt_groups$Maxflowsac, na.rm = TRUE)
cwt_flow_sd <- sd(cwt_groups$Maxflowsac, na.rm = TRUE)
cwt_median_flow <- median(cwt_groups$Maxflowsac, na.rm = TRUE)

cwt_release_distance_km <- c(
  Battle = 300.443064082,
  RBDD = 240.3631
)

# Load CJS/CovInd model.
cjs_loaded <- load_rdata(compare_config$cjs_fit_file)
cjs_fit <- if ("fit" %in% cjs_loaded$objects) {
  get("fit", envir = cjs_loaded$env)
} else {
  stanfit_objects <- cjs_loaded$objects[
    vapply(cjs_loaded$objects, function(x) inherits(get(x, envir = cjs_loaded$env), "stanfit"), logical(1))
  ]
  if (length(stanfit_objects) != 1) {
    stop("Could not identify a unique CJS stanfit object in ", compare_config$cjs_fit_file)
  }
  get(stanfit_objects, envir = cjs_loaded$env)
}

# Rebuild CovInd data objects in a separate environment because 02_GetData.R
# clears its target environment.
data_env <- new.env(parent = globalenv())
sys.source(here("scripts", "02_GetData.R"), envir = data_env)

FL <- data_env$FL
MaxflowSac <- data_env$MaxflowSac
Rmult <- data_env$Rmult
ReachKM_ind <- data_env$ReachKM_ind

cjs_post <- rstan::extract(
  cjs_fit,
  pars = c("S_bReach", "T_bReach", "S_bCov", "TT_bCov", "S_bSz", "T_bSz")
)

cjs_size_mu <- mean(FL, na.rm = TRUE)
cjs_size_sd <- sd(FL, na.rm = TRUE)

cjs_flow_mu <- mean(MaxflowSac, na.rm = TRUE)
cjs_flow_sd <- sd(MaxflowSac, na.rm = TRUE)
cjs_median_flow <- median(MaxflowSac, na.rm = TRUE) 

# Mainstem CJS is release to Sacramento: reaches 1:3.
cjs_rmult <- colMeans(Rmult[, 1:3, drop = FALSE], na.rm = TRUE)
cjs_reach_km <- colMeans(ReachKM_ind[, 1:3, drop = FALSE], na.rm = TRUE)

shared_median_flow <- median(c(cwt_groups$Maxflowsac, MaxflowSac), na.rm = TRUE)

shared_fl_grid <- seq(
  from = floor(min(cwt_groups$avg_length, na.rm = TRUE)),
  to = 100,
  length.out = 100
)

shared_flow_min <- max(
  min(cwt_groups$Maxflowsac, na.rm = TRUE),
  min(MaxflowSac, na.rm = TRUE)
)
shared_flow_max <- min(
  max(cwt_groups$Maxflowsac, na.rm = TRUE),
  max(MaxflowSac, na.rm = TRUE)
)

if (!is.finite(shared_flow_min) || !is.finite(shared_flow_max) || shared_flow_min >= shared_flow_max) {
  shared_flow_min <- min(c(cwt_groups$Maxflowsac, MaxflowSac), na.rm = TRUE)
  shared_flow_max <- max(c(cwt_groups$Maxflowsac, MaxflowSac), na.rm = TRUE)
}

shared_flow_grid <- seq(
  from = shared_flow_min,
  to = shared_flow_max,
  length.out = 100
)

shared_fl_mm <- median(c(cwt_groups$avg_length, FL), na.rm = TRUE)

# Prediction functions --------------------------------------------------------------
cwt_draws_one <- function(size_mm, area, flow_cfs = cwt_median_flow) {
  reach_km <- unname(cwt_release_distance_km[area])
  size_z <- (size_mm - cwt_size_mu) / cwt_size_sd
  flow_z <- (flow_cfs - cwt_flow_mu) / cwt_flow_sd
  rmult <- reach_km / 100

  surv_draw <- inv_logit(
    cwt_post$S_bReach +
      cwt_post$S_bCov * flow_z +
      cwt_post$S_bSz * size_z
  )^rmult

  tt_draw <- exp(
    cwt_post$T_bReach +
      cwt_post$TT_bCov * flow_z +
      cwt_post$T_bSz * size_z
  ) * reach_km / 100

  tibble(survival = surv_draw, travel_time = tt_draw)
}

cjs_sac_draws_one <- function(size_mm, flow_cfs = cjs_median_flow) {
  size_z <- (size_mm - cjs_size_mu) / cjs_size_sd
  flow_z <- (flow_cfs - cjs_flow_mu) / cjs_flow_sd

  surv_reach_draws <- lapply(seq_along(cjs_rmult), function(j) {
    inv_logit(
      cjs_post$S_bReach[, j] +
        cjs_post$S_bCov * flow_z +
        cjs_post$S_bSz * size_z
    )^cjs_rmult[j]
  })
  surv_draw <- Reduce(`*`, surv_reach_draws)

  tt_reach_draws <- lapply(seq_along(cjs_reach_km), function(j) {
    exp(
      cjs_post$T_bReach[, j] +
        cjs_post$TT_bCov * flow_z +
        cjs_post$T_bSz * size_z
    ) * cjs_reach_km[j] / 100
  })
  tt_draw <- Reduce(`+`, tt_reach_draws)

  tibble(survival = surv_draw, travel_time = tt_draw)
}

combine_location_draws <- function(draw_list) {
  surv_mat <- do.call(cbind, lapply(draw_list, `[[`, "survival"))
  tt_mat <- do.call(cbind, lapply(draw_list, `[[`, "travel_time"))

  list(
    survival = rowMeans(surv_mat, na.rm = TRUE),
    travel_time = rowMeans(tt_mat, na.rm = TRUE)
  )
}

cwt_predict_combined <- function(size_mm, flow_cfs = cwt_median_flow) {
  combined <- combine_location_draws(lapply(names(cwt_release_distance_km), function(area) {
    cwt_draws_one(size_mm, area, flow_cfs = flow_cfs)
  }))

  bind_cols(
    tibble(
      model = "CWT",
      location = "Combined CWT release locations",
      fork_length_mm = size_mm,
      max_flow_cfs = flow_cfs
    ),
    summarise_draws(combined$survival) %>% rename_with(~ paste0("survival_", .x)),
    summarise_draws(combined$travel_time) %>% rename_with(~ paste0("travel_time_", .x))
  )
}

cjs_predict_sacramento <- function(size_mm, flow_cfs = cjs_median_flow) {
  cjs_draws <- cjs_sac_draws_one(size_mm, flow_cfs = flow_cfs)
  bind_cols(
    tibble(
      model = "CJS",
      location = "Sacramento",
      fork_length_mm = size_mm,
      max_flow_cfs = flow_cfs
    ),
    summarise_draws(cjs_draws$survival) %>% rename_with(~ paste0("survival_", .x)),
    summarise_draws(cjs_draws$travel_time) %>% rename_with(~ paste0("travel_time_", .x))
  )
}

predict_over_grid <- function(grid, grid_col, predict_fn) {
  tibble(!!grid_col := grid) %>%
    mutate(pred = map(.data[[grid_col]], predict_fn)) %>%
    select(pred) %>%
    unnest(pred)
}

comparison_long <- function(comparison, x_col) {
  bind_rows(
    comparison %>%
      transmute(
        model,
        x = .data[[x_col]],
        metric = "Survival",
        mean = survival_mean,
        q2.5 = survival_q2.5,
        q97.5 = survival_q97.5
      ),
    comparison %>%
      transmute(
        model,
        x = .data[[x_col]],
        metric = "Travel time",
        mean = travel_time_mean,
        q2.5 = travel_time_q2.5,
        q97.5 = travel_time_q97.5
      )
  ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("Survival", "Travel time"),
        labels = c("Survival posterior prediction", "Travel time posterior prediction")
      )
    )
}

# Build comparison table ------------------------------------------------------------
cwt_comparison <- predict_over_grid(
  shared_fl_grid,
  "fork_length_mm",
  ~ cwt_predict_combined(.x, flow_cfs = shared_median_flow)
)

cjs_comparison <- predict_over_grid(
  shared_fl_grid,
  "fork_length_mm",
  ~ cjs_predict_sacramento(.x, flow_cfs = shared_median_flow)
)

fl_comparison <- bind_rows(cwt_comparison, cjs_comparison) %>%
  mutate(model = factor(model, levels = c("CWT", "CJS")))

write.csv(
  fl_comparison,
  table_file("fl_prediction_comparison"),
  row.names = FALSE
)

cwt_flow_comparison <- predict_over_grid(
  shared_flow_grid,
  "max_flow_cfs",
  ~ cwt_predict_combined(shared_fl_mm, flow_cfs = .x)
)

cjs_flow_comparison <- predict_over_grid(
  shared_flow_grid,
  "max_flow_cfs",
  ~ cjs_predict_sacramento(shared_fl_mm, flow_cfs = .x)
)

flow_comparison <- bind_rows(cwt_flow_comparison, cjs_flow_comparison) %>%
  mutate(model = factor(model, levels = c("CWT", "CJS")))

write.csv(
  flow_comparison,
  table_file("max_flow_prediction_comparison"),
  row.names = FALSE
)

# Plots -----------------------------------------------------------------------------
fl_comparison_long <- comparison_long(fl_comparison, "fork_length_mm")

plot_fl_comparison <- ggplot(
  fl_comparison_long,
  aes(x = x, y = mean, color = model, fill = model)
) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_x_continuous(
    breaks = seq(floor(min(shared_fl_grid) / 5) * 5, 100, by = 5)
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  labs(
    x = "Fork length (mm)",
    y = NULL,
    color = "Model",
    fill = "Model",
    title = "CWT vs CJS fork-length predictions",
    subtitle = paste0(
      "CWT release locations combined; CJS Sacramento only; FL range ",
      round(min(shared_fl_grid), 1), "-", round(max(shared_fl_grid), 1),
      " mm; Red Bluff max flow fixed at ", scales::comma(round(shared_median_flow, 0)), " cfs"
    )
  )

ggsave(
  fig_file("fl_survival_travel_time_comparison"),
  plot = plot_fl_comparison,
  width = 9,
  height = 9,
  dpi = 350
)

flow_comparison_long <- comparison_long(flow_comparison, "max_flow_cfs")

plot_flow_comparison <- ggplot(
  flow_comparison_long,
  aes(x = x, y = mean, color = model, fill = model)
) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    labels = scales::label_comma()
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  labs(
    x = "Red Bluff max flow (cfs)",
    y = NULL,
    color = "Model",
    fill = "Model",
    title = "CWT vs CJS Red Bluff flow predictions",
    subtitle = paste0(
      "CWT release locations combined; CJS Sacramento only; FL fixed at ",
      round(shared_fl_mm, 1), " mm"
    )
  )

ggsave(
  fig_file("max_flow_survival_travel_time_comparison"),
  plot = plot_flow_comparison,
  width = 9,
  height = 9,
  dpi = 350
)
