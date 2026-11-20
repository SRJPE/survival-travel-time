# Summarise the joint CJS + CWT simulation (draft, 30 Sep 2026)
# Reads every fit_*.rds written by 30_run_sim_joint.R and reports, by scenario x model:
# bias, relative bias, 95% interval coverage, interval width, Rhat, divergences, run time.

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

run_dir <- here("results", "sims", "joint_cjs_cwt", "run2")   # "run1" = draft test fits from 30 Sep
out_dir <- file.path(run_dir, "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fits <- list.files(run_dir, pattern = "^fit_.*\\.rds$", full.names = TRUE) |>
  map(readRDS) |>
  map(\(x) rename_with(x, \(n) sub("^posterior::", "", n))) |>
  list_rbind()

model_labels <- c(
  A_cwt_only_constp = "A: CWT only, constant pCap (like 06)",
  B_joint_constp    = "B: Joint, constant pCap + trials",
  C_joint_flowp     = "C: Joint, flow-dependent pCap + trials",
  D_cjs_only        = "D: CJS only (like 03)"
)

# Parameters each model can estimate (others are fixed or not identified in that model)
estimable <- tribble(
  ~model,               ~parameter,
  "A_cwt_only_constp",  c("bQ", "bL", "gQ", "gL", "sd_S", "S_KL", "TT_KL", "a0"),
  "B_joint_constp",     c("bQ", "bL", "delta_AT", "gQ", "gL", "sd_S", "S_KL", "S_KL_tagged", "TT_KL", "a0"),
  "C_joint_flowp",      c("bQ", "bL", "delta_AT", "gQ", "gL", "sd_S", "S_KL", "S_KL_tagged", "TT_KL", "a0", "bp"),
  "D_cjs_only",         c("bQ", "bL", "gQ", "gL", "sd_S", "S_KL", "S_KL_tagged", "TT_KL")
) |>
  unnest(parameter)

res <- fits |>
  semi_join(estimable, by = c("model", "parameter")) |>
  mutate(err = mean - truth,
         covered = truth >= q2.5 & truth <= q97.5,
         width = q97.5 - q2.5)

by_param <- res |>
  group_by(scenario, model, parameter) |>
  summarise(n_reps = n(), truth = first(truth), mean_est = mean(mean),
            bias = mean(err), rel_bias = mean(err) / abs(first(truth)),
            rmse = sqrt(mean(err^2)), coverage = mean(covered),
            mean_width = mean(width), max_rhat = max(rhat), min_ess = min(ess_bulk),
            .groups = "drop")

sampler <- fits |>
  distinct(scenario, rep, model, n_div, n_treedepth, minutes) |>
  group_by(scenario, model) |>
  summarise(fits = n(), mean_minutes = mean(minutes),
            fits_with_div = sum(n_div > 0), mean_div = mean(n_div),
            fits_with_treedepth = sum(n_treedepth > 0), .groups = "drop")

write.csv(by_param, file.path(out_dir, "summary_by_parameter.csv"), row.names = FALSE)
write.csv(sampler, file.path(out_dir, "sampler_summary.csv"), row.names = FALSE)
write.csv(res, file.path(out_dir, "all_fits.csv"), row.names = FALSE)

# Key comparison table: flow effect and survival to KL
key <- by_param |>
  filter(parameter %in% c("bQ", "S_KL", "delta_AT", "bp")) |>
  select(scenario, model, parameter, truth, mean_est, rel_bias, coverage, mean_width) |>
  arrange(parameter, scenario, model)
write.csv(key, file.path(out_dir, "key_results.csv"), row.names = FALSE)
print(as.data.frame(key |> mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)

# Figures (Okabe-Ito) ----------------------------------------------------------------
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
names(okabe_ito) <- model_labels
shapes <- c(16, 17, 15, 18)
names(shapes) <- model_labels

plot_pars <- c("bQ", "bL", "delta_AT", "S_KL", "TT_KL", "gQ", "bp", "a0")

p_err <- res |>
  filter(parameter %in% plot_pars) |>
  mutate(model = model_labels[model],
         parameter = factor(parameter, levels = plot_pars)) |>
  ggplot(aes(x = err, y = model, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2, alpha = 0.7, position = position_jitter(height = 0.15, seed = 1)) +
  facet_grid(scenario ~ parameter, scales = "free_x") +
  scale_colour_manual(values = okabe_ito) +
  scale_shape_manual(values = shapes) +
  labs(x = "Posterior mean minus truth", y = NULL, colour = NULL, shape = NULL,
       title = "Joint CJS + CWT simulation: estimation error by model",
       subtitle = "One point per replicate. S1: KL efficiency constant in flow; S2: KL efficiency falls with flow") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", axis.text.y = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "error_plot.png"), p_err, width = 14, height = 6, dpi = 200, bg = "white")

p_bq <- res |>
  filter(parameter == "bQ") |>
  mutate(model = model_labels[model]) |>
  ggplot(aes(x = rep, y = mean, ymin = q2.5, ymax = q97.5, colour = model, shape = model)) +
  geom_hline(aes(yintercept = truth), linetype = "dashed", colour = "grey40") +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) +
  facet_wrap(~scenario, ncol = 1) +
  scale_colour_manual(values = okabe_ito) +
  scale_shape_manual(values = shapes) +
  labs(x = "Replicate", y = "Flow effect on survival (bQ), 95% interval",
       colour = NULL, shape = NULL,
       title = "Flow effect on survival: does the joint model recover it?",
       subtitle = "Dashed line = truth (0.5)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "flow_effect_by_replicate.png"), p_bq, width = 10, height = 7, dpi = 200, bg = "white")
