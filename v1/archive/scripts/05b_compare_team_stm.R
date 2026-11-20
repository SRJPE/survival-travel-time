# ============================================================
# 05b_compare_team_stm.R
#
# What changes for the integrated model team: their current stm
# (built from fit_CovIndContv2_TT_clean, misaligned travel times,
# older data) vs ours (stm_MaxFlowFL_aligned, from the corrected
# long fit). Compares the forecasts they use, at the same flows in
# cfs, with random release-group effects off (_nore) for survival.
#
# Writes results/tables/compare_team_stm.csv and
# results/figures/compare_team_stm.png.
# ============================================================

library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

stm_files <- c(team = here("data", "stm.Rdata"),
               ours = here("results", "stm_MaxFlowFL_KLdelta.Rdata"))

# median and 95% interval of each grid column, with the grid in cfs
summ_cols <- function(draws, cfs) {
  tibble(cfs = cfs,
         median = apply(draws, 2, median),
         q2.5 = apply(draws, 2, quantile, 0.025),
         q97.5 = apply(draws, 2, quantile, 0.975))
}

extract_forecasts <- function(file, label) {
  e <- new.env(); load(file, envir = e); s <- e$stm
  # flow grids back to cfs using each stm's own scaling (what each pipeline assumes)
  cfs_sac <- s$Xvec * s$sdx + s$mux
  cfs_but <- s$XvecT[[1]] * s$sdxB + s$muxB
  cfs_fea <- s$XvecT[[2]] * s$sdxF + s$muxF
  out <- list(
    summ_cols(s$SurvForecast_nore, cfs_sac)          |> mutate(quantity = "Survival", river = "Sacramento"),
    summ_cols(s$TribSurvForecast_nore[, , 1], cfs_but) |> mutate(quantity = "Survival", river = "Butte"),
    summ_cols(s$TribSurvForecast_nore[, , 2], cfs_fea) |> mutate(quantity = "Survival", river = "Feather"),
    summ_cols(s$TTForecast, cfs_sac)                 |> mutate(quantity = "Travel time (d)", river = "Sacramento"),
    summ_cols(s$TTForecastT[, , 1], cfs_but)         |> mutate(quantity = "Travel time (d)", river = "Butte"),
    summ_cols(s$TTForecastT[, , 2], cfs_fea)         |> mutate(quantity = "Travel time (d)", river = "Feather")
  ) |> list_rbind() |> mutate(stm = label, draws = nrow(s$SurvForecast_nore))
  rm(e); gc()
  out
}

cmp <- imap(stm_files, extract_forecasts) |> list_rbind()

# side-by-side table at a few flows per river (every 6th grid point)
tab <- cmp |>
  group_by(stm, quantity, river) |>
  filter(row_number() %in% c(1, 7, 13, 19, 25)) |>
  ungroup() |>
  mutate(across(c(cfs, median, q2.5, q97.5), \(x) round(x, 2))) |>
  select(quantity, river, stm, cfs, median, q2.5, q97.5) |>
  arrange(quantity, river, cfs, stm)
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(tab, here("results", "tables", "compare_team_stm.csv"), row.names = FALSE)
print(as.data.frame(tab))

# figure: median and 95% interval vs flow; Okabe-Ito colours plus line types
p <- ggplot(cmp, aes(cfs, median, colour = stm, fill = stm, linetype = stm)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(quantity ~ river, scales = "free") +
  scale_colour_manual(values = c(team = "#D55E00", ours = "#0072B2"),
                      labels = c(team = "Team's current stm", ours = "Corrected (ours)")) +
  scale_fill_manual(values = c(team = "#D55E00", ours = "#0072B2"),
                    labels = c(team = "Team's current stm", ours = "Corrected (ours)")) +
  scale_linetype_manual(values = c(team = "dashed", ours = "solid"),
                        labels = c(team = "Team's current stm", ours = "Corrected (ours)")) +
  labs(x = "Flow (cfs)", y = NULL, colour = NULL, fill = NULL, linetype = NULL,
       title = "Forecasts the team uses: current stm vs corrected",
       subtitle = "Survival without release-group effect (_nore); travel time with a new release group drawn") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("results", "figures", "compare_team_stm.png"), p,
       width = 11, height = 6.5, dpi = 200, bg = "white")
p
