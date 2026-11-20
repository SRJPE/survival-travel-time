yr_tbl <- cwt_groups |>
  mutate(p_bt = p_bt$p_mean, flow_g = as.numeric(Maxflowsac), w = as.numeric(Nrel)) |>
  group_by(year) |>
  summarise(flow = weighted.mean(flow_g, w),
            btspas_capture = weighted.mean(p_bt, w)) |>
  inner_join(oe |> filter(c == "c1", exp_median >= 5) |> select(year, OE), by = "year")
cor.test(yr_tbl$OE, yr_tbl$flow, method = "spearman")
cor.test(yr_tbl$btspas_capture, yr_tbl$flow, method = "spearman")
yr_tbl |> arrange(flow) |> mutate(across(where(is.double), \(x) signif(x, 3))) |> as.data.frame()

p3 <- yr_tbl |>
  pivot_longer(c(OE, btspas_capture), names_to = "what", values_to = "value") |>
  mutate(what = recode(what, OE = "Observed / expected recoveries (Lorenzen c = 1)",
                       btspas_capture = "BT-SPAS-X capture probability (arrival-weighted)")) |>
  ggplot(aes(flow, value, label = year)) +
  geom_smooth(method = "lm", se = FALSE, colour = "#D55E00", linetype = "dashed", linewidth = 0.7) +
  geom_point(colour = "#0072B2", size = 2) +
  geom_text(size = 3, vjust = -0.8, colour = "grey30") +
  scale_y_log10() +
  facet_wrap(~ what, scales = "free_y", ncol = 1) +
  labs(x = "Knights Landing monthly max flow in release month (cfs, release-weighted)", y = NULL,
       title = "Does the year-to-year scatter track flow?") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(here("results", "figures", "lorenzen_OE_vs_flow.png"), p3,
       width = 8, height = 7, dpi = 200, bg = "white")
p3