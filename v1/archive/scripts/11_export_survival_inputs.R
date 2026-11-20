# ============================================================
# 11_export_survival_inputs.R
# Handoff files for the integrated model (runSRJPE), from the new stm:
#   1. Thin the stm to 1,500 draws (the integrated model samples 500; 8,000
#      would make the CSV too big for GitHub). Same draws kept in every array.
#   2. Write the survival-input format runSRJPE now expects (Liz's
#      extract_stm_inputs.R, same layout): survival_draws.csv.gz and
#      survival_sites.csv, with the discharge covariate columns filled in.
#   3. Round-trip check: the CSV reproduces the stm slices.
# Outputs: results/handoff/stm.Rdata, survival_draws.csv.gz, survival_sites.csv
# Memory: the full stm is roughly 10 GB in R. Restart R first. A few minutes.
# ============================================================
library(tidyverse); library(here)

stm_file <- here("results", "stm_MaxFlowFL_KLdelta.Rdata")
out_dir  <- here("results", "handoff")
n_keep   <- 1500
signif_digits <- 10
dir.create(out_dir, showWarnings = FALSE)

load(stm_file)                                   # object: stm (8,000 draws)
nsim <- dim(stm$SurvForecastSz_rst)[1]
keep <- unique(round(seq(1, nsim, length.out = n_keep)))

# ---- 1. thin: every draw-indexed object (first dimension = nsim) -------------
thin_one <- function(x) {
  if (is.array(x) && dim(x)[1] == nsim) {
    do.call(`[`, c(list(x, keep), rep(list(TRUE), length(dim(x)) - 1), drop = FALSE))
  } else if (is.null(dim(x)) && is.numeric(x) && length(x) == nsim) {
    x[keep]
  } else x
}
stm_full <- stm
stm <- map(stm_full, thin_one)
thinned <- names(stm)[map_lgl(names(stm), \(nm) !identical(dim(stm[[nm]]), dim(stm_full[[nm]])) ||
                                               length(stm[[nm]]) != length(stm_full[[nm]]))]
message("Thinned ", length(thinned), " objects to ", length(keep), " draws: ", paste(thinned, collapse = ", "))
stopifnot(identical(stm$SurvForecastSz_rst, stm_full$SurvForecastSz_rst[keep, , , , drop = FALSE]),
          identical(stm$TribSurvForecastSz_rst, stm_full$TribSurvForecastSz_rst[keep, , , , drop = FALSE]))
rm(stm_full); invisible(gc())
save(stm, file = file.path(out_dir, "stm.Rdata"))

# ---- 2. survival-input format (as extract_stm_inputs.R) ----------------------
# RST index -> site, unchanged from Flora's data script (02: dist_rstwoodson,
# dist_rstsac, TribForRST = 1, 2, 2). Redding and Red Bluff (mainstem 5, 6) dropped.
rst_map <- tribble(
  ~site,        ~stream,         ~stm_array,               ~stm_rst_index, ~flow_grid, ~flow_covariate_stream, ~flow_covariate_site_group,
  "ucc",        "clear creek",   "SurvForecastSz_rst",     1,              "mainstem", "sacramento river",     "knights landing",
  "ubc",        "battle creek",  "SurvForecastSz_rst",     2,              "mainstem", "sacramento river",     "knights landing",
  "mill creek", "mill creek",    "SurvForecastSz_rst",     3,              "mainstem", "sacramento river",     "knights landing",
  "deer creek", "deer creek",    "SurvForecastSz_rst",     4,              "mainstem", "sacramento river",     "knights landing",
  "okie dam",   "butte creek",   "TribSurvForecastSz_rst", 1,              "butte",    "butte creek",          "butte creek",
  "hallwood",   "yuba river",    "TribSurvForecastSz_rst", 2,              "feather",  "feather river",        "feather river",
  "lfc",        "feather river", "TribSurvForecastSz_rst", 3,              "feather",  "feather river",        "feather river"
)

flow_grids <- list(mainstem = stm$Xvec * stm$sdx + stm$mux,
                   butte    = stm$XvecT[, 1] * stm$sdxB + stm$muxB,
                   feather  = stm$XvecT[, 2] * stm$sdxF + stm$muxF)
fork_length_grid <- stm$Xsz * stm$sdsz + stm$musz

# expand.grid varies its first argument fastest, matching R's array order
site_draws <- function(site, stm_array, stm_rst_index, flow_grid, ...) {
  arr <- stm[[stm_array]][, , , stm_rst_index]
  expand.grid(sim = seq_len(dim(arr)[1]), fork_length_mm = fork_length_grid,
              flow_cfs = flow_grids[[flow_grid]]) |>
    mutate(site = site, endpoint = "delta entry", survival = as.vector(arr))
}
survival_draws <- pmap(rst_map, site_draws) |>
  bind_rows() |>
  mutate(across(c(flow_cfs, fork_length_mm, survival), \(x) signif(x, signif_digits))) |>
  select(site, endpoint, flow_cfs, fork_length_mm, sim, survival)

survival_sites <- rst_map |>
  transmute(site, stream,
            flow_covariate_name = "monthly_max_flow",          # monthly maximum in the release month
            flow_covariate_stream, flow_covariate_site_group,
            source = paste0("stm_MaxFlowFL_KLdelta.Rdata (05 on fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL, 03g; ",
                            "SRJPEdata d3584db), thinned to ", length(keep), " draws by 11_export_survival_inputs.R"),
            source_md5 = unname(tools::md5sum(stm_file)),
            created_on = as.character(Sys.Date()))

draws_path <- file.path(out_dir, "survival_draws.csv.gz")
write_csv(survival_draws, draws_path)
write_csv(survival_sites, file.path(out_dir, "survival_sites.csv"))

# ---- 3. round-trip check ------------------------------------------------------
draws_back <- read_csv(draws_path, col_types = "ccdddd")
round_trip_error <- pmap_dbl(rst_map, function(site, stm_array, stm_rst_index, ...) {
  original <- stm[[stm_array]][, , , stm_rst_index]
  read_back <- draws_back |> filter(site == !!site) |> arrange(flow_cfs, fork_length_mm, sim) |>
    pull(survival) |> array(dim = dim(original))
  max(abs(read_back - original) / abs(original))
})
stopifnot(round_trip_error < 10^-(signif_digits - 1))

message("Wrote ", nrow(survival_draws), " rows to ", draws_path, " (",
        round(file.size(draws_path) / 1e6, 1), " MB; round trip OK, max relative error ",
        signif(max(round_trip_error), 2), ")")
print(survival_draws |> summarise(flow_min = min(flow_cfs), flow_max = max(flow_cfs),
                                  surv_median = median(survival), .by = site))
