# ============================================================
# scripts/export/19_write_handoff.R      (new, 18 Nov 2026; dated folders 19 Nov 2026)
#
# Writes one dated handoff, handoff/<handoff_date>/, with two folders and a
# MANIFEST.csv in each:
#
#   runSRJPE/     small files a JPE run needs: parameter draws, site tables,
#                 prediction functions, discharge scenarios and the fit
#                 checksums. Goes into the runSRJPE repo (git).
#   sharepoint/   the fits (needed only to re-export without refitting) plus an
#                 archive copy of everything in runSRJPE/. Goes to SharePoint,
#                 in a folder named for the same date. Not in git (.gitignore).
#
# A handoff that has been sent is never changed. This script stops if
# handoff/<handoff_date>/ already exists; a new handoff gets a new date (or a
# suffix, e.g. "2026-11-19b"). Code changes after a handoff reach people through
# git, not by rewriting a sent folder.
#
# The stm is no longer handed off: travel time comes from
# travel_time_pars.csv.gz and predict_travel_time().
#
# Run after 09, 12, 12b, 13 and 15. About two minutes (copies and checksums
# the 6.5 GB telemetry fit). Set copy_fits <- FALSE to skip copying the fits
# (they are still listed, with checksums, in the manifest).
# ============================================================

library(tidyverse)
library(here)

copy_fits    <- TRUE
handoff_date <- format(Sys.Date(), "%Y-%m-%d")   # or a fixed label, e.g. "2026-11-19b"

out_root <- here("handoff", handoff_date)
if (dir.exists(out_root))
  stop("handoff/", handoff_date, " already exists. A sent handoff is never overwritten: ",
       "set handoff_date to a new label (e.g. \"", handoff_date, "b\").")

rj <- \(...) file.path(out_root, "runSRJPE", ...)
sp <- \(...) file.path(out_root, "sharepoint", ...)
ex <- \(...) here("results", "survival_export", ...)

small <- tribble(
  ~from,                                                ~file,                                   ~made_by,                         ~what,
  ex("survival_pars.csv.gz"),                           "survival_pars.csv.gz",                  "12_export_survival_pars.R",      "survival parameter draws (8,000), one row per draw",
  ex("survival_sites.csv"),                             "survival_sites.csv",                    "12_export_survival_pars.R",      "RST sites: curve, distances, gauge, scaling, fitted discharge range",
  ex("travel_time_pars.csv.gz"),                        "travel_time_pars.csv.gz",               "12b_export_travel_time_pars.R",  "travel-time parameter draws, same draw index as survival_pars",
  ex("travel_time_sites.csv"),                          "travel_time_sites.csv",                 "12b_export_travel_time_pars.R",  "survival_sites plus reach lengths km1-km3",
  ex("lorenzen_pars.csv.gz"),                           "lorenzen_pars.csv.gz",                  "13_export_lorenzen_pars.R",      "Lorenzen exponent c per draw (only for small_fish = 'lorenzen')",
  here("R", "survival_predict_functions.R"),            "survival_predict_functions.R",          "R/",                             "predict_survival(), predict_travel_time(), draw_re()",
  here("results", "tables", "12_scenario_seasons.csv"), "discharge/scenario_seasons.csv",        "15_discharge_scenarios.R",       "seasons defining dry / average / wet",
  here("results", "tables", "12_monthly_max_by_season.csv"), "discharge/monthly_max_by_season.csv", "15_discharge_scenarios.R",   "monthly max discharge by gauge, season, month",
  here("data", "fit_checksums.csv"),                    "fit_checksums.csv",                     "09_check_fits.R",                "md5 of the fits these files came from",
  here("handoff", "README.md"),                         "README.md",                             "handoff/README.md",              "what the files are and how to use them"
)

fits <- tribble(
  ~from,                                                                     ~file,                                                   ~made_by,                ~what,
  here("results", "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata"), "fits/fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata", "03_fit_telemetry.R", "telemetry fit; input to export steps 10 and 11",
  here("results", "fit_CWT_gamma_G1.Rdata"),                                 "fits/fit_CWT_gamma_G1.Rdata",                           "07_fit_cwt_gamma.R",    "CWT gamma fit G1; input to export step 13"
)

missing <- bind_rows(small, fits) |> filter(!file.exists(from))
if (nrow(missing) > 0) { print(select(missing, from, made_by)); stop("Run the steps above first") }

# Fits must still match the reference checksums (09)
ref <- read_csv(here("data", "fit_checksums.csv"), show_col_types = FALSE)
fit_md5 <- unname(tools::md5sum(fits$from))
stopifnot("fits differ from data/fit_checksums.csv: run 09" = setequal(fit_md5, ref$md5))

copy_set <- function(tab, dest) {
  to <- dest(tab$file)
  walk(unique(dirname(to)), \(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))
  stopifnot(all(file.copy(tab$from, to, overwrite = TRUE, copy.date = TRUE)))
  to
}
rj_to <- copy_set(small, rj)
sp_to <- copy_set(small, \(f) sp("runSRJPE_copy", f))
if (copy_fits) {
  cat("Copying fits to handoff/sharepoint/fits (about 6.5 GB)...\n")
  copy_set(fits, sp)
}

# ---- Provenance ----------------------------------------------------------------------------
git_sha <- tryCatch(system2("git", c("-C", shQuote(here()), "rev-parse", "--short", "HEAD"), stdout = TRUE),
                    error = \(e) NA_character_, warning = \(w) NA_character_)
cmdstan_v <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = \(e) NA_character_)
written <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M"), " (handoff ", handoff_date, ")")

provenance <- tibble(
  file = c("SRJPEdata commit", "R version", "CmdStan version", "survival-travel-time git commit"),
  made_by = c(tryCatch(packageDescription("SRJPEdata")$RemoteSha %||% "unknown", error = \(e) "not installed"),
              R.version.string, cmdstan_v, git_sha[1] %||% NA_character_),
  what = c("data version used by 01, 02 and 15", "", "", "uncommitted changes are not recorded"))

manifest_small <- small |>
  mutate(bytes = file.size(rj_to), md5 = unname(tools::md5sum(rj_to))) |>
  select(file, bytes, md5, made_by, what)
manifest_fits <- fits |>
  mutate(bytes = file.size(from), md5 = fit_md5) |>
  select(file, bytes, md5, made_by, what)

manifest_rj <- bind_rows(manifest_small, mutate(manifest_fits, what = paste(what, "(in SharePoint)")),
                         provenance) |>
  mutate(written = written)
manifest_sp <- bind_rows(mutate(manifest_small, file = paste0("runSRJPE_copy/", file)),
                         manifest_fits, provenance) |>
  mutate(written = written)

write_csv(manifest_rj, rj("MANIFEST.csv"))
write_csv(manifest_sp, sp("MANIFEST.csv"))
print(manifest_rj |> select(file, bytes, made_by), n = Inf, width = Inf)
cat("\nWrote", rj(), "(to runSRJPE)\nand", sp(),
    if (copy_fits) "(to SharePoint)\n" else "(to SharePoint; fits NOT copied, copy_fits = FALSE)\n")
