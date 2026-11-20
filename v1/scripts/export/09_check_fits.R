# ============================================================
# scripts/export/09_check_fits.R      (new, 18 Nov 2026)
#
# Run first in the export tier. Checks that the fit files in results/ are the
# ones the current handoff was made from, so a wrong or stale copy (for example
# downloaded from SharePoint) stops here instead of silently giving different
# numbers downstream.
#
# Reference checksums: data/fit_checksums.csv (in git).
#   - If it does not exist, this script writes it from the fits present now
#     (do this once, on the machine with the production fits).
#   - After a deliberate refit (03 or 07), delete it and run this again to
#     record the new fits; commit the new file with the new handoff.
#
# About a minute (md5 of the 6.5 GB telemetry fit).
# ============================================================

library(tidyverse)
library(here)

ref_file <- here("data", "fit_checksums.csv")
fits <- tribble(
  ~fit,           ~file,                                                    ~used_by,
  "telemetry",    "fit_CovIndCont_studymonth_ncRE_KLdelta_MaxFlow_FL.Rdata", "10, 11",
  "cwt_gamma_G1", "fit_CWT_gamma_G1.Rdata",                                 "13"
) |>
  mutate(path = here("results", file), present = file.exists(path))

if (!all(fits$present)) {
  print(select(filter(fits, !present), fit, file, used_by))
  stop("Fit files missing from results/ (above). Copy them from the SharePoint folder.")
}

cat("Checksumming fits...\n")
now <- fits |> mutate(bytes = file.size(path), md5 = unname(tools::md5sum(path)))

if (!file.exists(ref_file)) {
  write_csv(select(now, fit, file, bytes, md5), ref_file)
  cat("No reference yet: wrote", ref_file, "from the fits present now. Commit it.\n")
} else {
  ref <- read_csv(ref_file, show_col_types = FALSE)
  cmp <- now |>
    left_join(select(ref, fit, md5_ref = md5), by = "fit") |>
    mutate(match = !is.na(md5_ref) & md5 == md5_ref)
  print(select(cmp, fit, file, match) |> as.data.frame(), right = FALSE, row.names = FALSE)
  if (!all(cmp$match))
    stop("Fit file(s) differ from data/fit_checksums.csv. Wrong download, or a refit: ",
         "if you refit on purpose, delete data/fit_checksums.csv and run this again.")
  cat("Fits match data/fit_checksums.csv.\n")
}
