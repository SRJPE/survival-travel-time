# ============================================================
# scripts/refit/00_save_taggedfish.R      (new, 18 Nov 2026)
#
# Saves the ERDDAP TaggedFish table once, so 01 runs offline and gives the
# same data every time. In the old project helper_fxt.R pulled it live from
# ERDDAP whenever it was sourced.
#
# TaggedFish gives fish length, weight, type, release date and location for
# every acoustic-tagged fish. The detections themselves are already archived
# in data/*AllYears/ (pulled 8 Oct 2026).
#
# The pull happens today, not on 8 Oct, so the script checks it against the
# data behind the 03g fit (data/reference/*_20261008.csv): every fitted fish
# must be present with the same fork length and release date.
#
# Output: data/TaggedFish.csv.gz (all columns as character, as ERDDAP returns them)
# Needs internet. Run once; it stops if the file already exists.
# ============================================================

library(tidyverse)
library(here)

out <- here("data", "TaggedFish.csv.gz")
if (file.exists(out)) stop(out, " already exists. Delete it to pull again.")

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
tf <- rerddap::tabledap(rerddap::info("FED_JSATS_taggedfish", url = my_url), url = my_url) |>
  as_tibble() |>
  mutate(across(everything(), as.character))
cat("TaggedFish rows:", nrow(tf), "| studies:", n_distinct(tf$study_id),
    "| duplicated fish_id:", sum(duplicated(tf$fish_id)), "\n")

# ---- Check against the fitted fish (8 Oct 2026) ------------------------------------------
ref <- bind_rows(
  read_csv(here("data", "reference", "Sac_data_20261008.csv"), show_col_types = FALSE,
           col_types = cols(.default = col_character())),
  read_csv(here("data", "reference", "FeaBut_data_20261008.csv"), show_col_types = FALSE,
           col_types = cols(.default = col_character()))) |>
  transmute(FishID, fl_ref = as.numeric(fish_length), rel_ref = fish_release_date)

chk <- ref |>
  left_join(tf |> distinct(fish_id, .keep_all = TRUE) |>
              transmute(FishID = fish_id, fl_new = as.numeric(fish_length), rel_new = fish_release_date),
            by = "FishID") |>
  mutate(missing = is.na(rel_new),
         fl_diff = !missing & !(fl_ref == fl_new | (is.na(fl_ref) & is.na(fl_new))),
         rel_diff = !missing & rel_ref != rel_new)

cat("Fitted fish:", nrow(chk), "| missing from today's pull:", sum(chk$missing),
    "| fork length differs:", sum(chk$fl_diff, na.rm = TRUE),
    "| release date differs:", sum(chk$rel_diff, na.rm = TRUE), "\n")
if (any(chk$missing | chk$fl_diff | chk$rel_diff, na.rm = TRUE)) {
  print(filter(chk, missing | fl_diff | rel_diff), n = 20)
  warning("Today's TaggedFish differs from the data behind the 03g fit for the fish above. ",
          "It is saved anyway; 01's reference check will show whether the model data change.")
}

write_csv(tf, out)
cat("Wrote", out, "(", round(file.size(out) / 1e6, 1), "MB )\n")
