# 01b_filter_cwt_reading_years.R
#
# Drop CWT release groups released when Knights Landing was not reading CWTs.
#
# Why: the Knights Landing CWT recovery data (SRJPEdata::rst_cwt_recaptures)
# start in December 1999 and have essentially no reads in 2013 or 2022, even
# though the screw trap was fishing (SRJPEdata::weekly_hours_fished). Groups
# released then have zero recoveries because nobody was reading tags, not
# because no fish survived. The binomial CWT likelihood reads those zeros as
# "none survived to be caught", which inflates sd_cwtSre, drags pCap and
# baseline survival down, and (because 2013 and 2022 were dry years) can
# create a spurious positive flow effect that is soft-shared into CJS.
#
# Rule: keep a release group only if at least one CWT (any tag code, any
# release group) was read at Knights Landing in the window_days after its
# release. A group with its own recoveries is always kept, since reads prove
# tags were being processed. The rule never uses a zero-recovery group's own
# outcome to drop it.
#
# Inputs:  data/hatchery_relrec_all.csv (unfiltered, written by 01_PrepData.R)
#          or, the first time, data/hatchery_relrec.csv (copied to _all)
# Outputs: data/hatchery_relrec.csv (filtered; read by 02, 06, 07, 09)
#          results/tables/cwt_groups_kl_reading_filter.csv (every group, kept or not)
#          results/tables/cwt_groups_kl_reading_filter_by_year.csv
#
# Can be run on its own (no ERDDAP pull) or sourced at the end of 01_PrepData.R.

library(tidyverse)
library(here)

window_days <- 60 # matches the cwt_tt < 60 travel-time filter in 01_PrepData.R

all_file <- here("data", "hatchery_relrec_all.csv")
cur_file <- here("data", "hatchery_relrec.csv")

# Keep one unfiltered copy. 01_PrepData.R writes _all directly; the first
# standalone run creates it from the existing (unfiltered) hatchery_relrec.csv.
if (!file.exists(all_file)) {
  stopifnot(file.copy(cur_file, all_file))
  message("Saved the unfiltered CWT table as ", all_file)
}

drerelrec_all <- read.csv(all_file, stringsAsFactors = FALSE)

# Every CWT read at Knights Landing, any tag code, any release group
kl_read_days <- SRJPEdata::rst_cwt_recaptures |>
  filter(!is.na(date), !is.na(tag_code)) |>
  pull(date) |>
  as.Date() |>
  as.numeric() |>
  sort()

reads_between <- function(from, to) {
  findInterval(to, kl_read_days) - findInterval(from - 1, kl_read_days)
}

own_recoveries <- drerelrec_all |>
  filter(!is.na(tag_code)) |>
  distinct(release_group_id, tag_code, date, forklength) |>
  count(release_group_id, name = "own_recoveries")

group_filter <- drerelrec_all |>
  distinct(release_group_id, release_location_name, year, mid_release_date) |>
  mutate(release_day = as.numeric(as.Date(mid_release_date))) |>
  left_join(own_recoveries, by = "release_group_id") |>
  mutate(
    own_recoveries = replace_na(own_recoveries, 0L),
    kl_reads_in_window = reads_between(release_day, release_day + window_days),
    kl_reading = kl_reads_in_window > 0
  ) |>
  select(-release_day) |>
  arrange(year, mid_release_date)

stopifnot(
  !anyDuplicated(group_filter$release_group_id),          # one row per group
  all(group_filter$kl_reading[group_filter$own_recoveries > 0])  # own reads imply reading
)

by_year <- group_filter |>
  group_by(year) |>
  summarise(
    groups = n(),
    kept = sum(kl_reading),
    dropped = sum(!kl_reading),
    zero_recovery_groups = sum(own_recoveries == 0),
    zero_recovery_kept = sum(own_recoveries == 0 & kl_reading),
    .groups = "drop"
  )

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(group_filter, here("results", "tables", "cwt_groups_kl_reading_filter.csv"), row.names = FALSE)
write.csv(by_year, here("results", "tables", "cwt_groups_kl_reading_filter_by_year.csv"), row.names = FALSE)

drerelrec <- drerelrec_all |>
  semi_join(filter(group_filter, kl_reading), by = "release_group_id")

write.csv(drerelrec, cur_file, row.names = FALSE)

message(
  "CWT groups kept: ", sum(group_filter$kl_reading), " of ", nrow(group_filter),
  "; dropped (no Knights Landing CWT reads within ", window_days, " days of release): ",
  sum(!group_filter$kl_reading), "\nDropped release years: ",
  paste(by_year$year[by_year$dropped > 0], collapse = ", ")
)
print(by_year, n = Inf)
