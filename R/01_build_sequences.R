#!/usr/bin/env Rscript
# =============================================================================
# Step 1 - From event tables to event sequences, for each analysis cohort.
#
# Input  (see docs/data_format.md): <data_dir>/persons.csv, <data_dir>/events.csv
#        and data/event_definitions.csv
# Output (one folder per analysis in <out_dir>/<analysis>/):
#   sequences.csv       patient_id, Time, Mark   (input of the MATLAB estimation)
#   event_mapping.csv   Mark, event_id, event_name, group
#   cohort_summary.csv  descriptive statistics of the cohort
#
# Study design (config/settings.json):
#   - events recorded within the calendar window and between ages 40 and 80;
#   - process clock: age in months since age 40;
#   - recording rule 1 (chronic conditions): first recorded diagnosis only;
#     rules 2 and 3: all recorded occurrences;
#   - cohort: individuals with at least `min_events` events of the inclusion
#     event groups (CMD for the primary cohort, CMD or RISK for the extended
#     cohort); the sex-stratified cohorts are the primary cohort split by sex;
#   - event types renumbered 1..D (CMD first, then RISK).
#
# Usage (from the repository root):
#   Rscript R/01_build_sequences.R [data_dir] [out_dir] [analysis ...]
#   e.g. Rscript R/01_build_sequences.R data/synthetic results/sequences
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path("data", "synthetic")
out_dir  <- if (length(args) >= 2) args[2] else file.path("results", "sequences")

cfg      <- fromJSON(file.path("config", "settings.json"))
design   <- cfg$study_design
analyses <- if (length(args) >= 3) args[-(1:2)] else names(cfg$analyses)

defs    <- fread(file.path("data", "event_definitions.csv"))
persons <- fread(file.path(data_dir, "persons.csv"), colClasses = list(character = "person_id"))
events  <- fread(file.path(data_dir, "events.csv"),  colClasses = list(character = "person_id"))
persons[, birth_date := as.IDate(birth_date)]
events[,  event_date := as.IDate(event_date)]

# ---- Events common to all analyses -----------------------------------------
ev <- persons[events, on = "person_id", nomatch = 0]
ev <- defs[, .(event_id, group, recording_rule)][ev, on = "event_id", nomatch = 0]

# calendar window
ev <- ev[event_date >= as.IDate(design$calendar_start) & event_date <= as.IDate(design$calendar_end)]

# age in months; keep ages min_age..max_age; clock origin at min_age
ev[, Time := as.numeric(event_date - birth_date) / design$days_per_month]
ev <- ev[Time >= design$min_age_years * 12 & Time <= design$max_age_years * 12]
ev[, Time := Time - design$min_age_years * 12]

# identical records (same person, event type and date) counted once
ev <- unique(ev, by = c("person_id", "event_id", "Time"))

# recording rule 1: first recorded diagnosis only
setorder(ev, person_id, Time, event_id)
ev <- ev[recording_rule != 1 | !duplicated(ev[, .(person_id, event_id)])]

# ---- Build one analysis ----------------------------------------------------
build_analysis <- function(name, spec) {
  sel <- ev
  if (spec$sex != "all") sel <- sel[sex == spec$sex]

  # cohort inclusion: >= min_events events of the inclusion groups
  n_incl   <- sel[group %in% spec$inclusion_events, .N, by = person_id]
  included <- n_incl[N >= design$min_events, person_id]
  sel <- sel[person_id %in% included & group %in% spec$events]

  # renumber event types 1..D (CMD first, then RISK, in event_id order)
  map <- defs[group %in% spec$events][order(match(group, c("CMD", "RISK")), event_id)]
  map[, Mark := seq_len(.N)]
  sel <- map[, .(event_id, Mark)][sel, on = "event_id"]

  setorder(sel, person_id, Time, Mark)
  sel[, patient_id := match(person_id, unique(person_id))]

  # descriptive statistics (computed on the analysis events)
  per_person <- sel[, .(n_events = .N, age_first = (min(Time) + design$min_age_years * 12) / 12,
                        sex = sex[1]), by = patient_id]
  summary <- data.table(
    analysis             = name,
    n_individuals        = nrow(per_person),
    n_women              = sum(per_person$sex == "F"),
    pct_women            = round(100 * mean(per_person$sex == "F"), 2),
    n_events             = nrow(sel),
    events_mean          = round(mean(per_person$n_events), 2),
    events_median        = median(per_person$n_events),
    events_min           = min(per_person$n_events),
    events_max           = max(per_person$n_events),
    age_first_event_mean = round(mean(per_person$age_first), 1),
    age_first_event_sd   = round(sd(per_person$age_first), 1),
    n_event_types        = nrow(map),
    n_event_types_observed = uniqueN(sel$Mark))

  dir.create(file.path(out_dir, name), recursive = TRUE, showWarnings = FALSE)
  fwrite(sel[, .(patient_id, Time, Mark)], file.path(out_dir, name, "sequences.csv"))
  fwrite(map[, .(Mark, event_id, event_name, group)], file.path(out_dir, name, "event_mapping.csv"))
  fwrite(summary, file.path(out_dir, name, "cohort_summary.csv"))

  cat(sprintf("%-15s N = %6d | events = %7d | D = %d | women = %.1f%%\n",
              name, summary$n_individuals, summary$n_events, summary$n_event_types, summary$pct_women))
  if (summary$n_event_types_observed < summary$n_event_types)
    warning(sprintf("%s: %d of %d event types never observed", name,
                    summary$n_event_types - summary$n_event_types_observed, summary$n_event_types))
}

for (a in analyses) build_analysis(a, cfg$analyses[[a]])
