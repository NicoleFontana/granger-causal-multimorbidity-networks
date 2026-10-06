# Input data format

The pipeline starts from three tables. Registry data (Danish National Patient
Register, UK Biobank) cannot be shared; `synthetic/generate_synthetic_data.R`
creates fake tables with exactly this structure so that the code can be run
end to end.

## `data/event_definitions.csv`

One row per event type (29 cardiometabolic diseases, CMD, and 38 risk-factor
conditions, RISK), as in Supplementary Table S9 of the paper.

| column | type | description |
|---|---|---|
| `event_id` | integer | 1–29 for CMD, 30–67 for RISK |
| `event_name` | text | event label |
| `group` | `CMD` / `RISK` | event group |
| `icd10` | text | ICD-10 codes defining the event |
| `procedure_codes` | text | OPCS-4 (UK Biobank) and NOMESCO (Danish registers) procedure codes, if any |
| `recording_rule` | 1 / 2 / 3 | 1 = chronic, first recorded diagnosis only; 2 = chronic with re-exacerbation, each hospitalisation with the condition as main diagnosis; 3 = episodic, all occurrences |

Mapping registry codes to `event_id` is registry-specific and done before this
pipeline: `events.csv` already contains event types.

## `persons.csv`

| column | type | description |
|---|---|---|
| `person_id` | text | pseudonymous identifier |
| `sex` | `F` / `M` | used for the sex-stratified analysis |
| `birth_date` | date (YYYY-MM-DD) | used to compute age at each event |

## `events.csv`

One row per recorded event.

| column | type | description |
|---|---|---|
| `person_id` | text | links to `persons.csv` |
| `event_id` | integer | links to `event_definitions.csv` |
| `event_date` | date (YYYY-MM-DD) | date of the hospital contact |

## From tables to event sequences

The sequence-building step (see the main README) applies the study design:
calendar window, ages 40–80, age in months as the process clock (origin at age
40), recording rules, cohort inclusion criteria, and renumbering of event types.
Its output, used by the MATLAB estimation, is a long table `patient_id, Time,
Mark` (Time in months since age 40; Mark = event type 1..D), converted to the
toolkit's `Seqs` structure (`Time`, `Mark`, `Start`, `Stop`, `Feature`).
