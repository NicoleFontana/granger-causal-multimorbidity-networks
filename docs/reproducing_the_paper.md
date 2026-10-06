# Correspondence between the paper and the code

| Paper | Code | Notes |
|---|---|---|
| Methods – event definitions, recording rules (Supplementary Table S9) | `data/event_definitions.csv`, `R/01_build_sequences.R` | mapping of registry codes to event types is registry-specific and precedes this pipeline |
| Methods – ages 40–80, calendar window, age in months as process clock | `R/01_build_sequences.R`, `config/settings.json` (`study_design`) | the UK Biobank replication used 1995–2019 (`calendar_start`) |
| Methods – primary, extended and sex-stratified cohorts | `config/settings.json` (`analyses`), `R/01_build_sequences.R` | |
| Methods – multivariate Hawkes process, Gaussian basis kernels, MLE-SGL | `matlab/pipeline/fit_network.m` (toolkit: `Initialization_Basis`, `Learning_MLE_Basis`) | L = 10 landmarks (`default_landmarks.m`), w = 12 months, lags up to 120 months |
| Methods – BIC grid search (Supplementary Table S10) | `matlab/pipeline/run_grid_search.m` | 70/30 split of individuals, seed 42 |
| Edge strengths (integrated excitation kernels) | `fit_network.m` (`ImpactFunc`), `results/models/*/excitation_matrix.csv` | rows: cause, columns: effect |
| Weighted in/out-degree, ten strongest edges, excitation kernels (Figs 3–4, Table 1) | `reports/network_report.Rmd` | |
| Risk-factor adjustment: CMD-to-CMD comparison, Wilcoxon signed-rank, 95th-percentile pairs, risk factor → CMD edges, shared risk-factor attribution and rank-sum test (Fig. 5, Table 2, Supplementary Tables S4–S5) | `reports/comparison_report.Rmd` (`comparison_risk_factor_adjustment.html`) | all tests two-sided |
| Extended cohort comparison | `comparison_extended_cohort.html` | |
| Sex-stratified comparison (Fig. 6) | `comparison_sex_stratified.html` | differences are men minus women |
| UK Biobank replication: Pearson correlation, top-5% overlap proportion, Jaccard index, shared/exclusive edges (Table 3) | `comparison_replication.html` | |
| Simulation study (Supplementary Methods S2, Tables S11–S12) | `matlab/run_simulations.m('paper')`, `matlab/simulation/` | N = 500,000: HPC scale. The original runs were not seeded, so reruns give close but not identical numbers |
| Additional simulations (not in the paper): structural recovery across N and D | `matlab/run_simulations.m('nd_sweep')` | 17 scenarios, D = 10–75, N = 1,000–500,000, penalties fixed at 100/100 |
| Additional simulations (not in the paper): sensitivity to the kernel basis | `matlab/run_simulations.m('misspec')` | D = 25, N = 20,000; nine estimation bases (bandwidth, number and position of landmarks) under three ground truths |
