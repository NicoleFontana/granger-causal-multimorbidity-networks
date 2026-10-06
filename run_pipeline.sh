#!/usr/bin/env bash
# Run the full pipeline on the synthetic data (about 5 minutes):
#   1. sequences and cohorts (R)
#   2. grid search and network estimation (MATLAB)
#   3. analysis reports (R Markdown)
# for the discovery cohort and for a replication cohort generated from the
# same synthetic network.
#
# Usage (from the repository root):  bash run_pipeline.sh
# Set MATLAB to the MATLAB executable if `matlab` is not on the PATH, e.g.
#   MATLAB=/Applications/MATLAB_R2023a.app/bin/matlab bash run_pipeline.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
MATLAB="${MATLAB:-matlab}"

[ -d external/Hawkes-Process-Toolkit ] || bash setup/install_toolkit.sh

# synthetic data (included in the repository; regenerated only if missing)
[ -f data/synthetic/events.csv ] || Rscript synthetic/generate_synthetic_data.R 5000 2026 data/synthetic 2026
[ -f data/synthetic_replication/events.csv ] || Rscript synthetic/generate_synthetic_data.R 3000 99 data/synthetic_replication 2026

echo "== 1. Sequences and cohorts"
Rscript R/01_build_sequences.R data/synthetic results/sequences
Rscript R/01_build_sequences.R data/synthetic_replication results_replication/sequences baseline

echo "== 2. Estimation (MATLAB)"
(cd matlab && "${MATLAB}" -batch "run_all_estimations; run_estimation('baseline', fullfile('..','results_replication','sequences'), fullfile('..','results_replication','models'));")

echo "== 3. Reports"
Rscript R/03_render_reports.R results results_replication
echo "Reports written to results/reports/"
