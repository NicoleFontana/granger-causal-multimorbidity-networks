#!/usr/bin/env Rscript
# =============================================================================
# Step 3 - Analysis reports (HTML) for each network and for the comparisons in
# the paper. Run after the MATLAB estimation (matlab/run_all_estimations.m).
#
# Usage (from the repository root):
#   Rscript R/03_render_reports.R [results_dir] [replication_results_dir]
# =============================================================================

suppressPackageStartupMessages(library(rmarkdown))

args    <- commandArgs(trailingOnly = TRUE)
res     <- normalizePath(if (length(args) >= 1) args[1] else "results")
res_rep <- if (length(args) >= 2) normalizePath(args[2]) else NA
out     <- file.path(res, "reports")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

analyses <- list.dirs(file.path(res, "models"), full.names = FALSE, recursive = FALSE)
for (a in analyses) {
  render("reports/network_report.Rmd", output_dir = out, output_file = sprintf("network_%s.html", a),
         params = list(analysis = a, results_dir = res), quiet = TRUE, envir = new.env())
  cat("Rendered network report:", a, "\n")
}

comparisons <- list(
  risk_factor_adjustment = list(analysis_a = "baseline", analysis_b = "risk_adjusted",
                                label_a = "baseline", label_b = "risk-factor-adjusted", attribution = TRUE),
  extended_cohort        = list(analysis_a = "risk_adjusted", analysis_b = "extended",
                                label_a = "primary cohort", label_b = "extended cohort", attribution = FALSE),
  sex_stratified         = list(analysis_a = "baseline_women", analysis_b = "baseline_men",
                                label_a = "women", label_b = "men", attribution = FALSE))
if (!is.na(res_rep))
  comparisons$replication <- list(analysis_a = "baseline", analysis_b = "baseline",
                                  label_a = "discovery cohort", label_b = "replication cohort",
                                  results_dir_b = res_rep, attribution = FALSE, overlap = TRUE)

for (n in names(comparisons)) {
  p <- modifyList(list(results_dir_a = res, results_dir_b = res, overlap = FALSE), comparisons[[n]])
  if (!all(dir.exists(file.path(c(p$results_dir_a, p$results_dir_b), "models", c(p$analysis_a, p$analysis_b))))) {
    cat("Skipped comparison (models missing):", n, "\n"); next
  }
  render("reports/comparison_report.Rmd", output_dir = out, output_file = sprintf("comparison_%s.html", n),
         params = p, quiet = TRUE, envir = new.env())
  cat("Rendered comparison report:", n, "\n")
}
