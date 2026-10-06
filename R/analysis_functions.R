# =============================================================================
# Functions shared by the analysis reports (reports/*.Rmd).
# All excitation matrices are oriented A[cause, effect]; edge strength is the
# integrated excitation kernel.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(patchwork)
  library(R.matlab)
})

EDGE_THRESHOLD <- 1e-3   # strength above which a directed edge is considered present

# ---- Loading ----------------------------------------------------------------

#' Read a fitted model written by matlab/pipeline/run_estimation.m
read_model <- function(model_dir, seq_dir = NULL) {
  A <- as.matrix(read.csv(file.path(model_dir, "excitation_matrix.csv"), row.names = 1, check.names = FALSE))
  mu <- fread(file.path(model_dir, "baseline_intensity.csv"))
  names <- mu$event_name
  dimnames(A) <- list(names, names)
  m <- readMat(file.path(model_dir, "model.mat"))
  Phi <- m$Phi                                   # Phi[cause, lag, effect]
  out <- list(A = A, mu = setNames(mu$mu, names), Phi = Phi, names = names)
  if (!is.null(seq_dir)) {
    out$mapping <- fread(file.path(seq_dir, "event_mapping.csv"))
    out$cohort  <- fread(file.path(seq_dir, "cohort_summary.csv"))
  }
  out
}

#' CMD-to-CMD submatrix of a model with CMD and RISK events
cmd_submatrix <- function(model) {
  cmd <- model$mapping$event_name[model$mapping$group == "CMD"]
  model$A[cmd, cmd]
}

# ---- Single-network summaries ------------------------------------------------

network_summary <- function(A, threshold = EDGE_THRESHOLD) {
  n_edges <- sum(A > threshold)
  data.table(event_types = nrow(A), possible_pairs = length(A), edges = n_edges,
             density_pct = round(100 * n_edges / length(A), 1),
             median_strength = round(median(A[A > threshold]), 4),
             max_strength = round(max(A), 4))
}

edge_list <- function(A) {
  data.table(cause = rep(rownames(A), times = ncol(A)),
             effect = rep(colnames(A), each = nrow(A)),
             strength = as.vector(A))[, self := cause == effect]
}

top_edges <- function(A, k = 10, exclude_self = FALSE) {
  e <- edge_list(A)
  if (exclude_self) e <- e[self == FALSE]
  e <- e[order(-strength)][seq_len(k)]
  e[, rank := .I]
  e[, .(rank, cause, effect, strength = round(strength, 4))]
}

degree_table <- function(A, threshold = EDGE_THRESHOLD) {
  data.table(event = rownames(A),
             in_degree = colSums(A > threshold), out_degree = rowSums(A > threshold),
             weighted_in = round(colSums(A), 3), weighted_out = round(rowSums(A), 3),
             self_excitation = round(diag(A), 3))[
    , centrality := round((weighted_in + weighted_out) / 2, 3)][order(-centrality)]
}

# ---- Plots ---------------------------------------------------------------------

plot_heatmap <- function(A, mark_top = 10, legend = "Excitation\nstrength") {
  e <- edge_list(A)
  lv <- sort(rownames(A))
  e[, `:=`(cause = factor(cause, levels = rev(lv)), effect = factor(effect, levels = lv))]
  top <- e[order(-strength)][seq_len(mark_top)]
  ggplot(e, aes(effect, cause, fill = strength)) +
    geom_tile(colour = "grey80", linewidth = 0.1) +
    scale_fill_gradientn(colours = c("white", "#fcae91", "#fb6a4a", "#cb181d", "#67000d"),
                         name = legend) +
    geom_point(data = top, shape = 8, colour = "#00bfc4", size = 2) +
    labs(x = "Target disease", y = "Source disease") +
    coord_fixed() + theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank())
}

plot_difference <- function(D, legend) {
  e <- edge_list(D)
  lv <- sort(rownames(D))
  e[, `:=`(cause = factor(cause, levels = rev(lv)), effect = factor(effect, levels = lv))]
  lim <- max(abs(D))
  ggplot(e, aes(effect, cause, fill = strength)) +
    geom_tile(colour = "grey80", linewidth = 0.1) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0,
                         limits = c(-lim, lim), name = legend) +
    labs(x = "Target disease", y = "Source disease") +
    coord_fixed() + theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1), panel.grid = element_blank())
}

#' Directed network split at the 75th percentile of edge strength (graphical only)
plot_network <- function(A, threshold = EDGE_THRESHOLD, seed = 1) {
  e <- edge_list(A)[strength > threshold]
  cut <- quantile(e$strength, 0.75)
  nodes <- data.frame(name = rownames(A), weighted_in = colSums(A))
  set.seed(seed)
  layout_graph <- graph_from_data_frame(e[self == FALSE], vertices = nodes)
  xy <- create_layout(layout_graph, layout = "fr")[, c("x", "y")]
  panel <- function(sub, title) {
    g <- graph_from_data_frame(sub, vertices = nodes)
    ggraph(g, layout = "manual", x = xy$x, y = xy$y) +
      geom_edge_fan(aes(width = strength), alpha = 0.35, colour = "grey40",
                    arrow = arrow(length = unit(1.5, "mm"), type = "closed"), end_cap = circle(2, "mm")) +
      geom_edge_loop(aes(width = strength), alpha = 0.35, colour = "grey40") +
      geom_node_point(aes(size = weighted_in), colour = "#7b3294") +
      geom_node_text(aes(label = name), size = 2.2, repel = TRUE) +
      scale_edge_width(range = c(0.1, 1.5), limits = range(e$strength), name = "Excitation\nstrength") +
      scale_size(range = c(1, 6), name = "Weighted\nin-degree") +
      ggtitle(title) + theme_void(base_size = 9)
  }
  panel(e[strength <= cut], sprintf("(A) Weak connections (strength <= %.2f, n = %d)", cut, sum(e$strength <= cut))) +
    panel(e[strength > cut], sprintf("(B) Strong connections (strength > %.2f, n = %d)", cut, sum(e$strength > cut))) +
    plot_layout(guides = "collect")
}

#' Excitation kernels of the k strongest non-self relationships
plot_kernels <- function(model, k = 6) {
  top <- top_edges(model$A, k, exclude_self = TRUE)
  lags <- seq_len(dim(model$Phi)[2]) - 1
  d <- rbindlist(lapply(seq_len(k), function(i) {
    ci <- match(top$cause[i], model$names); ei <- match(top$effect[i], model$names)
    data.table(lag = lags, value = model$Phi[ci, , ei],
               panel = sprintf("%s -> %s\n(strength = %.3f)", top$cause[i], top$effect[i], top$strength[i]))
  }))
  d[, panel := factor(panel, levels = unique(panel))]
  ggplot(d, aes(lag, value)) +
    geom_area(fill = "#2166ac", alpha = 0.2) + geom_line(colour = "#2166ac") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~panel, ncol = 3) +
    scale_x_continuous(breaks = seq(0, 120, 24)) +
    labs(x = "Lag (months)", y = "Excitation kernel value") + theme_minimal(base_size = 9)
}

# ---- Comparison of two networks over the same event types --------------------

#' Element-wise comparison (B minus A). Substantially altered pairs: absolute
#' difference above the 95th percentile of the non-zero absolute differences.
#' Wilcoxon signed-rank test (two-sided) on the element-wise differences.
compare_networks <- function(A, B, pct = 0.95) {
  stopifnot(identical(rownames(A), rownames(B)), identical(colnames(A), colnames(B)))
  D <- B - A
  d <- as.vector(D)
  nz <- abs(d)[abs(d) > 0]
  q <- unname(quantile(nz, pct))
  wt <- wilcox.test(d, mu = 0, alternative = "two.sided", exact = FALSE)
  e <- edge_list(D)[, .(cause, effect, difference = strength)]
  e[, `:=`(strength_A = as.vector(A), strength_B = as.vector(B))]
  stats <- data.table(
    pearson_r = cor(as.vector(A), as.vector(B)),
    mean_difference = mean(d), mean_abs_difference = mean(abs(d)),
    n_positive = sum(d > 0), n_negative = sum(d < 0), n_zero = sum(d == 0),
    signed_rank_V = unname(wt$statistic), signed_rank_n = sum(d != 0), signed_rank_p = wt$p.value,
    threshold_abs_difference = q, n_altered = sum(abs(d) > q))
  list(diff = D, stats = stats,
       altered = e[abs(difference) > q][order(-abs(difference))])
}

#' Shared risk-factor attribution. For each CMD pair (cause -> effect), rank the
#' risk factors by their excitation of the cause and of the effect (rank 1 =
#' strongest); score = mean of the two ranks; the pair's exposure summary is the
#' minimum score over risk factors. Altered vs unaltered pairs are compared with
#' a two-sided Wilcoxon rank-sum test.
risk_attribution <- function(model_adj, altered) {
  map <- model_adj$mapping
  cmd <- map$event_name[map$group == "CMD"]; risk <- map$event_name[map$group == "RISK"]
  R <- model_adj$A[risk, cmd, drop = FALSE]          # risk factor -> CMD
  pairs <- CJ(cause = cmd, effect = cmd, sorted = FALSE)
  res <- pairs[, {
    ri <- rank(-abs(R[, cause]), ties.method = "average")
    rj <- rank(-abs(R[, effect]), ties.method = "average")
    s <- (ri + rj) / 2
    .(min_score = min(s), top_risk_factor = risk[which.min(s)])
  }, by = .(cause, effect)]
  res[, altered := paste(cause, effect) %in% paste(altered$cause, altered$effect)]
  wt <- wilcox.test(res[altered == TRUE, min_score], res[altered == FALSE, min_score],
                    alternative = "two.sided", exact = FALSE)
  list(pairs = res,
       test = data.table(n_altered = sum(res$altered), n_unaltered = sum(!res$altered),
                         median_altered = median(res[altered == TRUE, min_score]),
                         median_unaltered = median(res[altered == FALSE, min_score]),
                         rank_sum_W = unname(wt$statistic), p = wt$p.value))
}

#' Overlap of the top-p% strongest edges of two networks, and edge presence
top_overlap <- function(A, B, p = 0.05, threshold = EDGE_THRESHOLD) {
  k <- round(p * length(A))
  ea <- edge_list(A)[order(-strength)][seq_len(k), paste(cause, effect, sep = " -> ")]
  eb <- edge_list(B)[order(-strength)][seq_len(k), paste(cause, effect, sep = " -> ")]
  shared <- intersect(ea, eb)
  pa <- as.vector(A > threshold); pb <- as.vector(B > threshold)
  data.table(top_k = k, n_shared_top = length(shared),
             overlap_proportion = length(shared) / k,
             jaccard_top = length(shared) / length(union(ea, eb)),
             edges_any = sum(pa | pb), edges_shared = sum(pa & pb),
             edges_only_A = sum(pa & !pb), edges_only_B = sum(!pa & pb))
}
