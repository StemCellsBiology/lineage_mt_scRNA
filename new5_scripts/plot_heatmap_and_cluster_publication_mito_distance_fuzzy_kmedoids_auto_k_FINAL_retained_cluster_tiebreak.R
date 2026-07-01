#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cluster)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 15) {
  stop(
    "Usage: Rscript plot_heatmap_and_cluster_publication_mito_distance_fuzzy_kmedoids_auto_k.R ",
    "<filtered_long.tsv.gz> <out_prefix> <cluster_out_base> ",
    "[coverage_threshold=100] [min_shared_variants=1] [na_distance_action=one] ",
    "[require_alt_signal_cell=1] [min_alt_for_cell=2] ",
    "[auto_k_min=2] [auto_k_max=30] ",
    "[assignment_prob_threshold=0.95] [min_cluster_size=2] ",
    "[min_retained_clusters=2] [min_assigned_fraction=0] [membership_exponent=2]"
  )
}

filtered_long_file <- args[1]
out_prefix <- args[2]
cluster_out_base <- args[3]

coverage_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 100
min_shared_variants <- if (length(args) >= 5) as.integer(args[5]) else 1L
na_distance_action <- if (length(args) >= 6) args[6] else "one"
require_alt_signal_cell <- if (length(args) >= 7) as.integer(args[7]) else 1L
min_alt_for_cell <- if (length(args) >= 8) as.numeric(args[8]) else 2
auto_k_min <- if (length(args) >= 9) as.integer(args[9]) else 2L
auto_k_max <- if (length(args) >= 10) as.integer(args[10]) else 30L
assignment_prob_threshold <- if (length(args) >= 11) as.numeric(args[11]) else 0.95
min_cluster_size <- if (length(args) >= 12) as.integer(args[12]) else 2L
min_retained_clusters <- if (length(args) >= 13) as.integer(args[13]) else 2L
min_assigned_fraction <- if (length(args) >= 14) as.numeric(args[14]) else 0
membership_exponent <- if (length(args) >= 15) as.numeric(args[15]) else 2

if (is.na(coverage_threshold)) stop("coverage_threshold must be numeric.")
if (is.na(min_shared_variants) || min_shared_variants < 1) stop("min_shared_variants must be integer >= 1.")
if (!na_distance_action %in% c("max", "one", "stop", "drop")) stop("na_distance_action must be max, one, stop, or drop.")
if (is.na(require_alt_signal_cell) || !require_alt_signal_cell %in% c(0L, 1L)) stop("require_alt_signal_cell must be 0 or 1.")
if (is.na(min_alt_for_cell) || min_alt_for_cell < 0) stop("min_alt_for_cell must be numeric >= 0.")
if (is.na(auto_k_min) || auto_k_min < 2) stop("auto_k_min must be integer >= 2.")
if (is.na(auto_k_max) || auto_k_max < auto_k_min) stop("auto_k_max must be >= auto_k_min.")
if (is.na(assignment_prob_threshold) || assignment_prob_threshold <= 0 || assignment_prob_threshold > 1) {
  stop("assignment_prob_threshold must be in (0, 1].")
}
if (is.na(min_cluster_size) || min_cluster_size < 1) stop("min_cluster_size must be integer >= 1.")
if (is.na(min_retained_clusters) || min_retained_clusters < 1) stop("min_retained_clusters must be integer >= 1.")
if (is.na(min_assigned_fraction) || min_assigned_fraction < 0 || min_assigned_fraction > 1) {
  stop("min_assigned_fraction must be in [0, 1].")
}
if (is.na(membership_exponent) || membership_exponent <= 1) stop("membership_exponent must be > 1.")

message("============================================================")
message("plot_heatmap_and_cluster_publication_mito_distance_fuzzy_kmedoids_auto_k.R")
message("Input filtered long table: ", filtered_long_file)
message("Output prefix: ", out_prefix)
message("Cluster output base: ", cluster_out_base)
message("coverage_threshold: ", coverage_threshold)
message("coverage rule: depth > coverage_threshold")
message("min_shared_variants: ", min_shared_variants)
message("na_distance_action: ", na_distance_action)
message("require_alt_signal_cell: ", require_alt_signal_cell)
message("min_alt_for_cell: ", min_alt_for_cell)
message("auto_k_min: ", auto_k_min)
message("auto_k_max: ", auto_k_max)
message("assignment_prob_threshold: ", assignment_prob_threshold)
message("min_cluster_size: ", min_cluster_size)
message("min_retained_clusters: ", min_retained_clusters)
message("min_assigned_fraction: ", min_assigned_fraction)
message("membership_exponent: ", membership_exponent)
message("============================================================")

require_cols <- function(df, cols, description) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(description, " missing required column(s): ", paste(missing, collapse = ", "))
  }
}

write_matrix_tsv_gz <- function(mat, file) {
  out <- data.frame(barcode = rownames(mat), mat, check.names = FALSE)
  fwrite(out, file, sep = "\t")
}

safe_silhouette_summary <- function(cluster_labels, dist_matrix) {
  keep <- !is.na(cluster_labels) & cluster_labels != "unassigned"
  cluster_labels <- cluster_labels[keep]

  if (sum(keep) < 2 || length(unique(cluster_labels)) < 2) {
    return(list(
      mean = NA_real_,
      median = NA_real_,
      n_negative = NA_integer_
    ))
  }

  d_sub <- dist_matrix[keep, keep, drop = FALSE]
  sil <- cluster::silhouette(as.integer(factor(cluster_labels)), as.dist(d_sub))
  sil_values <- sil[, 3]

  list(
    mean = mean(sil_values, na.rm = TRUE),
    median = median(sil_values, na.rm = TRUE),
    n_negative = as.integer(sum(sil_values < 0, na.rm = TRUE))
  )
}

if (!file.exists(filtered_long_file) || file.info(filtered_long_file)$size == 0) {
  stop("Input filtered long table does not exist or is empty: ", filtered_long_file)
}

long <- fread(filtered_long_file, data.table = FALSE)
require_cols(long, c("barcode", "variant", "depth", "alt_count", "vaf"), "filtered long table")

long$barcode <- as.character(long$barcode)
long$variant <- as.character(long$variant)
long$depth <- as.numeric(long$depth)
long$alt_count <- as.numeric(long$alt_count)
long$vaf <- as.numeric(long$vaf)

long <- long[
  !is.na(long$barcode) & long$barcode != "" &
    !is.na(long$variant) & long$variant != "" &
    !is.na(long$depth) &
    !is.na(long$alt_count) &
    !is.na(long$vaf),
  ,
  drop = FALSE
]

if (nrow(long) == 0) stop("No valid rows after cleaning.")

message("Rows after cleaning: ", nrow(long))
message("Unique cells: ", length(unique(long$barcode)))
message("Unique variants: ", length(unique(long$variant)))

alt_signal_long <- long[
  long$depth > coverage_threshold &
    long$alt_count >= min_alt_for_cell &
    long$vaf > 0,
  ,
  drop = FALSE
]

fwrite(
  data.frame(
    metric = c(
      "coverage_threshold",
      "coverage_rule",
      "min_alt_for_cell",
      "n_rows_alt_signal",
      "n_cells_with_alt_signal",
      "n_variants_with_alt_signal"
    ),
    value = c(
      as.character(coverage_threshold),
      "depth > coverage_threshold",
      as.character(min_alt_for_cell),
      as.character(nrow(alt_signal_long)),
      as.character(length(unique(alt_signal_long$barcode))),
      as.character(length(unique(alt_signal_long$variant)))
    )
  ),
  paste0(out_prefix, ".alt_signal_cell_qc.tsv"),
  sep = "\t"
)

agg <- as.data.table(long)[
  ,
  .(
    vaf = mean(vaf, na.rm = TRUE),
    depth = max(depth, na.rm = TRUE),
    alt_count = max(alt_count, na.rm = TRUE)
  ),
  by = .(barcode, variant)
]

af_wide <- dcast(agg, barcode ~ variant, value.var = "vaf", fill = NA_real_)
cov_wide <- dcast(agg, barcode ~ variant, value.var = "depth", fill = NA_real_)

af_wide <- as.data.frame(af_wide[order(af_wide$barcode), , drop = FALSE])
cov_wide <- as.data.frame(cov_wide[order(cov_wide$barcode), , drop = FALSE])

if (!identical(af_wide$barcode, cov_wide$barcode)) stop("AF/COV barcode mismatch.")
if (!identical(colnames(af_wide), colnames(cov_wide))) stop("AF/COV column mismatch.")

variants <- setdiff(colnames(af_wide), "barcode")

AF <- as.matrix(af_wide[, variants, drop = FALSE])
COV <- as.matrix(cov_wide[, variants, drop = FALSE])
rownames(AF) <- af_wide$barcode
rownames(COV) <- cov_wide$barcode
mode(AF) <- "numeric"
mode(COV) <- "numeric"

cells_with_alt_signal <- unique(long$barcode[
  long$depth > coverage_threshold &
    long$alt_count >= min_alt_for_cell &
    long$vaf > 0
])

if (require_alt_signal_cell == 1) {
  keep_alt_cells <- rownames(AF) %in% cells_with_alt_signal
  if (sum(keep_alt_cells) < 2) stop("Fewer than 2 cells remain after ALT-signal cell filter.")
  AF <- AF[keep_alt_cells, , drop = FALSE]
  COV <- COV[keep_alt_cells, , drop = FALSE]
}

covered_by_cell <- rowSums(COV > coverage_threshold, na.rm = TRUE)
covered_by_variant <- colSums(COV > coverage_threshold, na.rm = TRUE)

fwrite(
  data.frame(
    barcode = names(covered_by_cell),
    n_variants_covered_above_threshold = as.integer(covered_by_cell),
    has_alt_signal_cell_filter = names(covered_by_cell) %in% cells_with_alt_signal
  ),
  paste0(out_prefix, ".coverage_per_cell_qc.tsv"),
  sep = "\t"
)

fwrite(
  data.frame(
    variant = names(covered_by_variant),
    n_cells_covered_above_threshold = as.integer(covered_by_variant)
  ),
  paste0(out_prefix, ".coverage_per_variant_qc.tsv"),
  sep = "\t"
)

cells_with_coverage <- covered_by_cell >= 1
if (sum(cells_with_coverage) < 2) stop("Fewer than 2 cells have at least one covered variant.")
AF <- AF[cells_with_coverage, , drop = FALSE]
COV <- COV[cells_with_coverage, , drop = FALSE]

variants_with_coverage <- colSums(COV > coverage_threshold, na.rm = TRUE) >= 1
if (sum(variants_with_coverage) < 1) stop("No variants covered above threshold in remaining cells.")
AF <- AF[, variants_with_coverage, drop = FALSE]
COV <- COV[, variants_with_coverage, drop = FALSE]

n_cells <- nrow(AF)
n_variants <- ncol(AF)

message("Cells used for D calculation: ", n_cells)
message("Variants used for D calculation: ", n_variants)

D <- matrix(
  NA_real_,
  nrow = n_cells,
  ncol = n_cells,
  dimnames = list(rownames(AF), rownames(AF))
)

N_shared <- matrix(
  0L,
  nrow = n_cells,
  ncol = n_cells,
  dimnames = list(rownames(AF), rownames(AF))
)

covered <- COV > coverage_threshold & !is.na(COV) & !is.na(AF)

for (i in seq_len(n_cells)) {
  if (i == 1 || i %% 100 == 0 || i == n_cells) {
    message("Processing cell ", i, " / ", n_cells)
  }
  D[i, i] <- 0
  N_shared[i, i] <- sum(covered[i, ])
  if (i < n_cells) {
    for (j in (i + 1):n_cells) {
      mask <- covered[i, ] & covered[j, ]
      denom <- sum(mask)
      N_shared[i, j] <- denom
      N_shared[j, i] <- denom
      if (denom >= min_shared_variants) {
        dij <- mean(sqrt(abs(AF[i, mask] - AF[j, mask])), na.rm = TRUE)
        D[i, j] <- dij
        D[j, i] <- dij
      }
    }
  }
}

diag(D) <- 0

Kmito <- 1 - D
diag(Kmito) <- 1

write_matrix_tsv_gz(D, paste0(out_prefix, ".publication_mito_distance.D.tsv.gz"))
write_matrix_tsv_gz(Kmito, paste0(out_prefix, ".publication_mito_relatedness.Kmito.tsv.gz"))
write_matrix_tsv_gz(N_shared, paste0(out_prefix, ".publication_mito_shared_variants.N.tsv.gz"))

D_clust <- D
offdiag <- row(D_clust) != col(D_clust)
n_na_offdiag <- sum(is.na(D_clust[offdiag]))
n_total_offdiag <- sum(offdiag)
cells_dropped <- character(0)

if (n_na_offdiag > 0) {
  if (na_distance_action == "stop") stop("NA distances present.")
  if (na_distance_action == "max") {
    finite_vals <- D_clust[is.finite(D_clust)]
    if (length(finite_vals) == 0) stop("No finite distances available.")
    D_clust[is.na(D_clust)] <- max(finite_vals, na.rm = TRUE)
    diag(D_clust) <- 0
  }
  if (na_distance_action == "one") {
    D_clust[is.na(D_clust)] <- 1
    diag(D_clust) <- 0
  }
  if (na_distance_action == "drop") {
    repeat {
      na_counts <- rowSums(is.na(D_clust))
      if (max(na_counts) == 0) break
      cell_to_drop <- names(which.max(na_counts))
      cells_dropped <- c(cells_dropped, cell_to_drop)
      D_clust <- D_clust[
        rownames(D_clust) != cell_to_drop,
        colnames(D_clust) != cell_to_drop,
        drop = FALSE
      ]
      if (nrow(D_clust) < 2) stop("Fewer than 2 cells remain after dropping.")
    }
  }
}

D_clust <- (D_clust + t(D_clust)) / 2
diag(D_clust) <- 0

if (any(!is.finite(D_clust))) stop("Clustering-ready D still contains non-finite values.")

write_matrix_tsv_gz(D_clust, paste0(out_prefix, ".publication_mito_distance.clustering_ready.D.tsv.gz"))

fwrite(
  data.frame(
    barcode = c(rownames(D_clust), cells_dropped),
    used_for_clustering = c(rep(TRUE, nrow(D_clust)), rep(FALSE, length(cells_dropped)))
  ),
  paste0(out_prefix, ".publication_mito_distance.cells_used_for_clustering.tsv"),
  sep = "\t"
)

n_clust_cells <- nrow(D_clust)
candidate_k_max <- min(auto_k_max, floor(n_clust_cells / 2))

if (candidate_k_max < auto_k_min) {
  candidate_k <- 2L
  if (candidate_k >= n_clust_cells) stop("Not enough cells for k=2 after filtering.")
  warning("Candidate k range collapsed. Testing only k=2.")
} else {
  candidate_k <- seq(auto_k_min, candidate_k_max)
}

message("Candidate K values: ", paste(candidate_k, collapse = ", "))

fanny_results <- list()

auto_results <- rbindlist(lapply(candidate_k, function(kk) {
  message("Running fuzzy k-medoids for k=", kk)

  fit <- tryCatch(
    cluster::fanny(
      D_clust,
      k = kk,
      diss = TRUE,
      memb.exp = membership_exponent
    ),
    error = function(e) {
      message("fanny failed for k=", kk, ": ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(data.frame(
      k = kk,
      fanny_ok = FALSE,
      mean_silhouette_assigned = NA_real_,
      median_silhouette_assigned = NA_real_,
      n_cells_total = nrow(D_clust),
      n_cells_initially_assigned = NA_integer_,
      n_cells_assigned = NA_integer_,
      n_cells_unassigned = NA_integer_,
      assigned_fraction = NA_real_,
      n_retained_clusters = NA_integer_,
      n_removed_small_clusters = NA_integer_,
      n_cells_moved_from_small_clusters_to_unassigned = NA_integer_,
      min_retained_cluster_size = NA_integer_,
      max_retained_cluster_size = NA_integer_,
      n_negative_silhouette_assigned_cells = NA_integer_,
      normalized_dunn = NA_real_,
      partition_coefficient = NA_real_,
      passes_min_retained_clusters = FALSE,
      passes_min_assigned_fraction = FALSE
    ))
  }

  membership <- fit$membership
  max_prob <- apply(membership, 1, max)
  raw_cluster <- apply(membership, 1, which.max)

  initially_assigned <- max_prob >= assignment_prob_threshold
  final_cluster <- rep("unassigned", length(raw_cluster))
  final_cluster[initially_assigned] <- paste0("C", raw_cluster[initially_assigned])
  names(final_cluster) <- rownames(D_clust)

  initial_cluster_sizes <- table(final_cluster[final_cluster != "unassigned"])
  small_clusters <- names(initial_cluster_sizes)[initial_cluster_sizes < min_cluster_size]
  moved_small <- final_cluster %in% small_clusters
  final_cluster[moved_small] <- "unassigned"

  retained_cluster_sizes <- table(final_cluster[final_cluster != "unassigned"])
  n_retained_clusters <- length(retained_cluster_sizes)
  n_assigned <- sum(final_cluster != "unassigned")
  n_unassigned <- sum(final_cluster == "unassigned")
  assigned_fraction <- n_assigned / length(final_cluster)

  sil_summary <- safe_silhouette_summary(final_cluster, D_clust)

  fanny_results[[as.character(kk)]] <<- list(
    fit = fit,
    membership = membership,
    max_prob = max_prob,
    raw_cluster = raw_cluster,
    initially_assigned = initially_assigned,
    final_cluster = final_cluster,
    small_clusters = small_clusters,
    moved_small = moved_small
  )

  data.frame(
    k = kk,
    fanny_ok = TRUE,
    mean_silhouette_assigned = sil_summary$mean,
    median_silhouette_assigned = sil_summary$median,
    n_cells_total = length(final_cluster),
    n_cells_initially_assigned = sum(initially_assigned),
    n_cells_assigned = n_assigned,
    n_cells_unassigned = n_unassigned,
    assigned_fraction = assigned_fraction,
    n_retained_clusters = as.integer(n_retained_clusters),
    n_removed_small_clusters = as.integer(length(small_clusters)),
    n_cells_moved_from_small_clusters_to_unassigned = as.integer(sum(moved_small)),
    min_retained_cluster_size = ifelse(n_retained_clusters > 0, as.integer(min(retained_cluster_sizes)), NA_integer_),
    max_retained_cluster_size = ifelse(n_retained_clusters > 0, as.integer(max(retained_cluster_sizes)), NA_integer_),
    n_negative_silhouette_assigned_cells = sil_summary$n_negative,
    normalized_dunn = as.numeric(fit$coeff["normalized dunn"]),
    partition_coefficient = as.numeric(fit$coeff["partition coefficient"]),
    passes_min_retained_clusters = n_retained_clusters >= min_retained_clusters,
    passes_min_assigned_fraction = assigned_fraction >= min_assigned_fraction
  )
}))

auto_results$rank_silhouette_assigned <- rank(
  -auto_results$mean_silhouette_assigned,
  ties.method = "first",
  na.last = "keep"
)

auto_results$selected <- FALSE

fwrite(auto_results, paste0(out_prefix, ".auto_k_selection.tsv"), sep = "\t")

# ------------------------------------------------------------
# Select best K after fuzzy k-medoids QC
# ------------------------------------------------------------
#
# For every candidate K, membership filtering and small-cluster
# filtering have already been applied above. Silhouette was then
# calculated only on retained assigned cells.
#
# A K is eligible only if:
#   1. fanny succeeded,
#   2. mean_silhouette_assigned is finite,
#   3. the number of retained clusters is >= min_retained_clusters,
#   4. assigned_fraction is >= min_assigned_fraction.
#
# Important: by default min_assigned_fraction = 0, so assigned_fraction
# is not used as a hard threshold unless explicitly set by the pipeline/user.
# In the default mode it is used only as a tie-breaker after silhouette.
#
# Ranking among eligible K values:
#   1. highest mean_silhouette_assigned,
#   2. if tied, highest assigned_fraction,
#   3. if tied, highest n_cells_assigned,
#   4. if tied, lowest n_retained_clusters,
#   5. if tied, lowest k.

valid_auto <- auto_results[
  fanny_ok == TRUE &
    is.finite(mean_silhouette_assigned) &
    passes_min_retained_clusters == TRUE &
    passes_min_assigned_fraction == TRUE,
  ,
  drop = FALSE
]

if (nrow(valid_auto) == 0) {
  stop(
    "No valid K passed all fuzzy STEP 12B QC filters: ",
    "fanny_ok, finite assigned-cell silhouette, ",
    "minimum retained clusters, and minimum assigned fraction. ",
    "No fuzzy cluster file will be produced for this sample."
  )
}

valid_auto <- valid_auto[
  order(
    -mean_silhouette_assigned,
    -assigned_fraction,
    -n_cells_assigned,
    n_retained_clusters,
    k
  ),
  ,
  drop = FALSE
]

best_row <- valid_auto[1, , drop = FALSE]
best_k <- as.integer(best_row$k)

auto_results$selected <- FALSE
auto_results$selected[auto_results$k == best_k] <- TRUE

auto_results$selection_rank <- NA_integer_
auto_results$selection_rank[
  match(valid_auto$k, auto_results$k)
] <- seq_len(nrow(valid_auto))

auto_results$best_k_selection_rule <- paste(
  "eligible: fanny_ok, finite mean_silhouette_assigned,",
  "n_retained_clusters >= min_retained_clusters,",
  "assigned_fraction >= min_assigned_fraction",
  "(default 0; no extra fraction threshold unless explicitly set);",
  "rank: highest silhouette, highest assigned_fraction,",
  "highest n_cells_assigned, smallest k"
)

fwrite(auto_results, paste0(out_prefix, ".auto_k_selection.tsv"), sep = "\t")

writeLines(as.character(best_k), paste0(out_prefix, ".best_k.txt"))

best <- fanny_results[[as.character(best_k)]]
if (is.null(best)) stop("Internal error: missing fanny result for best_k=", best_k)

membership <- best$membership
max_prob <- best$max_prob
raw_cluster <- best$raw_cluster
initially_assigned <- best$initially_assigned
final_cluster <- best$final_cluster
small_clusters <- best$small_clusters
moved_small <- best$moved_small

all_assignments <- data.frame(
  barcode = rownames(D_clust),
  mt_cluster_raw = paste0("C", raw_cluster),
  max_membership = as.numeric(max_prob),
  initially_assigned_by_membership = initially_assigned,
  moved_to_unassigned_due_to_small_cluster = moved_small,
  mt_cluster = as.character(final_cluster),
  assignment_status = ifelse(final_cluster == "unassigned", "unassigned", "assigned"),
  stringsAsFactors = FALSE
)

clusters_df <- all_assignments[
  all_assignments$assignment_status == "assigned",
  c("barcode", "mt_cluster"),
  drop = FALSE
]

if (nrow(clusters_df) == 0) {
  stop("best_k produced zero assigned cells after filtering small clusters.")
}

clusters_out <- paste0(cluster_out_base, ".k", best_k, ".tsv")
fwrite(clusters_df, clusters_out, sep = "\t")

fwrite(
  all_assignments,
  paste0(out_prefix, ".fuzzy_kmedoids.assignments.k", best_k, ".tsv"),
  sep = "\t"
)

membership_df <- data.frame(
  barcode = rownames(membership),
  membership,
  check.names = FALSE
)
fwrite(
  membership_df,
  paste0(out_prefix, ".fuzzy_kmedoids.membership.k", best_k, ".tsv"),
  sep = "\t"
)

retained_cluster_sizes <- as.data.frame(table(clusters_df$mt_cluster))
colnames(retained_cluster_sizes) <- c("mt_cluster", "n_cells")
fwrite(
  retained_cluster_sizes,
  paste0(out_prefix, ".fuzzy_kmedoids.cluster_sizes.k", best_k, ".tsv"),
  sep = "\t"
)

if (length(small_clusters) > 0) {
  fwrite(
    data.frame(removed_small_cluster = small_clusters),
    paste0(out_prefix, ".fuzzy_kmedoids.removed_small_clusters.k", best_k, ".tsv"),
    sep = "\t"
  )
} else {
  fwrite(
    data.frame(removed_small_cluster = character(0)),
    paste0(out_prefix, ".fuzzy_kmedoids.removed_small_clusters.k", best_k, ".tsv"),
    sep = "\t"
  )
}

heatmap_pdf <- paste0(out_prefix, ".fuzzy_kmedoids.distance_heatmap.k", best_k, ".pdf")

if (requireNamespace("pheatmap", quietly = TRUE)) {
  ann <- data.frame(mt_cluster = all_assignments$mt_cluster, row.names = all_assignments$barcode)
  hc_tmp <- hclust(as.dist(D_clust), method = "average")
  pdf(heatmap_pdf, width = 10, height = 10)
  pheatmap::pheatmap(
    D_clust,
    clustering_distance_rows = as.dist(D_clust),
    clustering_distance_cols = as.dist(D_clust),
    clustering_method = "average",
    annotation_row = ann,
    annotation_col = ann,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste0("Fuzzy k-medoids mtDNA distance, best_k=", best_k)
  )
  dev.off()
} else {
  hc_tmp <- hclust(as.dist(D_clust), method = "average")
  ord <- hc_tmp$order
  D_ord <- D_clust[ord, ord, drop = FALSE]
  hm_dt <- as.data.table(as.table(D_ord))
  colnames(hm_dt) <- c("cell_i", "cell_j", "distance")
  hm_dt$cell_i <- factor(hm_dt$cell_i, levels = rownames(D_ord))
  hm_dt$cell_j <- factor(hm_dt$cell_j, levels = colnames(D_ord))
  p <- ggplot(hm_dt, aes(x = cell_i, y = cell_j, fill = distance)) +
    geom_raster() +
    theme_void() +
    labs(title = paste0("Fuzzy k-medoids mtDNA distance, best_k=", best_k), fill = "D")
  ggsave(heatmap_pdf, p, width = 10, height = 10)
}

finite_D <- D[is.finite(D)]
finite_D_offdiag <- D[offdiag & is.finite(D)]

qc <- data.frame(
  input_file = filtered_long_file,
  out_prefix = out_prefix,
  cluster_out_base = cluster_out_base,
  clusters_out = clusters_out,
  best_k = best_k,
  best_mean_silhouette_assigned = as.numeric(best_row$mean_silhouette_assigned),
  best_median_silhouette_assigned = as.numeric(best_row$median_silhouette_assigned),
  best_n_cells_assigned = as.integer(best_row$n_cells_assigned),
  best_n_cells_unassigned = as.integer(best_row$n_cells_unassigned),
  best_assigned_fraction = as.numeric(best_row$assigned_fraction),
  best_n_retained_clusters = as.integer(best_row$n_retained_clusters),
  best_min_retained_cluster_size = as.integer(best_row$min_retained_cluster_size),
  best_k_selection_rule = paste(
    "eligible: fanny_ok, finite mean_silhouette_assigned,",
    "n_retained_clusters >= min_retained_clusters,",
    "assigned_fraction >= min_assigned_fraction",
    "(default 0; no extra fraction threshold unless explicitly set);",
    "rank: highest silhouette, highest assigned_fraction,",
    "highest n_cells_assigned, lowest n_retained_clusters,",
    "smallest k"
  ),
  coverage_threshold = coverage_threshold,
  coverage_rule = "depth > coverage_threshold",
  min_shared_variants = min_shared_variants,
  na_distance_action = na_distance_action,
  require_alt_signal_cell = require_alt_signal_cell,
  min_alt_for_cell = min_alt_for_cell,
  auto_k_min = auto_k_min,
  auto_k_max = auto_k_max,
  assignment_prob_threshold = assignment_prob_threshold,
  min_cluster_size = min_cluster_size,
  min_retained_clusters = min_retained_clusters,
  min_assigned_fraction = min_assigned_fraction,
  membership_exponent = membership_exponent,
  n_rows_long_after_basic_cleaning = nrow(long),
  n_rows_aggregated_barcode_variant = nrow(agg),
  n_cells_initial_after_basic_cleaning = length(unique(long$barcode)),
  n_variants_initial_after_basic_cleaning = length(unique(long$variant)),
  n_rows_alt_signal = nrow(alt_signal_long),
  n_cells_with_alt_signal = length(unique(alt_signal_long$barcode)),
  n_variants_with_alt_signal = length(unique(alt_signal_long$variant)),
  n_cells_used_for_distance = nrow(D),
  n_cells_used_for_clustering = nrow(D_clust),
  n_variants_used_for_distance = ncol(AF),
  n_cells_dropped_for_clustering = length(cells_dropped),
  n_offdiag_pairs = n_total_offdiag,
  n_offdiag_pairs_with_NA_distance = n_na_offdiag,
  fraction_offdiag_pairs_with_NA_distance = n_na_offdiag / n_total_offdiag,
  min_finite_D = ifelse(length(finite_D) > 0, min(finite_D, na.rm = TRUE), NA_real_),
  max_finite_D = ifelse(length(finite_D) > 0, max(finite_D, na.rm = TRUE), NA_real_),
  mean_finite_offdiag_D = ifelse(length(finite_D_offdiag) > 0, mean(finite_D_offdiag, na.rm = TRUE), NA_real_),
  median_finite_offdiag_D = ifelse(length(finite_D_offdiag) > 0, median(finite_D_offdiag, na.rm = TRUE), NA_real_)
)

fwrite(qc, paste0(out_prefix, ".fuzzy_kmedoids.qc.tsv"), sep = "\t")

message("============================================================")
message("DONE")
message("best_k: ", best_k)
message("assigned-only clusters_out: ", clusters_out)
message("full assignments: ", paste0(out_prefix, ".fuzzy_kmedoids.assignments.k", best_k, ".tsv"))
message("membership: ", paste0(out_prefix, ".fuzzy_kmedoids.membership.k", best_k, ".tsv"))
message("auto_k_selection: ", paste0(out_prefix, ".auto_k_selection.tsv"))
message("best_k file: ", paste0(out_prefix, ".best_k.txt"))
message("QC: ", paste0(out_prefix, ".fuzzy_kmedoids.qc.tsv"))
message("============================================================")
