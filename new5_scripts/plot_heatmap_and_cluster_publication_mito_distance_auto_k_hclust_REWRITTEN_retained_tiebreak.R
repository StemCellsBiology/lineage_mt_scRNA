#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cluster)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 11) {
  stop(
    "Usage: Rscript plot_heatmap_and_cluster_publication_mito_distance_auto_k.R ",
    "<filtered_long.tsv.gz> <out_prefix> <cluster_out_base> ",
    "[coverage_threshold=100] [min_shared_variants=1] [na_distance_action=one] ",
    "[require_alt_signal_cell=1] [min_alt_for_cell=2] ",
    "[auto_k_min=2] [auto_k_max=30] [auto_k_min_cluster_size=3]"
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
auto_k_min_cluster_size <- if (length(args) >= 11) as.integer(args[11]) else 3L

if (is.na(coverage_threshold)) stop("coverage_threshold must be numeric.")
if (is.na(min_shared_variants) || min_shared_variants < 1) stop("min_shared_variants must be integer >= 1.")
if (!na_distance_action %in% c("max", "one", "stop", "drop")) stop("na_distance_action must be max, one, stop, or drop.")
if (is.na(require_alt_signal_cell) || !require_alt_signal_cell %in% c(0L, 1L)) stop("require_alt_signal_cell must be 0 or 1.")
if (is.na(min_alt_for_cell) || min_alt_for_cell < 0) stop("min_alt_for_cell must be numeric >= 0.")
if (is.na(auto_k_min) || auto_k_min < 2) stop("auto_k_min must be integer >= 2.")
if (is.na(auto_k_max) || auto_k_max < auto_k_min) stop("auto_k_max must be >= auto_k_min.")
if (is.na(auto_k_min_cluster_size) || auto_k_min_cluster_size < 1) stop("auto_k_min_cluster_size must be integer >= 1.")

message("============================================================")
message("plot_heatmap_and_cluster_publication_mito_distance_auto_k.R")
message("Mode: hclust with small clusters moved to unassigned")
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
message("auto_k_min_cluster_size: ", auto_k_min_cluster_size)
message("small-cluster rule: clusters smaller than auto_k_min_cluster_size are moved to unassigned")
message(
  "best-k rule with finite silhouette: highest mean_silhouette_assigned; ",
  "ties resolved by highest assigned_fraction, highest n_cells_assigned, ",
  "lowest n_retained_clusters, then lowest k"
)
message(
  "fallback rule without finite silhouette: select retained candidate by ",
  "highest n_retained_clusters, highest n_cells_assigned, then lowest k"
)
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

safe_hclust_silhouette_summary <- function(cluster_labels, dist_matrix) {
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

order_cells_by_hclust_assignment <- function(assignments_df) {
  df <- assignments_df

  df$assignment_order <- ifelse(df$mt_cluster == "unassigned", 1L, 0L)

  df$cluster_number <- suppressWarnings(
    as.integer(gsub("^C", "", df$mt_cluster))
  )

  df$cluster_number[is.na(df$cluster_number)] <- Inf

  df <- df[
    order(
      df$assignment_order,
      df$cluster_number,
      df$barcode
    ),
    ,
    drop = FALSE
  ]

  df$barcode
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

dist_obj <- as.dist(D_clust)
hc <- hclust(dist_obj, method = "average")
saveRDS(hc, paste0(out_prefix, ".publication_mito_distance.hclust.rds"))

# ------------------------------------------------------------
# Save standard hclust dendrogram
# ------------------------------------------------------------

hclust_pdf <- paste0(
  out_prefix,
  ".publication_mito_distance.hclust.dendrogram.pdf"
)

pdf(hclust_pdf, width = 12, height = 8)
plot(
  hc,
  labels = FALSE,
  hang = -1,
  main = "Publication-like mtDNA distance hierarchical clustering",
  xlab = "Cells",
  sub = paste0("method = average; n = ", nrow(D_clust), " cells")
)
dev.off()

message("Saved standard hclust dendrogram:")
message(hclust_pdf)

n_clust_cells <- nrow(D_clust)
candidate_k_max <- min(auto_k_max, floor(n_clust_cells / 2))

if (candidate_k_max < auto_k_min) {
  candidate_k <- 2
  if (candidate_k > n_clust_cells) stop("Not enough cells for k=2 after filtering.")
  warning("Candidate k range collapsed. Testing only k=2.")
} else {
  candidate_k <- seq(auto_k_min, candidate_k_max)
}

message("Candidate k values: ", paste(candidate_k, collapse = ", "))

hclust_results <- list()

auto_results <- rbindlist(lapply(candidate_k, function(kk) {
  cl_raw <- cutree(hc, k = kk)
  tab_raw <- table(cl_raw)

  small_clusters <- names(tab_raw)[tab_raw < auto_k_min_cluster_size]

  cl_final <- paste0("C", as.integer(cl_raw))
  names(cl_final) <- names(cl_raw)

  moved_small <- as.character(cl_raw) %in% small_clusters
  cl_final[moved_small] <- "unassigned"

  retained_cluster_sizes <- table(cl_final[cl_final != "unassigned"])
  n_retained_clusters <- length(retained_cluster_sizes)
  n_assigned <- sum(cl_final != "unassigned")
  n_unassigned <- sum(cl_final == "unassigned")
  assigned_fraction <- n_assigned / length(cl_final)

  sil_summary <- safe_hclust_silhouette_summary(cl_final, D_clust)

  hclust_results[[as.character(kk)]] <<- list(
    cl_raw = cl_raw,
    cl_final = cl_final,
    small_clusters = small_clusters,
    moved_small = moved_small
  )

  data.frame(
    k = kk,
    mean_silhouette_assigned = sil_summary$mean,
    median_silhouette_assigned = sil_summary$median,
    n_cells_total = length(cl_final),
    n_cells_assigned = n_assigned,
    n_cells_unassigned = n_unassigned,
    assigned_fraction = assigned_fraction,
    n_raw_clusters = length(tab_raw),
    raw_min_cluster_size = as.integer(min(tab_raw)),
    raw_max_cluster_size = as.integer(max(tab_raw)),
    n_retained_clusters = as.integer(n_retained_clusters),
    n_removed_small_clusters = as.integer(length(small_clusters)),
    n_cells_moved_from_small_clusters_to_unassigned = as.integer(sum(moved_small)),
    min_retained_cluster_size = ifelse(
      n_retained_clusters > 0,
      as.integer(min(retained_cluster_sizes)),
      NA_integer_
    ),
    max_retained_cluster_size = ifelse(
      n_retained_clusters > 0,
      as.integer(max(retained_cluster_sizes)),
      NA_integer_
    ),
    n_negative_silhouette_assigned_cells = sil_summary$n_negative,
    passes_min_retained_clusters = n_retained_clusters >= 1
  )
}))

auto_results$rank_silhouette_assigned <- rank(
  -auto_results$mean_silhouette_assigned,
  ties.method = "first",
  na.last = "keep"
)

auto_results$selected <- FALSE

fwrite(auto_results, paste0(out_prefix, ".auto_k_selection.tsv"), sep = "\t")

valid_auto_with_silhouette <- auto_results[
  is.finite(mean_silhouette_assigned) &
    passes_min_retained_clusters == TRUE
]

valid_auto_retained <- auto_results[
  passes_min_retained_clusters == TRUE &
    n_cells_assigned > 0
]

if (nrow(valid_auto_with_silhouette) > 0) {

  # Normal mode:
  # finite silhouette exists for at least one candidate k.
  #
  # Selection rule:
  #   1. highest mean_silhouette_assigned
  #   2. if tied, highest assigned_fraction
  #   3. if tied, highest n_cells_assigned
  #   4. if tied, lowest n_retained_clusters
  #   5. if tied, lowest k

  valid_auto_with_silhouette <- valid_auto_with_silhouette[
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

  best_row <- valid_auto_with_silhouette[1, , drop = FALSE]

  best_k_selection_mode <- "ranked_mean_silhouette_assigned"

} else if (nrow(valid_auto_retained) > 0) {

  warning(
    "No candidate k had finite assigned-cell silhouette after moving small clusters to unassigned. ",
    "Selecting best k among candidates with at least one retained cluster. ",
    "This usually means that only one retained cluster remained after filtering small clusters."
  )

  # Fallback mode:
  # No finite silhouette exists.
  # Keep original fallback logic from the previous hclust script.
  valid_auto_retained <- valid_auto_retained[
    order(
      -n_retained_clusters,
      -n_cells_assigned,
      k
    ),
    ,
    drop = FALSE
  ]

  best_row <- valid_auto_retained[1, , drop = FALSE]

  best_k_selection_mode <- "retained_cluster_without_finite_silhouette"

} else {

  stop(
    "No valid k passed STEP 12 hclust QC filters after moving small clusters to unassigned: ",
    "no candidate retained any assigned cells. ",
    "No cluster file will be produced for this sample."
  )
}

best_k <- as.integer(best_row$k)

auto_results$selected <- FALSE
auto_results$selected[auto_results$k == best_k] <- TRUE

auto_results$selection_rank <- NA_integer_

if (best_k_selection_mode == "ranked_mean_silhouette_assigned") {
  auto_results$selection_rank[
    match(valid_auto_with_silhouette$k, auto_results$k)
  ] <- seq_len(nrow(valid_auto_with_silhouette))
} else if (best_k_selection_mode == "retained_cluster_without_finite_silhouette") {
  auto_results$selection_rank[
    match(valid_auto_retained$k, auto_results$k)
  ] <- seq_len(nrow(valid_auto_retained))
}

auto_results$best_k_selection_mode <- NA_character_
auto_results$best_k_selection_mode[auto_results$k == best_k] <- best_k_selection_mode

auto_results$best_k_selection_rule <- NA_character_
auto_results$best_k_selection_rule[auto_results$k == best_k] <- if (
  best_k_selection_mode == "ranked_mean_silhouette_assigned"
) {
  paste(
    "finite silhouette mode:",
    "highest mean_silhouette_assigned;",
    "ties resolved by highest assigned_fraction,",
    "highest n_cells_assigned,",
    "lowest n_retained_clusters,",
    "then lowest k"
  )
} else {
  paste(
    "fallback mode without finite silhouette:",
    "highest n_retained_clusters,",
    "highest n_cells_assigned,",
    "then lowest k"
  )
}

fwrite(auto_results, paste0(out_prefix, ".auto_k_selection.tsv"), sep = "\t")

writeLines(as.character(best_k), paste0(out_prefix, ".best_k.txt"))

# ------------------------------------------------------------
# Save hclust dendrogram with rectangles for selected best_k
# ------------------------------------------------------------

hclust_bestk_pdf <- paste0(
  out_prefix,
  ".publication_mito_distance.hclust.dendrogram.best_k",
  best_k,
  ".pdf"
)

pdf(hclust_bestk_pdf, width = 12, height = 8)
plot(
  hc,
  labels = FALSE,
  hang = -1,
  main = paste0("Publication-like mtDNA hclust, best_k = ", best_k),
  xlab = "Cells",
  sub = paste0("method = average; n = ", nrow(D_clust), " cells")
)
rect.hclust(hc, k = best_k, border = 2:(best_k + 1))
dev.off()

message("Saved hclust dendrogram with best_k rectangles:")
message(hclust_bestk_pdf)

best <- hclust_results[[as.character(best_k)]]
if (is.null(best)) stop("Internal error: missing hclust result for best_k=", best_k)

clusters_raw <- best$cl_raw
clusters_final <- best$cl_final
small_clusters <- best$small_clusters
moved_small <- best$moved_small

all_assignments <- data.frame(
  barcode = names(clusters_final),
  mt_cluster_raw = paste0("C", as.integer(clusters_raw)),
  moved_to_unassigned_due_to_small_cluster = moved_small,
  mt_cluster = as.character(clusters_final),
  assignment_status = ifelse(clusters_final == "unassigned", "unassigned", "assigned"),
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
  paste0(out_prefix, ".hclust.assignments.k", best_k, ".tsv"),
  sep = "\t"
)

if (length(small_clusters) > 0) {
  fwrite(
    data.frame(removed_small_cluster = small_clusters),
    paste0(out_prefix, ".hclust.removed_small_clusters.k", best_k, ".tsv"),
    sep = "\t"
  )
} else {
  fwrite(
    data.frame(removed_small_cluster = character(0)),
    paste0(out_prefix, ".hclust.removed_small_clusters.k", best_k, ".tsv"),
    sep = "\t"
  )
}

cluster_sizes <- as.data.frame(table(clusters_df$mt_cluster))
colnames(cluster_sizes) <- c("mt_cluster", "n_cells")
fwrite(cluster_sizes, paste0(out_prefix, ".publication_mito_distance.cluster_sizes.k", best_k, ".tsv"), sep = "\t")

all_cluster_sizes <- as.data.frame(table(all_assignments$mt_cluster))
colnames(all_cluster_sizes) <- c("mt_cluster", "n_cells")
fwrite(all_cluster_sizes, paste0(out_prefix, ".publication_mito_distance.cluster_sizes_with_unassigned.k", best_k, ".tsv"), sep = "\t")

heatmap_pdf <- paste0(out_prefix, ".publication_mito_distance.heatmap.k", best_k, ".pdf")

if (requireNamespace("pheatmap", quietly = TRUE)) {
  ann <- data.frame(mt_cluster = all_assignments$mt_cluster, row.names = all_assignments$barcode)
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
    main = paste0("Publication-like mtDNA distance, best_k=", best_k)
  )
  dev.off()
} else {
  ord <- hc$order
  D_ord <- D_clust[ord, ord, drop = FALSE]
  hm_dt <- as.data.table(as.table(D_ord))
  colnames(hm_dt) <- c("cell_i", "cell_j", "distance")
  hm_dt$cell_i <- factor(hm_dt$cell_i, levels = rownames(D_ord))
  hm_dt$cell_j <- factor(hm_dt$cell_j, levels = colnames(D_ord))
  p <- ggplot(hm_dt, aes(x = cell_i, y = cell_j, fill = distance)) +
    geom_raster() +
    theme_void() +
    labs(title = paste0("Publication-like mtDNA distance, best_k=", best_k), fill = "D")
  ggsave(heatmap_pdf, p, width = 10, height = 10)
}

# ------------------------------------------------------------
# Optional diagnostic heatmap ordered by final hclust assignment
# ------------------------------------------------------------

ordered_cells_hclust <- order_cells_by_hclust_assignment(all_assignments)
D_hclust_ordered <- D_clust[ordered_cells_hclust, ordered_cells_hclust, drop = FALSE]

distance_hclust_order_pdf <- paste0(
  out_prefix,
  ".hclust.distance_heatmap.ordered_by_assignment.k",
  best_k,
  ".pdf"
)

distance_long <- as.data.table(as.table(D_hclust_ordered))
colnames(distance_long) <- c("cell_i", "cell_j", "distance")

distance_long$cell_i <- factor(
  distance_long$cell_i,
  levels = ordered_cells_hclust
)

distance_long$cell_j <- factor(
  distance_long$cell_j,
  levels = ordered_cells_hclust
)

p_distance_hclust_order <- ggplot(
  distance_long,
  aes(x = cell_i, y = cell_j, fill = distance)
) +
  geom_raster() +
  theme_void() +
  labs(
    title = paste0("mtDNA distance ordered by hclust assignment, best_k = ", best_k),
    fill = "D"
  )

ggsave(
  distance_hclust_order_pdf,
  p_distance_hclust_order,
  width = 10,
  height = 10
)

message("Saved distance heatmap ordered by hclust assignment:")
message(distance_hclust_order_pdf)

finite_D <- D[is.finite(D)]
finite_D_offdiag <- D[offdiag & is.finite(D)]

best_k_selection_rule <- if (
  best_k_selection_mode == "ranked_mean_silhouette_assigned"
) {
  paste(
    "finite silhouette mode:",
    "highest mean_silhouette_assigned;",
    "ties resolved by highest assigned_fraction,",
    "highest n_cells_assigned,",
    "lowest n_retained_clusters,",
    "then lowest k"
  )
} else {
  paste(
    "fallback mode without finite silhouette:",
    "highest n_retained_clusters,",
    "highest n_cells_assigned,",
    "then lowest k"
  )
}

qc <- data.frame(
  input_file = filtered_long_file,
  out_prefix = out_prefix,
  cluster_out_base = cluster_out_base,
  clusters_out = clusters_out,
  assignments_out = paste0(out_prefix, ".hclust.assignments.k", best_k, ".tsv"),
  best_k = best_k,
  best_k_selection_mode = best_k_selection_mode,
  best_k_selection_rule = best_k_selection_rule,
  best_mean_silhouette_assigned = as.numeric(best_row$mean_silhouette_assigned),
  best_median_silhouette_assigned = as.numeric(best_row$median_silhouette_assigned),
  best_n_cells_assigned = as.integer(best_row$n_cells_assigned),
  best_n_cells_unassigned = as.integer(best_row$n_cells_unassigned),
  best_assigned_fraction = as.numeric(best_row$assigned_fraction),
  best_n_retained_clusters = as.integer(best_row$n_retained_clusters),
  best_min_retained_cluster_size = as.integer(best_row$min_retained_cluster_size),
  best_max_retained_cluster_size = as.integer(best_row$max_retained_cluster_size),
  best_raw_min_cluster_size = as.integer(best_row$raw_min_cluster_size),
  best_raw_max_cluster_size = as.integer(best_row$raw_max_cluster_size),
  best_n_removed_small_clusters = as.integer(best_row$n_removed_small_clusters),
  best_n_cells_moved_from_small_clusters_to_unassigned = as.integer(best_row$n_cells_moved_from_small_clusters_to_unassigned),
  coverage_threshold = coverage_threshold,
  coverage_rule = "depth > coverage_threshold",
  min_shared_variants = min_shared_variants,
  na_distance_action = na_distance_action,
  require_alt_signal_cell = require_alt_signal_cell,
  min_alt_for_cell = min_alt_for_cell,
  auto_k_min = auto_k_min,
  auto_k_max = auto_k_max,
  auto_k_min_cluster_size = auto_k_min_cluster_size,
  small_cluster_rule = "clusters smaller than auto_k_min_cluster_size moved to unassigned",
  best_k_rule = best_k_selection_rule,
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

fwrite(qc, paste0(out_prefix, ".publication_mito_distance.qc.tsv"), sep = "\t")

message("============================================================")
message("DONE")
message("best_k: ", best_k)
message("best_k_selection_mode: ", best_k_selection_mode)
message("best_k_selection_rule: ", best_k_selection_rule)
message("assigned-only clusters_out: ", clusters_out)
message("full assignments: ", paste0(out_prefix, ".hclust.assignments.k", best_k, ".tsv"))
message("auto_k_selection: ", paste0(out_prefix, ".auto_k_selection.tsv"))
message("best_k file: ", paste0(out_prefix, ".best_k.txt"))
message("QC: ", paste0(out_prefix, ".publication_mito_distance.qc.tsv"))
message("============================================================")
