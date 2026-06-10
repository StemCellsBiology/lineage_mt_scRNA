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

n_clust_cells <- nrow(D_clust)
candidate_k_max <- min(auto_k_max, floor(n_clust_cells / 2))

if (candidate_k_max < auto_k_min) {
  candidate_k <- 2
  if (candidate_k > n_clust_cells) stop("Not enough cells for k=2 after filtering.")
  warning("Candidate k range collapsed. Testing only k=2.")
} else {
  candidate_k <- seq(auto_k_min, candidate_k_max)
}

auto_results <- rbindlist(lapply(candidate_k, function(kk) {
  cl <- cutree(hc, k = kk)
  tab <- table(cl)
  min_size <- min(tab)
  max_size <- max(tab)
  sil <- cluster::silhouette(cl, dist_obj)
  sil_values <- sil[, 3]
  data.frame(
    k = kk,
    mean_silhouette = mean(sil_values, na.rm = TRUE),
    median_silhouette = median(sil_values, na.rm = TRUE),
    min_cluster_size = as.integer(min_size),
    max_cluster_size = as.integer(max_size),
    n_negative_silhouette_cells = as.integer(sum(sil_values < 0, na.rm = TRUE)),
    passes_min_cluster_size = min_size >= auto_k_min_cluster_size
  )
}))

auto_results$rank_silhouette <- rank(
  -auto_results$mean_silhouette,
  ties.method = "first",
  na.last = "keep"
)

fwrite(auto_results, paste0(out_prefix, ".auto_k_selection.tsv"), sep = "\t")

valid_auto <- auto_results[
  is.finite(mean_silhouette) &
    passes_min_cluster_size == TRUE
]

if (nrow(valid_auto) == 0) {
  warning("No k passed min cluster size filter. Selecting best k by silhouette without this filter.")
  valid_auto <- auto_results[is.finite(mean_silhouette)]
}
if (nrow(valid_auto) == 0) stop("No finite silhouette values available.")

best_row <- valid_auto[which.max(valid_auto$mean_silhouette)]
best_k <- as.integer(best_row$k)

writeLines(as.character(best_k), paste0(out_prefix, ".best_k.txt"))

clusters <- cutree(hc, k = best_k)

clusters_df <- data.frame(
  barcode = names(clusters),
  mt_cluster = paste0("C", as.integer(clusters)),
  stringsAsFactors = FALSE
)

clusters_out <- paste0(cluster_out_base, ".k", best_k, ".tsv")
fwrite(clusters_df, clusters_out, sep = "\t")

cluster_sizes <- as.data.frame(table(clusters_df$mt_cluster))
colnames(cluster_sizes) <- c("mt_cluster", "n_cells")
fwrite(cluster_sizes, paste0(out_prefix, ".publication_mito_distance.cluster_sizes.k", best_k, ".tsv"), sep = "\t")

heatmap_pdf <- paste0(out_prefix, ".publication_mito_distance.heatmap.k", best_k, ".pdf")

if (requireNamespace("pheatmap", quietly = TRUE)) {
  ann <- data.frame(mt_cluster = clusters_df$mt_cluster, row.names = clusters_df$barcode)
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

finite_D <- D[is.finite(D)]
finite_D_offdiag <- D[offdiag & is.finite(D)]

qc <- data.frame(
  input_file = filtered_long_file,
  out_prefix = out_prefix,
  cluster_out_base = cluster_out_base,
  clusters_out = clusters_out,
  best_k = best_k,
  best_mean_silhouette = as.numeric(best_row$mean_silhouette),
  best_median_silhouette = as.numeric(best_row$median_silhouette),
  best_min_cluster_size = as.integer(best_row$min_cluster_size),
  coverage_threshold = coverage_threshold,
  coverage_rule = "depth > coverage_threshold",
  min_shared_variants = min_shared_variants,
  na_distance_action = na_distance_action,
  require_alt_signal_cell = require_alt_signal_cell,
  min_alt_for_cell = min_alt_for_cell,
  auto_k_min = auto_k_min,
  auto_k_max = auto_k_max,
  auto_k_min_cluster_size = auto_k_min_cluster_size,
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
message("clusters_out: ", clusters_out)
message("auto_k_selection: ", paste0(out_prefix, ".auto_k_selection.tsv"))
message("best_k file: ", paste0(out_prefix, ".best_k.txt"))
message("QC: ", paste0(out_prefix, ".publication_mito_distance.qc.tsv"))
message("============================================================")
