#!/usr/bin/env Rscript

# ============================================================
# qc_mt_clusters_on_seurat_safe_terminal_v6_multik.R
#
# Multi-k version for STEP 15 of pipeline21.bash.
#
# Purpose:
#   Add mtDNA cluster assignments for multiple k values to one Seurat object
#   while loading the original Seurat RDS only once.
#
# Additional feature:
#   Add per-cell mtDNA variant information from:
#     *.filtered_variants_by_cell.long.tsv.gz
#
#   The information is added both to Seurat metadata and to the interactive
#   HTML hover text:
#
#     variant
#     depth
#     alt_count
#     vaf
#     alt_mean_bq
#
# Expected arguments:
#
# Rscript qc_mt_clusters_on_seurat_safe_terminal_v6_multik.R \
#   <seurat.rds> \
#   <mt_reads_per_cell.tsv> \
#   <filtered_variants_by_cell.long.tsv.gz> \
#   <step15_cluster_files.tsv> \
#   <out_prefix_base> \
#   <barcode_prefix_or_NONE> \
#   <mt_cluster_colname_prefix> \
#   <metadata_check_column_or_NONE> \
#   <expected_metadata_value_or_NONE>
#
# Example from pipeline21.bash:
#
# Rscript qc_mt_clusters_on_seurat_safe_terminal_v6_multik.R \
#   ~/KSzade/SCT_annotated_KSz.Rds \
#   sample.uniq_mt_reads_per_cell.tsv \
#   sample.filtered_variants_by_cell.long.tsv.gz \
#   sample.step15_cluster_files.10.txt \
#   sample.mt.multik \
#   sample_ \
#   sample.mt_cluster \
#   sample \
#   sample
#
# The cluster list file must have two tab-separated columns:
#
#   k    cluster_file
#
# Example:
#
#   20   /path/to/sample.clusters.10.k20.tsv
#   10   /path/to/sample.clusters.10.k10.tsv
#
# Main outputs:
#
#   <out_prefix_base>.seurat_with_mt_lineage_multik.rds
#   <out_prefix_base>.seurat_only_mt_assigned_cells_union_multik.rds
#   <out_prefix_base>.seurat_metadata_with_mt_lineage_multik.tsv
#   <out_prefix_base>.mt_variant_hover_per_cell.tsv
#   <out_prefix_base>.filtered_variants_by_cell.long.with_seurat_barcodes.tsv
#
# Per-k outputs:
#
#   <out_prefix_base>.k<K>.RNA_UMAP_harmony_by_<cluster_col>.pdf
#   <out_prefix_base>.k<K>.RNA_UMAP_harmony_<sample>_mt_clusters_highlighted.svg
#   <out_prefix_base>.k<K>.RNA_UMAP_harmony_<sample>_mt_clusters_highlighted.interactive.html
#   <out_prefix_base>.k<K>.<cluster_col>_vs_sample.tsv
#   <out_prefix_base>.k<K>.<cluster_col>_vs_patient.tsv
#   <out_prefix_base>.k<K>.<cluster_col>_vs_seurat_clusters.tsv
#   <out_prefix_base>.k<K>.<cluster_col>_vs_cluster_annotation.tsv
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

# ============================================================
# 1. ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  stop(
    "Usage:\n",
    "Rscript qc_mt_clusters_on_seurat_safe_terminal_v6_multik.R ",
    "<seurat.rds> ",
    "<mt_reads_per_cell.tsv> ",
    "<filtered_variants_by_cell.long.tsv.gz> ",
    "<step15_cluster_files.tsv> ",
    "<out_prefix_base> ",
    "<barcode_prefix_or_NONE> ",
    "<mt_cluster_colname_prefix> ",
    "<metadata_check_column_or_NONE> ",
    "<expected_metadata_value_or_NONE>\n"
  )
}

rds_file <- args[1]
mt_reads_file <- args[2]
mt_variants_long_file <- args[3]
cluster_files_list <- args[4]
out_prefix_base <- args[5]
barcode_prefix <- args[6]
mt_cluster_colname_prefix <- args[7]
metadata_check_column <- args[8]
expected_metadata_value <- args[9]
publication_coverage_threshold <- as.numeric(args[10])

if (is.na(publication_coverage_threshold)) {
  stop("publication_coverage_threshold must be numeric.")
}

message("Publication-like coverage threshold for hover: ",
        publication_coverage_threshold)
message("Hover coverage rule: depth > publication_coverage_threshold and vaf > 0")

message("============================================================")
message("qc_mt_clusters_on_seurat_safe_terminal_v6_multik.R")
message("============================================================")
message("Input Seurat RDS: ", rds_file)
message("Input mtDNA reads per cell: ", mt_reads_file)
message("Input mtDNA variants long table: ", mt_variants_long_file)
message("Input cluster list: ", cluster_files_list)
message("Output prefix base: ", out_prefix_base)
message("Barcode prefix: ", barcode_prefix)
message("mtDNA cluster metadata prefix: ", mt_cluster_colname_prefix)
message("Metadata check column: ", metadata_check_column)
message("Expected metadata value: ", expected_metadata_value)
message("============================================================")

# ============================================================
# 2. SMALL UTILITY FUNCTIONS
# ============================================================

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x
}

safe_num <- function(x, digits = 4) {
  x <- suppressWarnings(as.numeric(x))
  out <- ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
  out
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

read_required_table <- function(path, description) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(description, " does not exist or is empty: ", path)
  }
  fread(path, data.table = FALSE)
}

require_cols <- function(df, cols, description) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(
      description,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

make_variant_hover <- function(df, max_variants = 20) {

  if (nrow(df) == 0) {
    return("none")
  }

  n_total <- nrow(df)
  df <- df[seq_len(min(n_total, max_variants)), , drop = FALSE]

  lines <- paste0(
    safe_chr(df$variant),
    " | depth=", safe_chr(df$depth),
    " | alt_count=", safe_chr(df$alt_count),
    " | vaf=", safe_num(df$vaf, digits = 4),
    " | alt_mean_bq=", safe_num(df$alt_mean_bq, digits = 2)
  )

  out <- paste(lines, collapse = "<br>")

  if (n_total > max_variants) {
    out <- paste0(
      out,
      "<br>... ",
      n_total - max_variants,
      " additional variants not shown"
    )
  }

  out
}

add_or_replace_metadata_vector <- function(obj, values, colname) {
  if (colname %in% colnames(obj@meta.data)) {
    message("Metadata column already exists and will be overwritten: ", colname)
    obj@meta.data[[colname]] <- values
  } else {
    obj <- AddMetaData(
      object = obj,
      metadata = values,
      col.name = colname
    )
  }
  obj
}

save_vector <- function(x, file) {
  write.table(
    x,
    file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

make_table_against_metadata <- function(meta, cluster_col, group_col, out_file) {

  if (!cluster_col %in% colnames(meta)) {
    warning("Cluster column not found: ", cluster_col)
    return(invisible(NULL))
  }

  if (!group_col %in% colnames(meta)) {
    message("Metadata column not present, skipping table: ", group_col)
    return(invisible(NULL))
  }

  tab <- as.data.frame.matrix(
    table(meta[[cluster_col]], meta[[group_col]])
  )

  fwrite(
    data.table(mt_cluster = rownames(tab), tab),
    out_file,
    sep = "\t"
  )

  invisible(NULL)
}

# ============================================================
# 3. LOAD INPUT FILES ONLY ONCE
# ============================================================

message("Loading Seurat object...")
obj <- readRDS(rds_file)

message("Loading unique mtDNA reads table...")
mtreads <- read_required_table(mt_reads_file, "mtDNA reads file")

message("Loading filtered mtDNA variants-by-cell long table...")
mtlong <- read_required_table(mt_variants_long_file, "mtDNA variants long file")

message("Loading STEP 15 cluster file list...")
cluster_list <- read_required_table(cluster_files_list, "STEP 15 cluster list")

# ============================================================
# 4. CHECK INPUT COLUMNS
# ============================================================

require_cols(
  mtreads,
  c("barcode", "unique_mt_reads"),
  "mtDNA reads file"
)

require_cols(
  mtlong,
  c("barcode", "variant", "depth", "alt_count", "vaf", "alt_mean_bq"),
  "mtDNA variants long file"
)

if (ncol(cluster_list) < 2) {
  stop(
    "STEP 15 cluster list must contain at least two columns: k and cluster_file. ",
    "Problematic file: ", cluster_files_list
  )
}

cluster_list <- cluster_list[, 1:2, drop = FALSE]
colnames(cluster_list) <- c("k", "cluster_file")
cluster_list$k <- as.character(cluster_list$k)
cluster_list$cluster_file <- as.character(cluster_list$cluster_file)

cluster_list <- cluster_list[
  !is.na(cluster_list$k) &
    cluster_list$k != "" &
    !is.na(cluster_list$cluster_file) &
    cluster_list$cluster_file != "",
  ,
  drop = FALSE
]

if (nrow(cluster_list) == 0) {
  stop("No valid cluster files found in STEP 15 cluster list: ", cluster_files_list)
}

for (f in cluster_list$cluster_file) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop("Cluster file listed in STEP 15 cluster list does not exist or is empty: ", f)
  }
}

mtreads$barcode <- as.character(mtreads$barcode)
mtreads$unique_mt_reads <- as.numeric(mtreads$unique_mt_reads)

mtlong$barcode <- as.character(mtlong$barcode)
mtlong$variant <- as.character(mtlong$variant)
mtlong$depth <- as.numeric(mtlong$depth)
mtlong$alt_count <- as.numeric(mtlong$alt_count)
mtlong$vaf <- as.numeric(mtlong$vaf)
mtlong$alt_mean_bq <- as.numeric(mtlong$alt_mean_bq)

if (anyDuplicated(mtreads$barcode) > 0) {
  warning("Duplicated barcodes found in mtDNA reads file. Keeping first occurrence.")
  mtreads <- mtreads[!duplicated(mtreads$barcode), , drop = FALSE]
}

# ============================================================
# 5. ADD BARCODE PREFIX TO COMMON mtDNA TABLES
# ============================================================

if (barcode_prefix != "NONE") {
  message("Adding barcode prefix to mtDNA barcodes: ", barcode_prefix)
  mtreads$barcode <- paste0(barcode_prefix, mtreads$barcode)
  mtlong$barcode <- paste0(barcode_prefix, mtlong$barcode)
} else {
  message("No barcode prefix added.")
  warning(
    "barcode_prefix is NONE. This is only safe if mtDNA barcodes already ",
    "match Seurat cell names exactly. For multi-sample Seurat objects, ",
    "using a sample-specific prefix is strongly recommended."
  )
}

seurat_barcodes <- colnames(obj)

message("Cells in Seurat object: ", length(seurat_barcodes))
message("Cells in mtDNA reads file after prefix handling: ", nrow(mtreads))
message("Unique cells in mtDNA long table after prefix handling: ", length(unique(mtlong$barcode)))

matched_mtreads <- intersect(seurat_barcodes, mtreads$barcode)
matched_mtlong <- intersect(seurat_barcodes, unique(mtlong$barcode))

message("Matched cells with unique mtDNA reads: ", length(matched_mtreads))
message("Matched cells with mtDNA variants long table: ", length(matched_mtlong))

if (length(matched_mtreads) == 0) {
  warning(
    "No matching barcodes between Seurat object and mtDNA reads file. ",
    "The script will continue, but unique_mt_reads will be 0 for all cells."
  )
}

if (length(matched_mtlong) == 0) {
  warning(
    "No matching barcodes between Seurat object and mtDNA variants long table. ",
    "The script will continue, but variant hover information will be 'none'."
  )
}

# Save long table with Seurat-compatible barcode names.
fwrite(
  mtlong,
  paste0(out_prefix_base, ".filtered_variants_by_cell.long.with_seurat_barcodes.tsv"),
  sep = "\t"
)

# ============================================================
# 6. ADD unique_mt_reads TO SEURAT METADATA ONCE
# ============================================================

mtreads_vec <- setNames(
  mtreads$unique_mt_reads,
  mtreads$barcode
)

unique_mt_reads_to_add <- unname(mtreads_vec[seurat_barcodes])
unique_mt_reads_to_add[is.na(unique_mt_reads_to_add)] <- 0

obj <- add_or_replace_metadata_vector(
  obj = obj,
  values = unique_mt_reads_to_add,
  colname = "unique_mt_reads"
)

message("unique_mt_reads summary:")
print(summary(obj@meta.data[["unique_mt_reads"]]))

# ============================================================
# 7. PREPARE PER-CELL mtDNA VARIANT METADATA AND HOVER TEXT
# ============================================================

message("Preparing per-cell mtDNA variant metadata and hover text...")

# Keep only rows with true ALT signal for hover and summary.
# The long table may contain rows for variant-cell combinations where alt_count = 0.
mtlong_detected <- mtlong[
  !is.na(mtlong$depth) &
    mtlong$depth > publication_coverage_threshold &
    !is.na(mtlong$vaf) &
    mtlong$vaf > 0,
  ,
  drop = FALSE
]

message("Rows in mtDNA long table: ", nrow(mtlong))
message("Rows with detected ALT signal: ", nrow(mtlong_detected))
message("Cells with detected ALT signal: ", length(unique(mtlong_detected$barcode)))

if (nrow(mtlong_detected) > 0) {

  mtlong_detected <- mtlong_detected[
    order(
      mtlong_detected$barcode,
      -mtlong_detected$vaf,
      -mtlong_detected$alt_count,
      mtlong_detected$variant
    ),
    ,
    drop = FALSE
  ]

  split_by_cell <- split(mtlong_detected, mtlong_detected$barcode)

  mt_variant_hover <- do.call(
    rbind,
    lapply(
      split_by_cell,
      function(x) {
        data.frame(
          barcode = unique(x$barcode)[1],
          mt_n_variants_detected = length(unique(x$variant)),
          mt_variants_detected = paste(unique(x$variant), collapse = ";"),
          mt_variants_hover = make_variant_hover(x, max_variants = 20),
          mt_top_variant = x$variant[1],
          mt_top_variant_depth = x$depth[1],
          mt_top_variant_alt_count = x$alt_count[1],
          mt_top_variant_vaf = x$vaf[1],
          mt_top_variant_alt_mean_bq = x$alt_mean_bq[1],
          mt_sum_alt_count = sum(x$alt_count, na.rm = TRUE),
          mt_mean_vaf = mean(x$vaf, na.rm = TRUE),
          mt_max_vaf = max(x$vaf, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    )
  )

} else {

  mt_variant_hover <- data.frame(
    barcode = character(0),
    mt_n_variants_detected = integer(0),
    mt_variants_detected = character(0),
    mt_variants_hover = character(0),
    mt_top_variant = character(0),
    mt_top_variant_depth = numeric(0),
    mt_top_variant_alt_count = numeric(0),
    mt_top_variant_vaf = numeric(0),
    mt_top_variant_alt_mean_bq = numeric(0),
    mt_sum_alt_count = numeric(0),
    mt_mean_vaf = numeric(0),
    mt_max_vaf = numeric(0),
    stringsAsFactors = FALSE
  )
}

fwrite(
  mt_variant_hover,
  paste0(out_prefix_base, ".mt_variant_hover_per_cell.tsv"),
  sep = "\t"
)

# Build metadata data.frame aligned to all Seurat cells.
n_seurat_cells <- length(seurat_barcodes)

mt_variant_meta <- data.frame(
  mt_n_variants_detected = rep(0L, n_seurat_cells),
  mt_variants_detected = rep("", n_seurat_cells),
  mt_variants_hover = rep("none", n_seurat_cells),
  mt_top_variant = rep("none", n_seurat_cells),
  mt_top_variant_depth = rep(NA_real_, n_seurat_cells),
  mt_top_variant_alt_count = rep(NA_real_, n_seurat_cells),
  mt_top_variant_vaf = rep(NA_real_, n_seurat_cells),
  mt_top_variant_alt_mean_bq = rep(NA_real_, n_seurat_cells),
  mt_sum_alt_count = rep(0, n_seurat_cells),
  mt_mean_vaf = rep(NA_real_, n_seurat_cells),
  mt_max_vaf = rep(NA_real_, n_seurat_cells),
  stringsAsFactors = FALSE
)

rownames(mt_variant_meta) <- seurat_barcodes

if (nrow(mt_variant_hover) > 0) {

  rownames(mt_variant_hover) <- mt_variant_hover$barcode

  common_variant_cells <- intersect(seurat_barcodes, rownames(mt_variant_hover))
  cols_to_transfer <- setdiff(colnames(mt_variant_hover), "barcode")

  if (length(common_variant_cells) > 0) {
    mt_variant_meta[common_variant_cells, cols_to_transfer] <-
      mt_variant_hover[common_variant_cells, cols_to_transfer, drop = FALSE]
  }
}

obj <- AddMetaData(
  object = obj,
  metadata = mt_variant_meta
)

message("Added mtDNA variant summary metadata to Seurat object.")
message("Summary of mt_n_variants_detected:")
print(summary(obj@meta.data[["mt_n_variants_detected"]]))

message("Most frequent top variants:")
print(head(sort(table(obj@meta.data[["mt_top_variant"]]), decreasing = TRUE), 20))

# ============================================================
# 8. GLOBAL BARCODE QC FOR COMMON TABLES
# ============================================================

barcode_qc_global <- data.frame(
  n_cells_seurat = length(seurat_barcodes),
  n_cells_mt_reads_after_prefix = nrow(mtreads),
  n_cells_mtlong_unique_after_prefix = length(unique(mtlong$barcode)),
  n_matched_mt_reads = length(matched_mtreads),
  n_matched_mtlong = length(matched_mtlong),
  n_unmatched_seurat_vs_mtreads = length(setdiff(seurat_barcodes, mtreads$barcode)),
  n_unmatched_mtreads_vs_seurat = length(setdiff(mtreads$barcode, seurat_barcodes)),
  n_unmatched_seurat_vs_mtlong = length(setdiff(seurat_barcodes, unique(mtlong$barcode))),
  n_unmatched_mtlong_vs_seurat = length(setdiff(unique(mtlong$barcode), seurat_barcodes)),
  barcode_prefix = barcode_prefix,
  metadata_check_column = metadata_check_column,
  expected_metadata_value = expected_metadata_value,
  stringsAsFactors = FALSE
)

write.csv(
  barcode_qc_global,
  paste0(out_prefix_base, ".barcode_matching_qc_common_tables.csv"),
  row.names = FALSE
)

save_vector(
  setdiff(seurat_barcodes, mtreads$barcode),
  paste0(out_prefix_base, ".unmatched_seurat_barcodes_vs_mtreads.txt")
)

save_vector(
  setdiff(mtreads$barcode, seurat_barcodes),
  paste0(out_prefix_base, ".unmatched_mtread_barcodes_vs_seurat.txt")
)

save_vector(
  setdiff(seurat_barcodes, unique(mtlong$barcode)),
  paste0(out_prefix_base, ".unmatched_seurat_barcodes_vs_mtlong.txt")
)

save_vector(
  setdiff(unique(mtlong$barcode), seurat_barcodes),
  paste0(out_prefix_base, ".unmatched_mtlong_barcodes_vs_seurat.txt")
)

# ============================================================
# 9. PREPARE BASE UMAP EMBEDDING ONCE
# ============================================================

has_umap_harmony <- "umap_harmony" %in% names(obj@reductions)

if (!has_umap_harmony) {
  warning("umap_harmony reduction not found in Seurat object. UMAP_harmony plots will be skipped.")
} else {

  base_emb <- as.data.frame(
    Embeddings(obj, reduction = "umap_harmony")
  )

  base_emb$cell <- rownames(base_emb)

  xcol <- colnames(base_emb)[1]
  ycol <- colnames(base_emb)[2]

  message("UMAP_harmony coordinate columns used for SVG/HTML: ", xcol, ", ", ycol)

  if (metadata_check_column != "NONE") {
    if (!metadata_check_column %in% colnames(obj@meta.data)) {
      stop(
        "Metadata check column not found in Seurat metadata: ",
        metadata_check_column,
        ". Available columns are: ",
        paste(colnames(obj@meta.data), collapse = ", ")
      )
    }
    base_emb$highlight_group <- obj@meta.data[base_emb$cell, metadata_check_column]
  } else {
    base_emb$highlight_group <- NA
  }

  base_emb$unique_mt_reads <- obj@meta.data[base_emb$cell, "unique_mt_reads"]
  base_emb$unique_mt_reads[is.na(base_emb$unique_mt_reads)] <- 0

  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    base_emb$seurat_clusters <- obj@meta.data[base_emb$cell, "seurat_clusters"]
  } else {
    base_emb$seurat_clusters <- NA
  }

  if ("cluster_annotation" %in% colnames(obj@meta.data)) {
    base_emb$cluster_annotation <- obj@meta.data[base_emb$cell, "cluster_annotation"]
  } else {
    base_emb$cluster_annotation <- NA
  }

  # Add variant metadata for hover.
  variant_hover_cols <- c(
    "mt_n_variants_detected",
    "mt_variants_hover",
    "mt_top_variant",
    "mt_top_variant_depth",
    "mt_top_variant_alt_count",
    "mt_top_variant_vaf",
    "mt_top_variant_alt_mean_bq"
  )

  for (cc in variant_hover_cols) {
    base_emb[[cc]] <- obj@meta.data[base_emb$cell, cc]
  }

  base_emb$mt_n_variants_detected[is.na(base_emb$mt_n_variants_detected)] <- 0
  base_emb$mt_variants_hover[is.na(base_emb$mt_variants_hover)] <- "none"
  base_emb$mt_top_variant[is.na(base_emb$mt_top_variant)] <- "none"
}

# ============================================================
# 10. LOOP OVER k VALUES
# ============================================================

all_mt_cluster_cols <- character(0)

for (ii in seq_len(nrow(cluster_list))) {

  K <- as.character(cluster_list$k[ii])
  cluster_file <- as.character(cluster_list$cluster_file[ii])

  k_label <- paste0("k", K)
  out_prefix_k <- paste0(out_prefix_base, ".", k_label)
  mt_cluster_colname <- paste0(mt_cluster_colname_prefix, "_", k_label)

  message("============================================================")
  message("Processing STEP 15 for ", k_label)
  message("Cluster file: ", cluster_file)
  message("Output prefix for this k: ", out_prefix_k)
  message("Metadata column for this k: ", mt_cluster_colname)
  message("============================================================")

  clusters <- fread(cluster_file, data.table = FALSE)

  require_cols(
    clusters,
    c("barcode", "mt_cluster"),
    paste0("mtDNA cluster file for ", k_label)
  )

  clusters$barcode <- as.character(clusters$barcode)
  clusters$mt_cluster <- as.character(clusters$mt_cluster)

  if (anyDuplicated(clusters$barcode) > 0) {
    warning(
      "Duplicated barcodes found in mtDNA cluster file for ",
      k_label,
      ". Keeping first occurrence."
    )
    clusters <- clusters[!duplicated(clusters$barcode), , drop = FALSE]
  }

  if (barcode_prefix != "NONE") {
    clusters$barcode <- paste0(barcode_prefix, clusters$barcode)
  }

  matched_clusters <- intersect(seurat_barcodes, clusters$barcode)

  message("Cells in mtDNA cluster file after prefix handling: ", nrow(clusters))
  message("Matched cells with mtDNA clusters for ", k_label, ": ", length(matched_clusters))

  if (length(matched_clusters) == 0) {
    stop(
      "No matching barcodes between Seurat object and mtDNA cluster file for ",
      k_label,
      ". Check barcode_prefix and cluster file: ",
      cluster_file
    )
  }

  if (barcode_prefix != "NONE") {
    wrong_prefix_matches <- matched_clusters[!startsWith(matched_clusters, barcode_prefix)]

    if (length(wrong_prefix_matches) > 0) {
      stop(
        "Some matched mtDNA cluster barcodes for ",
        k_label,
        " do not start with the requested barcode_prefix. Check barcode handling."
      )
    }
  }

  # Per-k barcode QC.
  barcode_qc_k <- data.frame(
    k = K,
    cluster_file = cluster_file,
    n_cells_seurat = length(seurat_barcodes),
    n_cells_mt_clusters_after_prefix = nrow(clusters),
    n_matched_clusters = length(matched_clusters),
    n_unmatched_seurat_vs_clusters = length(setdiff(seurat_barcodes, clusters$barcode)),
    n_unmatched_clusters_vs_seurat = length(setdiff(clusters$barcode, seurat_barcodes)),
    barcode_prefix = barcode_prefix,
    metadata_check_column = metadata_check_column,
    expected_metadata_value = expected_metadata_value,
    stringsAsFactors = FALSE
  )

  write.csv(
    barcode_qc_k,
    paste0(out_prefix_k, ".barcode_matching_qc.csv"),
    row.names = FALSE
  )

  save_vector(
    setdiff(seurat_barcodes, clusters$barcode),
    paste0(out_prefix_k, ".unmatched_seurat_barcodes_vs_clusters.txt")
  )

  save_vector(
    setdiff(clusters$barcode, seurat_barcodes),
    paste0(out_prefix_k, ".unmatched_cluster_barcodes_vs_seurat.txt")
  )

  # Add cluster column to Seurat metadata.
  cluster_vec <- setNames(
    clusters$mt_cluster,
    clusters$barcode
  )

  mt_cluster_to_add <- unname(cluster_vec[seurat_barcodes])
  mt_cluster_to_add[is.na(mt_cluster_to_add)] <- "unassigned"

  obj <- add_or_replace_metadata_vector(
    obj = obj,
    values = mt_cluster_to_add,
    colname = mt_cluster_colname
  )

  all_mt_cluster_cols <- c(all_mt_cluster_cols, mt_cluster_colname)

  message("mtDNA cluster distribution for ", k_label, ":")
  print(table(obj@meta.data[[mt_cluster_colname]], useNA = "ifany"))

  # ============================================================
  # Per-k safety check
  # ============================================================

  assigned_cells <- colnames(obj)[obj@meta.data[[mt_cluster_colname]] != "unassigned"]

  message("Cells assigned to mtDNA clusters for ", k_label, ": ", length(assigned_cells))

  if (metadata_check_column != "NONE" && expected_metadata_value != "NONE") {

    if (!metadata_check_column %in% colnames(obj@meta.data)) {
      stop(
        "Metadata check column not found in Seurat metadata: ",
        metadata_check_column,
        ". Available columns are: ",
        paste(colnames(obj@meta.data), collapse = ", ")
      )
    }

    message(
      "Distribution of ",
      metadata_check_column,
      " among mtDNA-assigned cells for ",
      k_label,
      ":"
    )

    meta_tab <- table(
      obj@meta.data[assigned_cells, metadata_check_column],
      useNA = "ifany"
    )

    print(meta_tab)

    meta_check_df <- data.frame(
      k = K,
      metadata_column = metadata_check_column,
      value = names(meta_tab),
      n_cells = as.integer(meta_tab),
      stringsAsFactors = FALSE
    )

    fwrite(
      meta_check_df,
      paste0(out_prefix_k, ".", metadata_check_column, "_distribution_among_mt_assigned_cells.tsv"),
      sep = "\t"
    )

    assigned_values <- unique(
      as.character(obj@meta.data[assigned_cells, metadata_check_column])
    )
    assigned_values <- assigned_values[!is.na(assigned_values)]

    unexpected_values <- setdiff(assigned_values, expected_metadata_value)

    if (length(unexpected_values) > 0) {
      stop(
        "Safety check failed for ",
        k_label,
        ": mtDNA-assigned cells are present in unexpected values of metadata column '",
        metadata_check_column,
        "': ",
        paste(unexpected_values, collapse = ", "),
        ". Expected only: ",
        expected_metadata_value,
        ". Check barcode_prefix and input files."
      )
    }

    message(
      "Safety check passed for ",
      k_label,
      ": all mtDNA-assigned cells have ",
      metadata_check_column,
      " = ",
      expected_metadata_value
    )

  } else {
    message(
      "Metadata safety check skipped for ",
      k_label,
      " because metadata_check_column or expected_metadata_value is NONE."
    )
  }

  # ============================================================
  # Per-k UMAP/PDF/SVG/HTML outputs
  # ============================================================

  if (!has_umap_harmony) {
    message("Skipping UMAP_harmony outputs for ", k_label, " because reduction is missing.")
  } else {

    message("Generating RNA UMAP_harmony colored by mtDNA cluster for ", k_label, "...")

    p1 <- DimPlot(
      obj,
      reduction = "umap_harmony",
      group.by = mt_cluster_colname,
      pt.size = 0.4
    ) +
      ggtitle(paste0("RNA UMAP_harmony colored by ", mt_cluster_colname))

    ggsave(
      paste0(out_prefix_k, ".RNA_UMAP_harmony_by_", mt_cluster_colname, ".pdf"),
      p1,
      width = 8,
      height = 6
    )

    # QC variable UMAPs for this k.
    qc_vars <- intersect(
      c(
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt",
        "unique_mt_reads",
        "mt_n_variants_detected",
        "mt_top_variant_vaf",
        "mt_top_variant_alt_count",
        "orig.ident",
        "sample",
        "patient",
        "sample_type",
        "seurat_clusters",
        "cluster_annotation",
        mt_cluster_colname
      ),
      colnames(obj@meta.data)
    )

    for (v in qc_vars) {

      message("Plotting RNA UMAP_harmony by: ", v, " for ", k_label)

      p <- if (is.numeric(obj@meta.data[[v]])) {
        FeaturePlot(
          obj,
          features = v,
          reduction = "umap_harmony",
          pt.size = 0.3
        ) +
          ggtitle(paste0("RNA UMAP_harmony colored by ", v))
      } else {
        DimPlot(
          obj,
          reduction = "umap_harmony",
          group.by = v,
          pt.size = 0.3
        ) +
          ggtitle(paste0("RNA UMAP_harmony colored by ", v))
      }

      ggsave(
        paste0(out_prefix_k, ".RNA_UMAP_harmony_by_", safe_filename(v), ".pdf"),
        p,
        width = 8,
        height = 6
      )
    }

    if (metadata_check_column != "NONE" && expected_metadata_value != "NONE") {

      emb <- base_emb

      emb$mt_cluster <- obj@meta.data[emb$cell, mt_cluster_colname]
      emb$mt_cluster[is.na(emb$mt_cluster)] <- "unassigned"

      emb$hover_text <- paste0(
        "cell: ", safe_chr(emb$cell),
        "<br>", metadata_check_column, ": ", safe_chr(emb$highlight_group),
        "<br>mt_cluster_column: ", mt_cluster_colname,
        "<br>mt_cluster: ", safe_chr(emb$mt_cluster),
        "<br>unique_mt_reads: ", safe_chr(emb$unique_mt_reads),
        "<br>seurat_clusters: ", safe_chr(emb$seurat_clusters),
        "<br>cluster_annotation: ", safe_chr(emb$cluster_annotation),
        "<br>mt_n_variants_detected: ", safe_chr(emb$mt_n_variants_detected),
        "<br>mt_top_variant: ", safe_chr(emb$mt_top_variant),
        "<br>mt_top_variant_depth: ", safe_chr(emb$mt_top_variant_depth),
        "<br>mt_top_variant_alt_count: ", safe_chr(emb$mt_top_variant_alt_count),
        "<br>mt_top_variant_vaf: ", safe_num(emb$mt_top_variant_vaf, digits = 4),
        "<br>mt_top_variant_alt_mean_bq: ", safe_num(emb$mt_top_variant_alt_mean_bq, digits = 2),
        "<br><br>mtDNA variants:",
        "<br>", safe_chr(emb$mt_variants_hover)
      )

      df_background <- emb[
        emb$highlight_group != expected_metadata_value | is.na(emb$highlight_group),
        ,
        drop = FALSE
      ]

      df_sample <- emb[
        emb$highlight_group == expected_metadata_value & !is.na(emb$highlight_group),
        ,
        drop = FALSE
      ]

      df_sample_unassigned <- df_sample[
        df_sample$mt_cluster == "unassigned",
        ,
        drop = FALSE
      ]

      df_sample_assigned <- df_sample[
        df_sample$mt_cluster != "unassigned",
        ,
        drop = FALSE
      ]

      message("Cells in background for ", k_label, ": ", nrow(df_background))
      message("Cells in highlighted group for ", k_label, ": ", nrow(df_sample))
      message("Highlighted cells unassigned for ", k_label, ": ", nrow(df_sample_unassigned))
      message("Highlighted cells assigned to mtDNA clusters for ", k_label, ": ", nrow(df_sample_assigned))

      p_highlight <- ggplot() +
        geom_point(
          data = df_background,
          aes(x = .data[[xcol]], y = .data[[ycol]]),
          color = "grey85",
          alpha = 0.05,
          size = 0.05
        ) +
        geom_point(
          data = df_sample_unassigned,
          aes(x = .data[[xcol]], y = .data[[ycol]]),
          color = "grey65",
          alpha = 0.30,
          size = 0.30
        ) +
        geom_point(
          data = df_sample_assigned,
          aes(x = .data[[xcol]], y = .data[[ycol]], color = mt_cluster),
          alpha = 0.95,
          size = 0.90
        ) +
        theme_classic() +
        labs(
          title = paste0(
            "UMAP Harmony: ",
            expected_metadata_value,
            " mtDNA clusters highlighted, ",
            k_label
          ),
          x = "UMAP Harmony 1",
          y = "UMAP Harmony 2",
          color = paste0(expected_metadata_value, " mtDNA cluster")
        )

      svg_file <- paste0(
        out_prefix_k,
        ".RNA_UMAP_harmony_",
        safe_filename(expected_metadata_value),
        "_mt_clusters_highlighted.svg"
      )

      if (requireNamespace("svglite", quietly = TRUE)) {

        ggsave(
          svg_file,
          p_highlight,
          width = 8,
          height = 6,
          device = svglite::svglite
        )

        message("Saved scalable SVG UMAP_harmony highlight plot:")
        message(svg_file)

      } else {

        warning(
          "Package 'svglite' is not installed. ",
          "Skipping SVG export for highlighted UMAP_harmony. ",
          "Install it in R with: install.packages('svglite')"
        )
      }

      html_file <- paste0(
        out_prefix_k,
        ".RNA_UMAP_harmony_",
        safe_filename(expected_metadata_value),
        "_mt_clusters_highlighted.interactive.html"
      )

      if (
        requireNamespace("plotly", quietly = TRUE) &&
        requireNamespace("htmlwidgets", quietly = TRUE)
      ) {

        message("Generating interactive HTML UMAP_harmony highlight plot for ", k_label, "...")

        p_interactive <- plotly::plot_ly()

        if (nrow(df_background) > 0) {
          p_interactive <- plotly::add_markers(
            p_interactive,
            x = df_background[[xcol]],
            y = df_background[[ycol]],
            text = df_background$hover_text,
            hoverinfo = "text",
            name = "other cells",
            marker = list(
              color = "rgba(180,180,180,0.25)",
              size = 3
            )
          )
        }

        if (nrow(df_sample_unassigned) > 0) {
          p_interactive <- plotly::add_markers(
            p_interactive,
            x = df_sample_unassigned[[xcol]],
            y = df_sample_unassigned[[ycol]],
            text = df_sample_unassigned$hover_text,
            hoverinfo = "text",
            name = paste0(expected_metadata_value, " unassigned"),
            marker = list(
              color = "rgba(120,120,120,0.65)",
              size = 4
            )
          )
        }

        if (nrow(df_sample_assigned) > 0) {
          p_interactive <- plotly::add_markers(
            p_interactive,
            data = df_sample_assigned,
            x = df_sample_assigned[[xcol]],
            y = df_sample_assigned[[ycol]],
            text = df_sample_assigned$hover_text,
            hoverinfo = "text",
            color = df_sample_assigned$mt_cluster,
            name = "mtDNA clusters",
            marker = list(
              size = 5,
              opacity = 0.95
            )
          )
        }

        p_interactive <- plotly::layout(
          p_interactive,
          title = paste0(
            "UMAP Harmony: ",
            expected_metadata_value,
            " mtDNA clusters highlighted, ",
            k_label
          ),
          xaxis = list(title = "UMAP Harmony 1"),
          yaxis = list(title = "UMAP Harmony 2"),
          legend = list(title = list(text = "Group / mtDNA cluster"))
        )

        htmlwidgets::saveWidget(
          widget = p_interactive,
          file = html_file,
          selfcontained = TRUE
        )

        message("Saved interactive HTML UMAP_harmony highlight plot:")
        message(html_file)

      } else {

        warning(
          "Packages 'plotly' and/or 'htmlwidgets' are not installed. ",
          "Skipping interactive HTML export. ",
          "Install them in R with: install.packages(c('plotly', 'htmlwidgets'))"
        )
      }

    } else {

      message(
        "Skipping highlighted SVG/HTML UMAP_harmony for ",
        k_label,
        " because metadata_check_column or expected_metadata_value is NONE."
      )
    }
  }

  # ============================================================
  # Per-k metadata tables
  # ============================================================

  meta_now <- obj@meta.data
  meta_now$barcode <- rownames(meta_now)

  make_table_against_metadata(
    meta_now,
    mt_cluster_colname,
    "sample",
    paste0(out_prefix_k, ".", mt_cluster_colname, "_vs_sample.tsv")
  )

  make_table_against_metadata(
    meta_now,
    mt_cluster_colname,
    "patient",
    paste0(out_prefix_k, ".", mt_cluster_colname, "_vs_patient.tsv")
  )

  make_table_against_metadata(
    meta_now,
    mt_cluster_colname,
    "orig.ident",
    paste0(out_prefix_k, ".", mt_cluster_colname, "_vs_orig_ident.tsv")
  )

  make_table_against_metadata(
    meta_now,
    mt_cluster_colname,
    "seurat_clusters",
    paste0(out_prefix_k, ".", mt_cluster_colname, "_vs_seurat_clusters.tsv")
  )

  make_table_against_metadata(
    meta_now,
    mt_cluster_colname,
    "cluster_annotation",
    paste0(out_prefix_k, ".", mt_cluster_colname, "_vs_cluster_annotation.tsv")
  )

  # Boxplots by mtDNA cluster for this k.
  qc_numeric <- intersect(
    c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      "nCount_SCT",
      "nFeature_SCT",
      "mapping.score",
      "unique_mt_reads",
      "mt_n_variants_detected",
      "mt_top_variant_vaf",
      "mt_top_variant_alt_count"
    ),
    colnames(meta_now)
  )

  for (v in qc_numeric) {

    if (!is.numeric(meta_now[[v]])) {
      next
    }

    message("Plotting boxplot: ", v, " by ", mt_cluster_colname, " for ", k_label)

    p <- ggplot(
      meta_now,
      aes(x = .data[[mt_cluster_colname]], y = .data[[v]])
    ) +
      geom_boxplot(outlier.size = 0.4) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(
        x = mt_cluster_colname,
        y = v,
        title = paste0(v, " by ", mt_cluster_colname)
      )

    ggsave(
      paste0(out_prefix_k, ".boxplot_", safe_filename(v), "_by_", mt_cluster_colname, ".pdf"),
      p,
      width = 10,
      height = 5
    )
  }
}

# ============================================================
# 11. SAVE FINAL MULTI-k OUTPUTS
# ============================================================

all_mt_cluster_cols <- unique(all_mt_cluster_cols)

message("Saving final multi-k Seurat object...")

out_rds <- paste0(out_prefix_base, ".seurat_with_mt_lineage_multik.rds")
saveRDS(obj, out_rds)

message("Saved Seurat object with multi-k mtDNA lineage information:")
message(out_rds)

meta <- obj@meta.data
meta$barcode <- rownames(meta)

fwrite(
  meta,
  paste0(out_prefix_base, ".seurat_metadata_with_mt_lineage_multik.tsv"),
  sep = "\t"
)

message("Saved full multi-k metadata table:")
message(paste0(out_prefix_base, ".seurat_metadata_with_mt_lineage_multik.tsv"))

# Save union of cells assigned to any mtDNA cluster across all k.
if (length(all_mt_cluster_cols) > 0) {

  assigned_union <- rep(FALSE, nrow(obj@meta.data))
  names(assigned_union) <- rownames(obj@meta.data)

  for (cc in all_mt_cluster_cols) {
    assigned_union <- assigned_union | (!is.na(obj@meta.data[[cc]]) & obj@meta.data[[cc]] != "unassigned")
  }

  assigned_union_cells <- names(assigned_union)[assigned_union]

  message("Cells assigned to mtDNA clusters in at least one k: ", length(assigned_union_cells))

  if (length(assigned_union_cells) > 0) {

    obj_mt_assigned_union <- subset(
      obj,
      cells = assigned_union_cells
    )

    out_rds_assigned <- paste0(out_prefix_base, ".seurat_only_mt_assigned_cells_union_multik.rds")

    saveRDS(obj_mt_assigned_union, out_rds_assigned)

    message("Saved Seurat object containing union of cells assigned to mtDNA clusters:")
    message(out_rds_assigned)

  } else {
    warning("No cells assigned to mtDNA clusters in any k. Skipping assigned-union RDS.")
  }
}

# Summary of generated cluster columns.
cluster_col_summary <- data.frame(
  mt_cluster_column = all_mt_cluster_cols,
  stringsAsFactors = FALSE
)

if (length(all_mt_cluster_cols) > 0) {
  cluster_col_summary$n_assigned <- vapply(
    all_mt_cluster_cols,
    function(cc) {
      sum(!is.na(obj@meta.data[[cc]]) & obj@meta.data[[cc]] != "unassigned")
    },
    numeric(1)
  )
}

fwrite(
  cluster_col_summary,
  paste0(out_prefix_base, ".mt_cluster_columns_summary.tsv"),
  sep = "\t"
)

# ============================================================
# 12. FINAL SUMMARY
# ============================================================

message("============================================================")
message("DONE")
message("Main outputs:")
message("  ", out_rds)
message("  ", paste0(out_prefix_base, ".seurat_metadata_with_mt_lineage_multik.tsv"))
message("  ", paste0(out_prefix_base, ".mt_variant_hover_per_cell.tsv"))
message("  ", paste0(out_prefix_base, ".filtered_variants_by_cell.long.with_seurat_barcodes.tsv"))
message("  ", paste0(out_prefix_base, ".barcode_matching_qc_common_tables.csv"))
message("  ", paste0(out_prefix_base, ".mt_cluster_columns_summary.tsv"))

if (length(all_mt_cluster_cols) > 0) {
  message("Added mtDNA cluster metadata columns:")
  for (cc in all_mt_cluster_cols) {
    message("  ", cc)
  }
}

message("============================================================")
