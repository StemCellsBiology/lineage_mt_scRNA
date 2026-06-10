# Materials and Methods

## Single-cell mitochondrial lineage inference

The pipeline is based on the methodology from the paper Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics. https://doi.org/10.1016/j.cell.2019.01.022.
(https://www.sciencedirect.com/science/article/pii/S0092867419300558).

Mitochondrial lineage structure was inferred from sample-specific single-cell BAM files using a custom workflow implemented in Bash, Python and R. The analysis was performed independently for each sample. The workflow extracted uniquely mapped mitochondrial reads, quantified mitochondrial read depth per cell, identified and filtered candidate mitochondrial variants, retained variant information in a long-format per-cell table, computed a coverage-aware pairwise mitochondrial distance between cells, selected the optimal number of mitochondrial clusters by silhouette analysis, and integrated the resulting mtDNA cluster assignments with an annotated Seurat object.

Final lineage clustering was performed directly from the filtered long-format variant table, which contained barcode, variant identity, mitochondrial position, reference-like allele, alternative allele, depth, ALT read count, variant allele fraction (VAF), and ALT mean base quality. A cell-by-variant heteroplasmy matrix was not used for final clustering. This design preserved per-cell/per-variant coverage information during pairwise distance calculation and allowed distances to be computed only over variants sufficiently covered in both cells.

The main workflow used the mitochondrial chromosome identifier `chrM`, a minimum mapping quality threshold of 30, and a minimum per-cell mitochondrial read threshold of 200 for inclusion in allele counting. Publication-like distance calculation used a strict coverage rule of `depth > 9`, equivalent to requiring at least 10 reads, and required at least one shared covered variant between two cells. Cells used for mitochondrial distance calculation were required to carry detectable ALT signal, defined as `depth > 9`, `ALT count >= 2`, and `VAF > 0` for at least one variant. Candidate numbers of mitochondrial clusters were evaluated from `k = 2` up to the smaller of 30 and half the number of cells available for clustering. The final `best_k` was selected by mean silhouette score after applying a minimum cluster-size criterion of three cells where possible.

---

## Barcode selection and sample-level preprocessing

For each sample, the pipeline located the corresponding BAM file and associated Cell Ranger barcode file. Cell Ranger barcodes were converted to plain text when necessary and intersected with sample-specific barcodes extracted from the annotated Seurat metadata table. The Seurat barcode table was required to contain the columns `sample` and `barcode`. Only barcodes present both in the Cell Ranger barcode list and in the Seurat-derived sample barcode table were retained for mitochondrial analysis. Samples with no overlapping barcodes were excluded, as this indicated inconsistent sample naming, barcode suffixes, or sample assignment.

The final intersected barcode list defined the sample-specific cell set used during read counting and downstream mtDNA variant detection.

---

## STEP 1. Extraction of uniquely mapped mitochondrial reads

**Script:** `BAMs/scripts/filter_unique_mt_reads.py`

Uniquely mapped mitochondrial reads were extracted from each sample-specific BAM file using `pysam`. The script retained only reads satisfying all of the following criteria: the read was mapped, was not secondary, was not supplementary, mapped to `chrM`, had mapping quality ≥ 30, contained an `NH` tag, and had `NH = 1`. The `NH = 1` requirement was used as a strict aligner-derived criterion for unique mapping.

The output was a mitochondrial-only BAM file containing primary, uniquely mapped mitochondrial reads. This BAM file was indexed with `samtools` and used for all subsequent mitochondrial read counting and allele-counting steps.

---

## STEP 2. Quantification of mitochondrial reads per cell

**Script:** `BAMs/scripts/count_unique_mt_reads_per_cell.py`

Per-cell mitochondrial read depth was calculated from the uniquely mapped mitochondrial BAM file. The script iterated through the mitochondrial reads and counted only primary, mapped reads carrying a valid `CB` cell-barcode tag. Reads were counted only if the `CB` value matched a barcode from the sample-specific intersected barcode list.

The output was a tab-separated table containing one row per valid barcode and two columns: `barcode` and `unique_mt_reads`. This table was used both for mitochondrial coverage-based cell selection and for later addition of mitochondrial read-depth metadata to the Seurat object.

---

## STEP 3. Selection of cells for mitochondrial allele counting

**Script:** `BAMs/scripts/select_barcodes_by_mt_reads.py`

Cells were selected for allele counting based on the number of uniquely mapped mitochondrial reads. Barcodes with `unique_mt_reads >= 200` were retained. This threshold was used to reduce noise from cells with insufficient mitochondrial coverage while maintaining a defined set of cells for pileup-based variant discovery.

The output was a one-column barcode file containing cells that passed the mitochondrial read-depth threshold.

---

## STEP 4. Per-cell mitochondrial allele counting

**Script:** `BAMs/scripts/per_cell_per_allele_counts1.py`

Per-cell, per-position mitochondrial allele counts were generated from the uniquely mapped mitochondrial BAM file using the barcode list selected in STEP 3. The script performed a pileup across mitochondrial positions and counted A, C, G, and T separately for each cell barcode. For each cell-position combination with non-zero depth, the script recorded nucleotide counts, total depth, and mean base quality for each observed nucleotide.

No base-quality filtering was applied during this step. Pileup was performed with minimum base quality set to zero, allowing the empirical base-quality distribution to be modelled downstream rather than imposing an arbitrary threshold at the counting stage.

The resulting table contained the columns `barcode`, `pos`, `A`, `C`, `G`, `T`, `depth`, `A_mean_bq`, `C_mean_bq`, `G_mean_bq`, and `T_mean_bq`. The script also generated a global base-quality distribution plot from all observed mitochondrial bases passing the read-level filters.

---

## STEP 5. Empirical modelling of base-quality distributions

**Script:** `BAMs/scripts/gaussian_plot.py`

Base-quality structure was estimated from the per-cell allele-count table. The script aggregated counts across cells for each mitochondrial position and nucleotide, producing per-position/per-allele total counts and weighted mean base qualities. A three-component Gaussian mixture model was then fitted to the distribution of per-position/per-allele mean base qualities.

The Gaussian components were ordered by their mean base quality, and the component with the highest mean was treated as the high-confidence component. For every position-allele observation, the posterior probability of belonging to this high-confidence component was calculated. Observations with posterior probability greater than 0.99 were marked as high confidence.

The script generated a fitted density plot and a QC table containing, for each position-allele observation, the mitochondrial position, allele, total read count, mean base quality, posterior probability of high-confidence assignment, and pass/fail status for the 99% high-confidence criterion.

---

## STEP 6. Selection of the operational base-quality threshold

**Implementation:** Bash/awk within the main pipeline

The operational base-quality threshold was selected from the QC table generated in STEP 5. The pipeline identified all position-allele observations marked as high confidence and selected the minimum mean base-quality value among them. This value was converted to an integer and stored as the sample-specific threshold used in downstream filtering.

The selected threshold was written to a sample-specific file and reused when restarting the workflow from later steps, ensuring that base-quality filtering remained reproducible across partial reruns.

---

## STEP 7. Position-level filtering of low-quality alleles

**Script:** `BAMs/scripts/filter_per_cell_allele_counts_by_bq_plus.py`

Low-quality mitochondrial positions were filtered using the empirical base-quality threshold selected in STEP 6. The script first aggregated the per-cell allele-count table across cells for each position-allele pair. For each allele at each position, it calculated total allele count, total position depth, allele fraction, and weighted mean base quality.

A position-allele was considered problematic when its mean base quality was below the empirical threshold and its allele fraction was at least 0.002. Alleles below the base-quality threshold but with allele fraction < 0.002 were treated as low-impact observations and did not cause removal of the position. Filtering was applied conservatively at the position level: if any allele at a position was classified as problematic, all per-cell rows for that position were removed.

This step produced a BQ-filtered per-cell allele-count table and QC outputs listing retained positions, removed positions, and position-allele filtering statistics.

---

## STEP 8. Identification of candidate mitochondrial variants

**Script:** `BAMs/scripts/call_candidate_alleles.py`

Candidate mitochondrial variants were identified from the BQ-filtered per-cell allele-count table. For each mitochondrial position, counts of A, C, G, and T were pooled across all cells. The most abundant nucleotide was treated as the reference-like allele at that position. Each remaining nucleotide with non-zero total count was considered a candidate alternative allele.

For each candidate variant, the script calculated total depth, total ALT count, population-level heteroplasmy, number of ALT-positive cells, and mean base quality of ALT-supporting observations. Candidate variants were labelled using the format `position_refLike>alt`.

In parallel, the script generated a long-format per-cell variant table containing `barcode`, `variant`, `pos`, `ref_like`, `alt`, `depth`, `alt_count`, `vaf`, and `alt_mean_bq`. Here, VAF was defined as `alt_count / depth`. This long-format table preserved per-cell/per-variant coverage and allele-fraction information for downstream filtering and clustering.

---

## STEP 9. Filtering of candidate mitochondrial variants

**Script:** `BAMs/scripts/filter_variants_no_bulk_after_bq_plus.py`

Candidate mitochondrial variants were filtered using population-level and cell-level evidence thresholds. Because base-quality filtering had already been applied upstream, this step did not perform additional base-quality filtering and did not require a bulk mitochondrial reference profile.

In the main workflow, candidate variants were retained if they satisfied all of the following criteria: population heteroplasmy ≥ 0.025, at least one ALT-positive cell, total ALT count ≥ 1, and total depth ≥ 1. Variants failing any of these criteria were removed. The script generated a filtered candidate variant table, a filtered long-format per-cell variant table, and a QC table recording pass/fail status for each filtering criterion.

The filtered long-format table was the principal input for mitochondrial distance calculation and lineage clustering.

---

## STEP 10. Not used

No active analytical operation was implemented as STEP 10 in the final workflow.

---

## STEP 12. Publication-like mitochondrial distance calculation and clustering

**Script:** `BAMs/scripts_new/plot_heatmap_and_cluster_publication_mito_distance_auto_k.R`

Mitochondrial lineage clustering was performed directly from the filtered long-format per-cell variant table produced in STEP 9. The script first cleaned the table and required the columns `barcode`, `variant`, `depth`, `alt_count`, and `vaf`. The table was aggregated to obtain barcode-by-variant VAF and coverage matrices while preserving depth information for each cell-variant observation.

Cells were restricted to those carrying detectable mitochondrial ALT signal when the ALT-signal filter was enabled. A cell was considered to carry ALT signal if it had at least one variant observation with `depth > 9`, `alt_count >= 2`, and `vaf > 0`.

Pairwise mitochondrial distances were computed only over variants sufficiently covered in both cells. For a given pair of cells, a variant contributed to the distance only if:

```text
coverage_i > 9 and coverage_j > 9
```

For each pair of cells, the mitochondrial distance was calculated as:

```text
D_ij = mean sqrt(abs(VAF_i - VAF_j))
```

where the mean was taken across variants satisfying the shared coverage requirement. At least one shared covered variant was required. The corresponding mitochondrial relatedness matrix was defined as:

```text
K_mito = 1 - D
```

Missing pairwise distances were replaced with 1 before clustering, according to the default pipeline setting. The distance matrix was symmetrized, its diagonal was set to zero, and average-linkage hierarchical clustering was performed using:

```text
hclust(D, method = "average")
```

The optimal number of mitochondrial clusters was selected automatically. Candidate values of `k` were evaluated from 2 to the smaller of 30 and half the number of cells available for clustering. For each candidate `k`, cluster assignments were generated by cutting the hierarchical tree, and silhouette widths were calculated from the mitochondrial distance matrix. Candidate solutions were summarized by mean silhouette, median silhouette, minimum cluster size, maximum cluster size, and number of cells with negative silhouette values. Candidate values with minimum cluster size below three cells were excluded where possible. The selected `best_k` was the valid candidate with the highest mean silhouette score. If no candidate passed the minimum cluster-size criterion, the `k` with the highest mean silhouette score was selected without that filter.

The script wrote the selected `best_k`, final mitochondrial cluster assignments, pairwise distance matrix, relatedness matrix, shared-variant matrix, hierarchical clustering object, auto-k selection table, cluster-size table, heatmap, and clustering QC table.

---

## STEP 13. Disabled multi-k comparison

STEP 13 was disabled in the final workflow. Earlier versions compared multiple manually selected values of `k`; in the final version, this was replaced by automatic `best_k` selection in STEP 12.

---

## STEP 14. Disabled mitochondrial tSNE

STEP 14 was disabled in the final workflow. Mitochondrial tSNE was not used for final lineage assignment. Final mtDNA lineage structure was inferred from the coverage-aware mitochondrial distance matrix and average-linkage hierarchical clustering.

---

## STEP 15. Integration of mitochondrial clusters with the Seurat object

**Script:** `BAMs/scripts_new/qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R`

Final mitochondrial cluster assignments were integrated with the annotated Seurat object. Although the script was designed to support multiple values of `k`, the final workflow passed only the automatically selected `best_k` from STEP 12.

The script loaded the Seurat RDS object, the per-cell mitochondrial read-count table, the filtered long-format per-cell mitochondrial variant table, and the cluster assignment file. Sample-specific barcode prefixes were added to mitochondrial barcodes to match Seurat cell names. Barcode matching was evaluated between the Seurat object and each mitochondrial table, including the read-count table, variant table, and cluster-assignment table. Matched and unmatched barcode counts were saved as QC outputs.

The number of uniquely mapped mitochondrial reads was added to Seurat metadata as `unique_mt_reads`. Per-cell mitochondrial variant summaries were then generated from the filtered long-format table. Variant observations used for these summaries required `depth > 9` and `vaf > 0`. For each cell, the script recorded the number of detected mitochondrial variants, detected variant IDs, top variant, top-variant depth, top-variant ALT count, top-variant VAF, top-variant ALT mean base quality, total ALT count, mean VAF, and maximum VAF.

The selected mtDNA cluster assignment was added to Seurat metadata as a sample-specific mitochondrial cluster column. Cells without mtDNA cluster assignment were labelled as `unassigned`. A sample-identity safety check confirmed that all mtDNA-assigned cells belonged to the expected sample according to the Seurat metadata column `sample`. This prevented accidental cross-sample assignment caused by barcode-prefix or sample-label inconsistencies.

The script generated a Seurat object containing mtDNA lineage metadata, a full metadata table, barcode-matching QC tables, per-cell mitochondrial variant hover tables, UMAP visualizations coloured by mtDNA cluster and selected QC variables, contingency tables comparing mtDNA clusters with sample, patient, Seurat clusters and cell annotations, and boxplots of RNA and mitochondrial QC variables across mtDNA clusters. When the required R packages were available, SVG and interactive HTML UMAP visualizations were also generated, with per-cell mitochondrial variant information included in hover text.

---

## Output management and reproducibility

Each sample was processed in a sample-specific working directory. Complete runs ending at STEP 15 were moved to the final results directory, whereas partial runs were left in the working directory to support controlled reruns. The pipeline stored sample-specific thresholds, including the empirical base-quality threshold, allowing downstream steps to be restarted without changing filtering parameters. Existing final result directories were not overwritten, preventing accidental replacement of completed analyses.

This workflow therefore implemented a reproducible strategy for single-cell mtDNA lineage inference in which mitochondrial reads were selected under strict unique-mapping criteria, base-quality filtering was estimated empirically from the data, candidate variants were filtered using population-level and cell-level evidence, mitochondrial distances were calculated using per-cell/per-variant coverage-aware VAF comparisons, and final mtDNA lineages were interpreted in the context of annotated single-cell transcriptomic data.
