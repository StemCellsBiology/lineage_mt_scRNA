#!/bin/bash
set -euo pipefail

source ~/venvs/scmt_env/bin/activate

SAMTOOLS=~/software/samtools-1.23.1/samtools1231bin/bin/samtools
SEURAT_RDS=~/KSzade/sctanno_ksz/rzeszow/SCT_annotated_KSz.Rds
SEURAT_BARCODE_TABLE=~/KSzade/sctanno_ksz/sctanno_ksz_all_cells_with_sample.tsv

SCRIPT_DIR="BAMs/scripts"
SCRIPT_NEW_DIR="BAMs/scripts_new"

RESULTS_BASE="/data/KSzade/lineages/results"
mkdir -p "$RESULTS_BASE"

if [[ ! -s "$SEURAT_BARCODE_TABLE" ]]; then
    echo "ERROR: Seurat barcode table does not exist or is empty:"
    echo "$SEURAT_BARCODE_TABLE"
    exit 1
fi

#TSNE_PERPLEXITY=30
MIN_MT_READS=200
MAPQ=30
CHR_M="chrM"

# ------------------------------------------------------------
# Publication-like mitochondrial distance threshold
# ------------------------------------------------------------
# The publication used C > 100.
# Here we use a lower threshold because per-cell mtDNA coverage is lower.
#
# Important:
#   The publication formula uses strict C > threshold.
#   Therefore, to obtain "at least 10 reads", use threshold 9:
#       depth > 9  ==  depth >= 10

PUBLICATION_COVERAGE_THRESHOLD="${PUBLICATION_COVERAGE_THRESHOLD:-9}"
PUBLICATION_MIN_SHARED_VARIANTS="${PUBLICATION_MIN_SHARED_VARIANTS:-1}"
PUBLICATION_NA_DISTANCE_ACTION="${PUBLICATION_NA_DISTANCE_ACTION:-one}"

PUBLICATION_REQUIRE_ALT_SIGNAL_CELL="${PUBLICATION_REQUIRE_ALT_SIGNAL_CELL:-1}"
PUBLICATION_MIN_ALT_FOR_CELL="${PUBLICATION_MIN_ALT_FOR_CELL:-2}"

AUTO_K_MIN="${AUTO_K_MIN:-2}"
AUTO_K_MAX="${AUTO_K_MAX:-30}"
AUTO_K_MIN_CLUSTER_SIZE="${AUTO_K_MIN_CLUSTER_SIZE:-2}"

KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD="${KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD:-0.95}"
KMEDOIDS_MIN_CLUSTER_SIZE="${KMEDOIDS_MIN_CLUSTER_SIZE:-2}"
KMEDOIDS_MIN_RETAINED_CLUSTERS="${KMEDOIDS_MIN_RETAINED_CLUSTERS:-1}"
KMEDOIDS_MIN_ASSIGNED_FRACTION="${KMEDOIDS_MIN_ASSIGNED_FRACTION:-0}"
KMEDOIDS_MEMBERSHIP_EXPONENT="${KMEDOIDS_MEMBERSHIP_EXPONENT:-2}"

# ------------------------------------------------------------
# Legacy labels used mainly in output file names
# ------------------------------------------------------------
# In the publication-like STEP 12 / STEP 12B logic, real filtering is controlled by:
#   PUBLICATION_COVERAGE_THRESHOLD
#   PUBLICATION_MIN_ALT_FOR_CELL
#
# use_DPTH and use_ALT are kept mostly as stable output filename labels.
# By default, make them consistent with the publication-like thresholds.

DEFAULT_USE_DPTH="${DEFAULT_USE_DPTH:-$((PUBLICATION_COVERAGE_THRESHOLD + 1))}"
DEFAULT_USE_ALT="${DEFAULT_USE_ALT:-$PUBLICATION_MIN_ALT_FOR_CELL}"

# ------------------------------------------------------------
# Sample control
# ------------------------------------------------------------
#
# DEFAULT_SAMPLES:
#   all samples known to this pipeline.
#
# SAMPLES:
#   selected samples to which START_STEP should apply.
#
# REST_MODE:
#   SKIP:
#       analyze only samples listed in SAMPLES.
#       This is the default when SAMPLES is provided.
#
#   AUTO:
#       samples listed in SAMPLES use START_STEP;
#       all remaining DEFAULT_SAMPLES use AUTO mode.
#
# Examples:
#
#   bash pipeline10.bash
#       Analyze all DEFAULT_SAMPLES with START_STEP=AUTO.
#
#   START_STEP=12 SAMPLES="BM114" bash pipeline10.bash
#       Analyze only BM114 from STEP 12.
#       BM119 and BM120 are skipped.
#
#   START_STEP=12 SAMPLES="BM114" REST_MODE=AUTO bash pipeline10.bash
#       BM114 starts from STEP 12.
#       BM119 and BM120 start in AUTO mode.
#
#   START_STEP=12 SAMPLES="BM114 BM119" REST_MODE=AUTO bash pipeline10.bash
#       BM114 and BM119 start from STEP 12.
#       BM120 starts in AUTO mode.
#
# ------------------------------------------------------------

DEFAULT_SAMPLES=(3 4 95 99 26ctrl 153ctri 173ctrl 180ctrl A05234 A00103 A03211 Y06102 BM106 BM107 BM108 BM102 BM114 BM115 BM119 BM120)

START_STEP="${START_STEP:-AUTO}"
STOP_STEP="${STOP_STEP:-999}"
REST_MODE="${REST_MODE:-SKIP}"

# Optional per-sample START_STEP override.
# Format examples:
#   SAMPLE_START_STEPS="BM106:15"
#   SAMPLE_START_STEPS="BM106:15 BM107:2 BM108:11"
SAMPLE_START_STEPS="${SAMPLE_START_STEPS:-}"

if [[ "$START_STEP" != "AUTO" && ! "$START_STEP" =~ ^[0-9]+$ ]]; then
    echo "ERROR: START_STEP must be AUTO or an integer, got: $START_STEP"
    exit 1
fi

if [[ ! "$STOP_STEP" =~ ^[0-9]+$ ]]; then
    echo "ERROR: STOP_STEP must be an integer, got: $STOP_STEP"
    exit 1
fi

if [[ "$REST_MODE" != "SKIP" && "$REST_MODE" != "AUTO" ]]; then
    echo "ERROR: REST_MODE must be SKIP or AUTO, got: $REST_MODE"
    exit 1
fi

SELECTED_SAMPLES=()

if [[ -n "${SAMPLES:-}" ]]; then
    read -r -a SELECTED_SAMPLES <<< "$SAMPLES"
else
    # If SAMPLES is not provided, all samples use START_STEP.
    SELECTED_SAMPLES=("${DEFAULT_SAMPLES[@]}")
fi

sample_is_selected() {
    local query="$1"
    local s

    for s in "${SELECTED_SAMPLES[@]}"
    do
        if [[ "$s" == "$query" ]]; then
            return 0
        fi
    done

    return 1
}

get_sample_specific_start_step() {
    local query="$1"
    local pair
    local sample
    local step

    for pair in $SAMPLE_START_STEPS
    do
        sample="${pair%%:*}"
        step="${pair##*:}"

        if [[ "$sample" == "$query" ]]; then
            if [[ "$step" != "AUTO" && ! "$step" =~ ^[0-9]+$ ]]; then
                echo "ERROR: invalid step for sample $sample in SAMPLE_START_STEPS: $step" >&2
                exit 1
            fi

            echo "$step"
            return 0
        fi
    done

    return 1
}

decide_auto_start_step() {
    local outdir="$1"
    local filename="$2"

    if [[ -s "${outdir}/${filename}.uniq_mt.bam" && -s "${outdir}/${filename}.uniq_mt.bam.bai" ]]; then
        echo "2"
    else
        echo "1"
    fi
}

SAMPLE_START_STEP=""

RUN_ONLY_STEPS="${RUN_ONLY_STEPS:-}"

run_step() {
    local step="$1"
    local s

    if [[ -n "$RUN_ONLY_STEPS" ]]; then
        for s in $RUN_ONLY_STEPS
        do
            if [[ "$s" == "$step" ]]; then
                return 0
            fi
        done
        return 1
    fi

    (( step >= SAMPLE_START_STEP && step <= STOP_STEP ))
}

echo "Requested START_STEP: $START_STEP"
echo "STOP_STEP: $STOP_STEP"
echo "REST_MODE: $REST_MODE"
echo "SAMPLE_START_STEPS: ${SAMPLE_START_STEPS:-none}"
echo "DEFAULT_USE_DPTH: $DEFAULT_USE_DPTH"
echo "DEFAULT_USE_ALT: $DEFAULT_USE_ALT"
echo "DEFAULT_SAMPLES: ${DEFAULT_SAMPLES[*]}"
echo "SELECTED_SAMPLES: ${SELECTED_SAMPLES[*]}"
echo "Final results directory: $RESULTS_BASE"

pwd

for filename in "${DEFAULT_SAMPLES[@]}"
do
    # ------------------------------------------------------------
    # Decide whether this sample should be processed and with which
    # start mode.
    # ------------------------------------------------------------

    SAMPLE_START_MODE=""

    if sample_specific_step=$(get_sample_specific_start_step "$filename"); then
        SAMPLE_START_MODE="$sample_specific_step"
        echo "Sample $filename has sample-specific START_STEP=$SAMPLE_START_MODE."

    elif sample_is_selected "$filename"; then
        SAMPLE_START_MODE="$START_STEP"
        echo "Sample $filename is selected. It will use global START_STEP=$START_STEP."

    else
        if [[ "$REST_MODE" == "AUTO" ]]; then
            SAMPLE_START_MODE="AUTO"
            echo "Sample $filename is not selected, but REST_MODE=AUTO. It will use AUTO mode."
        else
            echo "Sample $filename is not selected and REST_MODE=SKIP. Skipping."
            echo
            continue
        fi
    fi

    file="BAMs/all_samples_a2m95/${filename}.bam"

    echo "============================================================"
    echo "Processing BAM: $file"
    echo "============================================================"

    bam_dir=$(dirname "$file")

    echo "BAM directory: $bam_dir"
    echo "Sample name: $filename"

    # ------------------------------------------------------------
    # Skip sample if final results already exist
    # ------------------------------------------------------------

    FINAL_OUTDIR="${RESULTS_BASE}/${filename}_mt_pipeline"

    if [[ -d "$FINAL_OUTDIR" ]]; then
        echo "Final results already exist for $filename:"
        echo "$FINAL_OUTDIR"
        echo "Skipping this sample to avoid overwriting completed analysis."
        echo
        continue
    fi

    # ------------------------------------------------------------
    # Output directory next to BAM
    # ------------------------------------------------------------

    OUTDIR="${bam_dir}/${filename}_mt_pipeline"
    mkdir -p "$OUTDIR"

    echo "Output directory: $OUTDIR"

    SKIP_REASON_FILE="${OUTDIR}/${filename}.skip_reason.txt"

    # ------------------------------------------------------------
    # Define/recover fixed thresholds used from STEP 11 onward
    # ------------------------------------------------------------

    # NOTE:
    # use_DPTH and use_ALT are legacy filename labels in STEP 12/12B/15/15B.
    # They do not control publication-like filtering in STEP 12/12B.
    # Actual filtering uses PUBLICATION_COVERAGE_THRESHOLD and PUBLICATION_MIN_ALT_FOR_CELL.

    DPTH_FILE="${OUTDIR}/${filename}.use_DPTH.txt"
    ALT_FILE="${OUTDIR}/${filename}.use_ALT.txt"

    if [[ -s "$DPTH_FILE" ]]; then
        use_DPTH=$(cat "$DPTH_FILE")
    else
        use_DPTH="$DEFAULT_USE_DPTH"
        echo "$use_DPTH" > "$DPTH_FILE"
    fi

    if [[ -s "$ALT_FILE" ]]; then
        use_ALT=$(cat "$ALT_FILE")
    else
        use_ALT="$DEFAULT_USE_ALT"
        echo "$use_ALT" > "$ALT_FILE"
    fi

    if [[ -z "$use_DPTH" || ! "$use_DPTH" =~ ^[0-9]+$ ]]; then
        echo "ERROR: use_DPTH is empty or not an integer:"
        echo "$use_DPTH"
        echo "File: $DPTH_FILE"
        exit 1
    fi

    if [[ -z "$use_ALT" || ! "$use_ALT" =~ ^[0-9]+$ ]]; then
        echo "ERROR: use_ALT is empty or not an integer:"
        echo "$use_ALT"
        echo "File: $ALT_FILE"
        exit 1
    fi

    echo "Using use_DPTH: $use_DPTH"
    echo "Using use_ALT: $use_ALT"
    echo "use_DPTH file: $DPTH_FILE"
    echo "use_ALT file: $ALT_FILE"

    # ------------------------------------------------------------
    # Decide sample-specific START_STEP
    # ------------------------------------------------------------

    if [[ "$SAMPLE_START_MODE" == "AUTO" ]]; then

        SAMPLE_START_STEP=$(decide_auto_start_step "$OUTDIR" "$filename")

        if [[ "$SAMPLE_START_STEP" == "2" ]]; then
            echo "AUTO START_STEP: found existing uniq_mt.bam and .bai"
            echo "Starting sample $filename from STEP 2"
        else
            echo "AUTO START_STEP: uniq_mt.bam or .bai missing"
            echo "Starting sample $filename from STEP 1"
        fi

    else
        SAMPLE_START_STEP="$SAMPLE_START_MODE"
        echo "Manual START_STEP for sample $filename: $SAMPLE_START_STEP"
    fi

    echo "Effective SAMPLE_START_STEP for $filename: $SAMPLE_START_STEP"

    MT_DOWNSTREAM_OK=0
    MT_FUZZY_DOWNSTREAM_OK=0

    # ------------------------------------------------------------
    # Recover STEP 12 / STEP 12B results if downstream QC steps
    # are requested without rerunning clustering.
    # ------------------------------------------------------------

    if run_step 15; then
        BEST_K_FILE_RECOVERY="${OUTDIR}/${filename}.best_k.${use_DPTH}.txt"

        if [[ -s "$BEST_K_FILE_RECOVERY" ]]; then
            BEST_K_RECOVERY=$(cat "$BEST_K_FILE_RECOVERY")
            CLUSTER_FILE_RECOVERY="${OUTDIR}/${filename}.clusters.${use_DPTH}.k${BEST_K_RECOVERY}.tsv"

            if [[ -s "$CLUSTER_FILE_RECOVERY" ]]; then
                echo "Recovered successful STEP 12 result for STEP 15."
                echo "Recovered best_k: $BEST_K_RECOVERY"
                echo "Recovered cluster file: $CLUSTER_FILE_RECOVERY"
                MT_DOWNSTREAM_OK=1
            else
                echo "No valid STEP 12 cluster file found for STEP 15 recovery."
                echo "Expected: $CLUSTER_FILE_RECOVERY"
                MT_DOWNSTREAM_OK=0
            fi
        else
            echo "No STEP 12 best_k file found for STEP 15 recovery."
            echo "Expected: $BEST_K_FILE_RECOVERY"
            MT_DOWNSTREAM_OK=0
        fi
    fi

    if run_step 151; then
        FUZZY_BEST_K_FILE_RECOVERY="${OUTDIR}/${filename}.best_k_fuzzy.${use_DPTH}.txt"

        if [[ -s "$FUZZY_BEST_K_FILE_RECOVERY" ]]; then
            FUZZY_BEST_K_RECOVERY=$(cat "$FUZZY_BEST_K_FILE_RECOVERY")
            FUZZY_CLUSTER_FILE_RECOVERY="${OUTDIR}/${filename}.clusters_fuzzy.${use_DPTH}.k${FUZZY_BEST_K_RECOVERY}.tsv"

            if [[ -s "$FUZZY_CLUSTER_FILE_RECOVERY" ]]; then
                echo "Recovered successful STEP 12B fuzzy result for STEP 15B."
                echo "Recovered fuzzy best_k: $FUZZY_BEST_K_RECOVERY"
                echo "Recovered fuzzy cluster file: $FUZZY_CLUSTER_FILE_RECOVERY"
                MT_FUZZY_DOWNSTREAM_OK=1
            else
                echo "No valid fuzzy cluster file found for STEP 15B recovery."
                echo "Expected: $FUZZY_CLUSTER_FILE_RECOVERY"
                MT_FUZZY_DOWNSTREAM_OK=0
            fi
        else
            echo "No STEP 12B fuzzy best_k file found for STEP 15B recovery."
            echo "Expected: $FUZZY_BEST_K_FILE_RECOVERY"
            MT_FUZZY_DOWNSTREAM_OK=0
        fi
    fi

    # ------------------------------------------------------------
    # Barcode file from the same directory as BAM
    # ------------------------------------------------------------

    barcode_file=""

    if [[ -f "${bam_dir}/barcodes.${filename}.tsv.gz" ]]; then
        barcode_file="${bam_dir}/barcodes.${filename}.tsv.gz"

    elif [[ -f "${bam_dir}/barcodes.${filename}.tsv" ]]; then
        barcode_file="${bam_dir}/barcodes.${filename}.tsv"

    elif [[ -f "${bam_dir}/${filename}.barcodes.tsv.gz" ]]; then
        barcode_file="${bam_dir}/${filename}.barcodes.tsv.gz"

    elif [[ -f "${bam_dir}/${filename}.barcodes.tsv" ]]; then
        barcode_file="${bam_dir}/${filename}.barcodes.tsv"

    elif [[ -f "${bam_dir}/barcodes.tsv.gz" ]]; then
        barcode_file="${bam_dir}/barcodes.tsv.gz"

    elif [[ -f "${bam_dir}/barcodes.tsv" ]]; then
        barcode_file="${bam_dir}/barcodes.tsv"

    elif [[ -f "${bam_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz" ]]; then
        barcode_file="${bam_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz"

    elif [[ -f "${bam_dir}/raw_feature_bc_matrix/barcodes.tsv.gz" ]]; then
        barcode_file="${bam_dir}/raw_feature_bc_matrix/barcodes.tsv.gz"
    fi

    if [[ -z "$barcode_file" ]]; then
        echo "ERROR: Could not find barcode file for sample $filename in directory $bam_dir"
        echo "Checked:"
        echo "  ${bam_dir}/barcodes.${filename}.tsv.gz"
        echo "  ${bam_dir}/barcodes.${filename}.tsv"
        echo "  ${bam_dir}/${filename}.barcodes.tsv.gz"
        echo "  ${bam_dir}/${filename}.barcodes.tsv"
        echo "  ${bam_dir}/barcodes.tsv.gz"
        echo "  ${bam_dir}/barcodes.tsv"
        echo "  ${bam_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz"
        echo "  ${bam_dir}/raw_feature_bc_matrix/barcodes.tsv.gz"
        exit 1
    fi

    echo "Barcode file: $barcode_file"

    echo "$filename" > "${OUTDIR}/${filename}.name.txt"
    echo "$barcode_file" > "${OUTDIR}/${filename}.barcode_file.txt"

    # ------------------------------------------------------------
    # Prepare barcode file as plain text
    # ------------------------------------------------------------

    barcode_file_plain_all="${OUTDIR}/${filename}.barcodes.plain.all_from_cellranger.tsv"
    barcode_file_plain="${OUTDIR}/${filename}.barcodes.from_seurat.tsv"

    if [[ "$barcode_file" == *.gz ]]; then
        echo "Barcode file is gzipped. Decompressing to: $barcode_file_plain_all"
        zcat "$barcode_file" > "$barcode_file_plain_all"
    else
        echo "Barcode file is plain text. Copying to: $barcode_file_plain_all"
        cp "$barcode_file" "$barcode_file_plain_all"
    fi

    echo "Filtering barcodes using Seurat metadata table:"
    echo "$SEURAT_BARCODE_TABLE"
    echo "Sample: $filename"

    awk -v sample="$filename" '
        BEGIN { FS = OFS = "\t" }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "sample") sample_col = i
                if ($i == "barcode") barcode_col = i
            }

            if (sample_col == "" || barcode_col == "") {
                print "ERROR: required columns sample and/or barcode not found" > "/dev/stderr"
                exit 1
            }

            next
        }

        $sample_col == sample {
            print $barcode_col
        }
    ' "$SEURAT_BARCODE_TABLE" | sort -u > "${OUTDIR}/${filename}.barcodes.from_seurat.raw.tsv"

    echo "Intersecting Seurat barcodes with Cell Ranger barcode list..."

    awk '
        NR == FNR {
            keep[$1] = 1
            next
        }

        ($1 in keep) {
            print $1
        }
    ' "${OUTDIR}/${filename}.barcodes.from_seurat.raw.tsv" "$barcode_file_plain_all" \
        > "$barcode_file_plain"

    n_seurat_barcodes=$(wc -l < "${OUTDIR}/${filename}.barcodes.from_seurat.raw.tsv")
    n_cellranger_barcodes=$(wc -l < "$barcode_file_plain_all")
    n_final_barcodes=$(wc -l < "$barcode_file_plain")

    echo "Barcode filtering summary for $filename:"
    echo "  Seurat barcodes for sample:      $n_seurat_barcodes"
    echo "  Cell Ranger barcodes available:  $n_cellranger_barcodes"
    echo "  Final intersected barcodes:      $n_final_barcodes"

    if (( n_final_barcodes == 0 )); then
        echo "ERROR: No overlapping barcodes found for sample $filename"
        echo "This usually means that sample names or barcode suffixes do not match."
        echo "Example Seurat barcodes for this sample:"
        head "${OUTDIR}/${filename}.barcodes.from_seurat.raw.tsv" || true
        echo "Example Cell Ranger barcodes:"
        head "$barcode_file_plain_all" || true
        exit 1
    fi

    echo "$barcode_file_plain" > "${OUTDIR}/${filename}.barcode_file_plain.txt"
    echo "$barcode_file_plain_all" > "${OUTDIR}/${filename}.barcode_file_plain_all_from_cellranger.txt"

    # ------------------------------------------------------------
    # STEP 1. Filter reads mapping uniquely only to MT with min MAPQ
    # ------------------------------------------------------------

    if run_step 1; then
        echo "STEP 1: filter_unique_mt_reads.py"

        if [[ ! -s "$file" ]]; then
            echo "ERROR: BAM file is required for STEP 1 but does not exist or is empty:"
            echo "$file"
            exit 1
        fi

        python3 "${SCRIPT_DIR}/filter_unique_mt_reads.py" \
            "$file" \
            "$CHR_M" \
            "${OUTDIR}/${filename}.uniq_mt.bam" \
            "$MAPQ" \
            > "${OUTDIR}/${filename}.filter_unique_mt_reads.log" 2>&1

        "$SAMTOOLS" index -@ 24 "${OUTDIR}/${filename}.uniq_mt.bam"
    else
        echo "Skipping STEP 1"
    fi

    # ------------------------------------------------------------
    # STEP 2. Count unique MT reads per cell
    # ------------------------------------------------------------

    if run_step 2; then
        echo "STEP 2: count_unique_mt_reads_per_cell.py"

        python3 "${SCRIPT_DIR}/count_unique_mt_reads_per_cell.py" \
            "${OUTDIR}/${filename}.uniq_mt.bam" \
            "$barcode_file_plain" \
            "${OUTDIR}/${filename}.uniq_mt_reads_per_cell.tsv" \
            > "${OUTDIR}/${filename}.count_unique_mt_reads_per_cell.log" 2>&1
    else
        echo "Skipping STEP 2"
    fi

    # ------------------------------------------------------------
    # STEP 3. Choose barcodes with min unique mt reads
    # ------------------------------------------------------------

    if run_step 3; then
        echo "STEP 3: select_barcodes_by_mt_reads.py"

        python3 "${SCRIPT_DIR}/select_barcodes_by_mt_reads.py" \
            "${OUTDIR}/${filename}.uniq_mt_reads_per_cell.tsv" \
            "$MIN_MT_READS" \
            "${OUTDIR}/${filename}.mt${MIN_MT_READS}_barcodes.tsv" \
            > "${OUTDIR}/${filename}.select_barcodes_by_mt_reads.log" 2>&1
    else
        echo "Skipping STEP 3"
    fi

    # ------------------------------------------------------------
    # STEP 4. Count per-cell, per-position, per-allele counts
    # ------------------------------------------------------------

    if run_step 4; then
        echo "STEP 4: per_cell_per_allele_counts1.py"

        python3 "${SCRIPT_DIR}/per_cell_per_allele_counts1.py" \
            "${OUTDIR}/${filename}.uniq_mt.bam" \
            "${OUTDIR}/${filename}.mt${MIN_MT_READS}_barcodes.tsv" \
            "${OUTDIR}/${filename}.per_cell_per_allele_counts1.tsv.gz" \
            > "${OUTDIR}/${filename}.per_cell_per_allele_counts1.log" 2>&1
    else
        echo "Skipping STEP 4"
    fi

    # ------------------------------------------------------------
    # STEP 5. Gaussian plot after counting
    # ------------------------------------------------------------

    if run_step 5; then
        echo "STEP 5: gaussian_plot.py"

        python3 "${SCRIPT_DIR}/gaussian_plot.py" \
            "${OUTDIR}/${filename}.per_cell_per_allele_counts1.tsv.gz" \
            "${OUTDIR}/${filename}" \
            > "${OUTDIR}/${filename}.gaussian_plot.log" 2>&1
    else
        echo "Skipping STEP 5"
    fi

    # ------------------------------------------------------------
    # STEP 6. Check BQ threshold
    # ------------------------------------------------------------

    if run_step 6; then
        echo "STEP 6: choose BQ threshold"

        QC_FILE="${OUTDIR}/${filename}.GMM_meanBQ_QC.tsv.gz"
        BQ_FILE="${OUTDIR}/${filename}.BQth.txt"
        USE_BQ_FILE="${OUTDIR}/${filename}.use_BQth.txt"

        if [[ ! -s "$QC_FILE" ]]; then
            echo "ERROR: QC file does not exist or is empty:"
            echo "$QC_FILE"
            exit 1
        fi

        bqth=$(
            zcat "$QC_FILE" | awk '
                $6 == "YES" {
                    if (found == 0 || $4 < min) {
                        min = $4
                        found = 1
                    }
                }
                END {
                    if (found == 1) {
                        print min
                    }
                }
            '
        )

        if [[ -z "$bqth" ]]; then
            echo "ERROR: BQ threshold is empty for $filename"
            echo "No rows with column 6 == YES were found in:"
            echo "$QC_FILE"
            echo
            echo "First lines of QC file:"
            zcat "$QC_FILE" | head || true
            exit 1
        fi

        echo "$bqth" > "$BQ_FILE"
        echo "BQ threshold raw value: $bqth"

        bqth_int=$(awk -v x="$bqth" 'BEGIN { printf "%d", x }')
        use_bqth="$bqth_int"

        if (( use_bqth < 0 )); then
            echo "ERROR: calculated use_bqth < 0: $use_bqth"
            exit 1
        fi

        echo "$use_bqth" > "$USE_BQ_FILE"

        echo "BQ threshold integer value: $bqth_int"
        echo "Using BQ threshold: $use_bqth"
    else
        echo "Skipping STEP 6"
    fi

    # Recover use_bqth only when downstream steps that still require BQ-filtered
    # intermediate files may be executed.
    #
    # STEP 12 and STEP 12B do not need use_bqth directly. They start from:
    #   *.filtered_variants_by_cell.long.tsv.gz
    #
    # Therefore, if START_STEP is 12 or higher, especially 121, do not require
    # use_BQth.txt.

    if (( SAMPLE_START_STEP > 6 && SAMPLE_START_STEP < 12 )); then
        bq_file="${OUTDIR}/${filename}.use_BQth.txt"

        if [[ ! -f "$bq_file" ]]; then
            echo "ERROR: SAMPLE_START_STEP > 6 and < 12 but BQ threshold file is missing:"
            echo "$bq_file"
            exit 1
        fi

        use_bqth=$(cat "$bq_file")

        if [[ -z "$use_bqth" ]]; then
            echo "ERROR: use_bqth is empty in $bq_file"
            exit 1
        fi

        echo "Recovered use_bqth from file: $use_bqth"

    elif (( SAMPLE_START_STEP >= 12 )); then
        echo "START_STEP >= 12, so use_BQth.txt is not required."
        echo "STEP 12/12B use filtered_variants_by_cell.long.tsv.gz directly."
    fi

    # ------------------------------------------------------------
    # STEP 7. Filtering by BQ threshold
    # ------------------------------------------------------------

    if run_step 7; then
        echo "STEP 7: filter_per_cell_allele_counts_by_bq_plus.py"

        python3 "${SCRIPT_DIR}/filter_per_cell_allele_counts_by_bq_plus.py" \
            "${OUTDIR}/${filename}.per_cell_per_allele_counts1.tsv.gz" \
            "$use_bqth" \
            "${OUTDIR}/${filename}.filtered_plus" \
            --low-impact-af 0.002 \
            > "${OUTDIR}/${filename}.filter_per_cell_allele_counts_by_bq_plus.log" 2>&1
    else
        echo "Skipping STEP 7"
    fi

    # ------------------------------------------------------------
    # STEP 8. Extract candidate alleles
    # ------------------------------------------------------------

    if run_step 8; then
        echo "STEP 8: call_candidate_alleles.py"

        python3 "${SCRIPT_DIR}/call_candidate_alleles.py" \
            "${OUTDIR}/${filename}.filtered_plus.BQfiltered.per_cell_alleles.tsv.gz" \
            "${OUTDIR}/${filename}.candidate_alleles.tsv" \
            "${OUTDIR}/${filename}.candidate_alleles_by_cell.long.tsv.gz" \
            > "${OUTDIR}/${filename}.call_candidate_alleles.log" 2>&1
    else
        echo "Skipping STEP 8"
    fi

    # ------------------------------------------------------------
    # STEP 9. Filter variants
    # ------------------------------------------------------------

    if run_step 9; then
        echo "STEP 9: filter_variants_no_bulk_after_bq_plus.py"

        python3 "${SCRIPT_DIR}/filter_variants_no_bulk_after_bq_plus.py" \
            "${OUTDIR}/${filename}.candidate_alleles.tsv" \
            "${OUTDIR}/${filename}.candidate_alleles_by_cell.long.tsv.gz" \
            "${OUTDIR}/${filename}.filtered_variants.no_bulk.tsv" \
            "${OUTDIR}/${filename}.filtered_variants_by_cell.long.tsv.gz" \
            --min-heteroplasmy 0.025 \
            --min-cells-alt-positive 1 \
            --min-total-alt-count 1 \
            > "${OUTDIR}/${filename}.filter_variants_no_bulk_after_bq_plus.log" 2>&1
    else
        echo "Skipping STEP 9"
    fi

    # ------------------------------------------------------------
    # STEP 11. Make heteroplasmy matrix
    # ------------------------------------------------------------

    if run_step 11; then
        echo "STEP 11: disabled in publication-like pipeline version"
        echo "Skipping STEP 11 intentionally."
    else
        echo "Skipping STEP 11"
    fi

    # ------------------------------------------------------------
    # STEP 12. Publication-like mitochondrial distance clustering
    #          with automatic single-k selection
    # ------------------------------------------------------------
    #
    # This step:
    #   1. uses the long table directly:
    #        *.filtered_variants_by_cell.long.tsv.gz
    #
    #   2. selects cells with real ALT signal if:
    #        PUBLICATION_REQUIRE_ALT_SIGNAL_CELL=1
    #
    #      cell-level inclusion rule:
    #        depth > PUBLICATION_COVERAGE_THRESHOLD
    #        alt_count >= PUBLICATION_MIN_ALT_FOR_CELL
    #        vaf > 0
    #
    #   3. computes publication-like pairwise mitochondrial distance:
    #
    #        D_ij = mean sqrt(abs(VAF_i - VAF_j))
    #
    #      only over variants sufficiently covered in both cells:
    #
    #        depth_i > PUBLICATION_COVERAGE_THRESHOLD
    #        depth_j > PUBLICATION_COVERAGE_THRESHOLD
    #
    #   4. performs:
    #
    #        hclust(D, method = "average")
    #
    #   5. automatically selects one best_k using mean silhouette.
    #
    #   6. writes only one best-k cluster file:
    #
    #        ${filename}.clusters.${use_DPTH}.k<best_k>.tsv
    #
    # ------------------------------------------------------------

    if run_step 12; then
        echo "STEP 12: publication-like mitochondrial distance clustering with automatic best_k selection"

        FILTERED_LONG_FILE="${OUTDIR}/${filename}.filtered_variants_by_cell.long.tsv.gz"

        VALID_K_FILE="${OUTDIR}/${filename}.valid_k.${use_DPTH}.txt"
        BEST_K_FILE="${OUTDIR}/${filename}.best_k.${use_DPTH}.txt"
        MAX_K_FILE="${OUTDIR}/${filename}.max_possible_k.${use_DPTH}.txt"
        AUTO_K_SELECTION_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.auto.auto_k_selection.tsv"

        PUBLICATION_CLUSTER_SCRIPT="${SCRIPT_NEW_DIR}/plot_heatmap_and_cluster_publication_mito_distance_auto_k_hclust_REWRITTEN_retained_tiebreak.R"

        if [[ ! -s "$FILTERED_LONG_FILE" ]]; then
            echo "ERROR: filtered variants long file does not exist or is empty:"
            echo "$FILTERED_LONG_FILE"
            exit 1
        fi

        if [[ ! -s "$PUBLICATION_CLUSTER_SCRIPT" ]]; then
            echo "ERROR: auto-k publication-like STEP 12 script does not exist or is empty:"
            echo "$PUBLICATION_CLUSTER_SCRIPT"
            exit 1
        fi

        echo "Publication-like coverage threshold: $PUBLICATION_COVERAGE_THRESHOLD"
        echo "Coverage rule used by script: depth > threshold"
        echo "Publication-like min shared variants: $PUBLICATION_MIN_SHARED_VARIANTS"
        echo "Publication-like NA distance action: $PUBLICATION_NA_DISTANCE_ACTION"
        echo "Require ALT-positive cell for publication-like clustering: $PUBLICATION_REQUIRE_ALT_SIGNAL_CELL"
        echo "Minimum ALT count for cell inclusion: $PUBLICATION_MIN_ALT_FOR_CELL"
        echo "AUTO_K_MIN: $AUTO_K_MIN"
        echo "AUTO_K_MAX: $AUTO_K_MAX"
        echo "AUTO_K_MIN_CLUSTER_SIZE: $AUTO_K_MIN_CLUSTER_SIZE"

        echo "Removing old downstream cluster files for this depth..."

        rm -f "${OUTDIR}/${filename}.clusters.${use_DPTH}.k"*.tsv
        rm -f "${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.k"*.log
        rm -f "${OUTDIR}/${filename}.valid_k.${use_DPTH}.txt"
        rm -f "${OUTDIR}/${filename}.best_k.${use_DPTH}.txt"
        rm -f "${OUTDIR}/${filename}.max_possible_k.${use_DPTH}.txt"
        rm -f "${OUTDIR}/${filename}.dynamic_k_list.${use_DPTH}.txt"
        rm -f "${OUTDIR}/${filename}.publication_mito.${use_DPTH}."*
        rm -f "${OUTDIR}/${filename}.mt_compare"*
        rm -f "${OUTDIR}/${filename}.mito_tsne.${use_DPTH}."*
        rm -rf "${OUTDIR}/${filename}.mt.k"*
        rm -rf "${OUTDIR}/${filename}.mt.multik"*
        rm -rf "${OUTDIR}/${filename}.mt.publication_assigned_only"*

        echo "Running publication-like mitochondrial distance clustering with automatic best_k selection"

        if Rscript "$PUBLICATION_CLUSTER_SCRIPT" \
            "$FILTERED_LONG_FILE" \
            "${OUTDIR}/${filename}.publication_mito.${use_DPTH}.auto" \
            "${OUTDIR}/${filename}.clusters.${use_DPTH}" \
            "$PUBLICATION_COVERAGE_THRESHOLD" \
            "$PUBLICATION_MIN_SHARED_VARIANTS" \
            "$PUBLICATION_NA_DISTANCE_ACTION" \
            "$PUBLICATION_REQUIRE_ALT_SIGNAL_CELL" \
            "$PUBLICATION_MIN_ALT_FOR_CELL" \
            "$AUTO_K_MIN" \
            "$AUTO_K_MAX" \
            "$AUTO_K_MIN_CLUSTER_SIZE" \
            > "${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.auto_k.log" 2>&1
        then
            R_BEST_K_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.auto.best_k.txt"

            if [[ ! -s "$R_BEST_K_FILE" ]]; then
                echo "ERROR: auto-k script finished but did not create best_k file:"
                echo "$R_BEST_K_FILE"
                exit 1
            fi

            BEST_K=$(cat "$R_BEST_K_FILE")

            if [[ -z "$BEST_K" || ! "$BEST_K" =~ ^[0-9]+$ ]]; then
                echo "ERROR: BEST_K is empty or not an integer:"
                echo "$BEST_K"
                exit 1
            fi

            BEST_CLUSTER_FILE="${OUTDIR}/${filename}.clusters.${use_DPTH}.k${BEST_K}.tsv"

            if [[ ! -s "$BEST_CLUSTER_FILE" ]]; then
                echo "ERROR: best-k cluster file does not exist or is empty:"
                echo "$BEST_CLUSTER_FILE"
                exit 1
            fi

            echo "$BEST_K" > "$BEST_K_FILE"
            echo "$BEST_K" > "$VALID_K_FILE"

            n_best_cells=$(
                awk 'NR > 1 { n++ } END { print n + 0 }' "$BEST_CLUSTER_FILE"
            )

            echo "$n_best_cells" > "$MAX_K_FILE"

            echo "Finished STEP 12 successfully"
            echo "Selected best_k: $BEST_K"
            echo "Cells in best-k cluster file: $n_best_cells"
            echo "Best-k cluster file: $BEST_CLUSTER_FILE"
            echo "Auto-k selection table:"
            echo "$AUTO_K_SELECTION_FILE"

            MT_DOWNSTREAM_OK=1

        else
            echo "ERROR: publication-like auto-k clustering failed for $filename."
            echo "See log:"
            echo "${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.auto_k.log"

            {
                echo "sample: $filename"
                echo "reason: publication-like auto-k clustering failed"
                echo "log: ${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.auto_k.log"
                echo "filtered_long_file: $FILTERED_LONG_FILE"
            } > "$SKIP_REASON_FILE"

            MT_DOWNSTREAM_OK=0
        fi

    else
        echo "Skipping STEP 12"
    fi
    
    # ------------------------------------------------------------
    # STEP 12B. Publication-like mitochondrial distance clustering
    #           using fuzzy k-medoids with high-confidence assignment
    # ------------------------------------------------------------
    #
    # This step is an additional clustering mode and does not replace STEP 12.
    # STEP 12 remains the hclust + auto-k version.
    #
    # STEP 12B:
    #   1. uses the same long table as STEP 12:
    #        *.filtered_variants_by_cell.long.tsv.gz
    #
    #   2. applies the same publication-like cell and coverage filters:
    #        depth > PUBLICATION_COVERAGE_THRESHOLD
    #        alt_count >= PUBLICATION_MIN_ALT_FOR_CELL
    #        vaf > 0
    #
    #   3. computes the same publication-like pairwise mitochondrial distance:
    #        D_ij = mean sqrt(abs(VAF_i - VAF_j))
    #
    #   4. performs fuzzy k-medoids on the distance matrix for candidate K values.
    #
    #   5. assigns cells only if:
    #        max fuzzy membership >= KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD
    #
    #   6. clusters with fewer than KMEDOIDS_MIN_CLUSTER_SIZE assigned cells
    #      are not retained as high-confidence clusters; their cells are marked
    #      as unassigned.
    #
    #   7. selects best_k among candidate K values using mean silhouette among
    #      retained assigned cells. K values with fewer than
    #      KMEDOIDS_MIN_RETAINED_CLUSTERS retained clusters are not valid
    #      candidates for best_k.
    #
    #   8. writes an assigned-only cluster file compatible with STEP 15:
    #        ${filename}.clusters_fuzzy.${use_DPTH}.k<best_k>.tsv
    #
    #   9. also writes full assignment and membership tables:
    #        *.fuzzy_kmedoids.assignments.k<best_k>.tsv
    #        *.fuzzy_kmedoids.membership.k<best_k>.tsv
    #
    # Recommended variables to define near the publication-like threshold block:
    #
    #   KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD="${KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD:-0.95}"
    #   KMEDOIDS_MIN_CLUSTER_SIZE="${KMEDOIDS_MIN_CLUSTER_SIZE:-2}"
    #   KMEDOIDS_MIN_RETAINED_CLUSTERS="${KMEDOIDS_MIN_RETAINED_CLUSTERS:-2}"
    #   KMEDOIDS_MIN_ASSIGNED_FRACTION="${KMEDOIDS_MIN_ASSIGNED_FRACTION:-0}"
    #   KMEDOIDS_MEMBERSHIP_EXPONENT="${KMEDOIDS_MEMBERSHIP_EXPONENT:-2}"
    #
    # ------------------------------------------------------------
    
    if run_step 121; then
        echo "STEP 12B: publication-like mitochondrial distance clustering with fuzzy k-medoids"
    
        FILTERED_LONG_FILE="${OUTDIR}/${filename}.filtered_variants_by_cell.long.tsv.gz"
    
        FUZZY_BEST_K_FILE="${OUTDIR}/${filename}.best_k_fuzzy.${use_DPTH}.txt"
        FUZZY_VALID_K_FILE="${OUTDIR}/${filename}.valid_k_fuzzy.${use_DPTH}.txt"
        FUZZY_AUTO_K_SELECTION_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy.auto_k_selection.tsv"
    
        FUZZY_CLUSTER_SCRIPT="${SCRIPT_NEW_DIR}/plot_heatmap_and_cluster_publication_mito_distance_fuzzy_kmedoids_auto_k_FINAL_retained_cluster_tiebreak.R"
    
        if [[ ! -s "$FILTERED_LONG_FILE" ]]; then
            echo "ERROR: filtered variants long file does not exist or is empty:"
            echo "$FILTERED_LONG_FILE"
            exit 1
        fi
    
        if [[ ! -s "$FUZZY_CLUSTER_SCRIPT" ]]; then
            echo "ERROR: fuzzy k-medoids STEP 12B script does not exist or is empty:"
            echo "$FUZZY_CLUSTER_SCRIPT"
            exit 1
        fi
    
        echo "Fuzzy k-medoids coverage threshold: $PUBLICATION_COVERAGE_THRESHOLD"
        echo "Coverage rule used by script: depth > threshold"
        echo "Fuzzy k-medoids min shared variants: $PUBLICATION_MIN_SHARED_VARIANTS"
        echo "Fuzzy k-medoids NA distance action: $PUBLICATION_NA_DISTANCE_ACTION"
        echo "Require ALT-positive cell for fuzzy k-medoids clustering: $PUBLICATION_REQUIRE_ALT_SIGNAL_CELL"
        echo "Minimum ALT count for cell inclusion: $PUBLICATION_MIN_ALT_FOR_CELL"
        echo "AUTO_K_MIN: $AUTO_K_MIN"
        echo "AUTO_K_MAX: $AUTO_K_MAX"
        echo "KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD: $KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD"
        echo "KMEDOIDS_MIN_CLUSTER_SIZE: $KMEDOIDS_MIN_CLUSTER_SIZE"
        echo "KMEDOIDS_MIN_RETAINED_CLUSTERS: $KMEDOIDS_MIN_RETAINED_CLUSTERS"
        echo "KMEDOIDS_MIN_ASSIGNED_FRACTION: $KMEDOIDS_MIN_ASSIGNED_FRACTION"
        echo "KMEDOIDS_MEMBERSHIP_EXPONENT: $KMEDOIDS_MEMBERSHIP_EXPONENT"
    
        echo "Removing old STEP 12B fuzzy k-medoids files for this depth..."
    
        rm -f "${OUTDIR}/${filename}.clusters_fuzzy.${use_DPTH}.k"*.tsv
        rm -f "${OUTDIR}/${filename}.best_k_fuzzy.${use_DPTH}.txt"
        rm -f "${OUTDIR}/${filename}.valid_k_fuzzy.${use_DPTH}.txt"
        rm -f "${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy."*
        rm -f "${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.fuzzy_kmedoids_auto_k.log"
    
        echo "Running publication-like fuzzy k-medoids clustering with automatic best_k selection"
    
        if Rscript "$FUZZY_CLUSTER_SCRIPT" \
            "$FILTERED_LONG_FILE" \
            "${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy" \
            "${OUTDIR}/${filename}.clusters_fuzzy.${use_DPTH}" \
            "$PUBLICATION_COVERAGE_THRESHOLD" \
            "$PUBLICATION_MIN_SHARED_VARIANTS" \
            "$PUBLICATION_NA_DISTANCE_ACTION" \
            "$PUBLICATION_REQUIRE_ALT_SIGNAL_CELL" \
            "$PUBLICATION_MIN_ALT_FOR_CELL" \
            "$AUTO_K_MIN" \
            "$AUTO_K_MAX" \
            "$KMEDOIDS_ASSIGNMENT_PROB_THRESHOLD" \
            "$KMEDOIDS_MIN_CLUSTER_SIZE" \
            "$KMEDOIDS_MIN_RETAINED_CLUSTERS" \
            "$KMEDOIDS_MIN_ASSIGNED_FRACTION" \
            "$KMEDOIDS_MEMBERSHIP_EXPONENT" \
            > "${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.fuzzy_kmedoids_auto_k.log" 2>&1
        then
            R_FUZZY_BEST_K_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy.best_k.txt"
    
            if [[ ! -s "$R_FUZZY_BEST_K_FILE" ]]; then
                echo "ERROR: fuzzy k-medoids script finished but did not create best_k file:"
                echo "$R_FUZZY_BEST_K_FILE"
                exit 1
            fi
    
            FUZZY_BEST_K=$(cat "$R_FUZZY_BEST_K_FILE")
    
            if [[ -z "$FUZZY_BEST_K" || ! "$FUZZY_BEST_K" =~ ^[0-9]+$ ]]; then
                echo "ERROR: FUZZY_BEST_K is empty or not an integer:"
                echo "$FUZZY_BEST_K"
                exit 1
            fi
    
            FUZZY_BEST_CLUSTER_FILE="${OUTDIR}/${filename}.clusters_fuzzy.${use_DPTH}.k${FUZZY_BEST_K}.tsv"
            FUZZY_ASSIGNMENT_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy.fuzzy_kmedoids.assignments.k${FUZZY_BEST_K}.tsv"
    
            if [[ ! -s "$FUZZY_BEST_CLUSTER_FILE" ]]; then
                echo "ERROR: fuzzy best-k assigned-only cluster file does not exist or is empty:"
                echo "$FUZZY_BEST_CLUSTER_FILE"
                exit 1
            fi
    
            echo "$FUZZY_BEST_K" > "$FUZZY_BEST_K_FILE"
            echo "$FUZZY_BEST_K" > "$FUZZY_VALID_K_FILE"
    
            n_fuzzy_assigned_cells=$(
                awk 'NR > 1 { n++ } END { print n + 0 }' "$FUZZY_BEST_CLUSTER_FILE"
            )
    
            echo "Finished STEP 12B successfully"
            echo "Selected fuzzy best_k: $FUZZY_BEST_K"
            echo "Fuzzy assigned cells in best-k cluster file: $n_fuzzy_assigned_cells"
            echo "Fuzzy assigned-only cluster file: $FUZZY_BEST_CLUSTER_FILE"
            echo "Fuzzy full assignment file: $FUZZY_ASSIGNMENT_FILE"
            echo "Fuzzy auto-k selection table:"
            echo "$FUZZY_AUTO_K_SELECTION_FILE"
    
            MT_FUZZY_DOWNSTREAM_OK=1
    
        else
            echo "ERROR: fuzzy k-medoids auto-k clustering failed for $filename."
            echo "See log:"
            echo "${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.fuzzy_kmedoids_auto_k.log"
    
            {
                echo "sample: $filename"
                echo "reason: fuzzy k-medoids auto-k clustering failed"
                echo "log: ${OUTDIR}/${filename}.plot_heatmap_and_cluster.${use_DPTH}.fuzzy_kmedoids_auto_k.log"
                echo "filtered_long_file: $FILTERED_LONG_FILE"
            } > "$SKIP_REASON_FILE"
    
            MT_FUZZY_DOWNSTREAM_OK=0
        fi
    
    else
        echo "Skipping STEP 12B"
    fi
    
    # ------------------------------------------------------------
    # STEP 13. Compare mtDNA clusterings across k values
    # ------------------------------------------------------------
    #
    # Disabled in auto-k version because STEP 12 now selects a single best_k.
    #
    # ------------------------------------------------------------

    if run_step 13; then
        echo "STEP 13: disabled in auto-k publication-like pipeline version"
        echo "Reason: STEP 12 selects one best_k, so there are no multiple k values to compare."
    else
        echo "Skipping STEP 13"
    fi

    # ------------------------------------------------------------
    # STEP 14. mito tSNE
    # ------------------------------------------------------------

    if run_step 14; then
        echo "STEP 14: disabled in publication-like pipeline version"
        echo "Skipping STEP 14 intentionally."
    else
        echo "Skipping STEP 14"
    fi

    # ------------------------------------------------------------
    # STEP 15. QC mt clusters on Seurat for one auto-selected best_k
    # ------------------------------------------------------------
    #
    # STEP 12 now selects a single best_k using silhouette score.
    # Therefore STEP 15 must not skip the first k value.
    #
    # Required R script:
    #   ${SCRIPT_NEW_DIR}/qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R
    #
    # The R script is still the multi-k script, but here we pass only one k:
    #   best_k
    #
    # ------------------------------------------------------------

    if run_step 15 && [[ "$MT_DOWNSTREAM_OK" == "1" ]]; then
        echo "STEP 15: qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R for one auto-selected best_k"

        VALID_K_FILE="${OUTDIR}/${filename}.valid_k.${use_DPTH}.txt"
        BEST_K_FILE="${OUTDIR}/${filename}.best_k.${use_DPTH}.txt"
        STEP15_CLUSTER_LIST="${OUTDIR}/${filename}.step15_cluster_files.${use_DPTH}.txt"
        FILTERED_LONG_FILE="${OUTDIR}/${filename}.filtered_variants_by_cell.long.tsv.gz"
        STEP15_R_SCRIPT="${SCRIPT_NEW_DIR}/qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R"

        if [[ ! -s "$BEST_K_FILE" ]]; then
            echo "WARNING: best_k file does not exist or is empty:"
            echo "$BEST_K_FILE"
            echo "Skipping STEP 15 for $filename."

        elif [[ ! -s "$FILTERED_LONG_FILE" ]]; then
            echo "WARNING: filtered variants long file does not exist or is empty:"
            echo "$FILTERED_LONG_FILE"
            echo "Skipping STEP 15 for $filename."

        elif [[ ! -s "$STEP15_R_SCRIPT" ]]; then
            echo "ERROR: STEP 15 R script does not exist or is empty:"
            echo "$STEP15_R_SCRIPT"
            exit 1

        else
            BEST_K=$(cat "$BEST_K_FILE")

            if [[ -z "$BEST_K" || ! "$BEST_K" =~ ^[0-9]+$ ]]; then
                echo "ERROR: BEST_K is empty or not an integer:"
                echo "$BEST_K"
                exit 1
            fi

            MAIN_CLUSTER_FILE="${OUTDIR}/${filename}.clusters.${use_DPTH}.k${BEST_K}.tsv"

            if [[ ! -s "$MAIN_CLUSTER_FILE" ]]; then
                echo "ERROR: best-k cluster file does not exist or is empty:"
                echo "$MAIN_CLUSTER_FILE"
                exit 1
            fi

            echo "$BEST_K" > "$VALID_K_FILE"

            echo -n "" > "$STEP15_CLUSTER_LIST"
            printf "%s\t%s\n" "$BEST_K" "$MAIN_CLUSTER_FILE" >> "$STEP15_CLUSTER_LIST"

            echo "Best k selected for STEP 15:"
            echo "$BEST_K"
            echo "Cluster file selected for STEP 15:"
            cat "$STEP15_CLUSTER_LIST"

            if Rscript "$STEP15_R_SCRIPT" \
                "$SEURAT_RDS" \
                "${OUTDIR}/${filename}.uniq_mt_reads_per_cell.tsv" \
                "$FILTERED_LONG_FILE" \
                "$STEP15_CLUSTER_LIST" \
                "${OUTDIR}/${filename}.mt.best_k" \
                "${filename}_" \
                "${filename}.mt_cluster" \
                sample \
                "$filename" \
                "$PUBLICATION_COVERAGE_THRESHOLD" \
                > "${OUTDIR}/${filename}.qc_mt_clusters_on_seurat.best_k.log" 2>&1
            then
                echo "Finished STEP 15 for auto-selected best_k=$BEST_K for $filename"
            else
                echo "WARNING: STEP 15 failed for $filename."
                echo "See log:"
                echo "${OUTDIR}/${filename}.qc_mt_clusters_on_seurat.best_k.log"

                {
                    echo "sample: $filename"
                    echo "reason: STEP 15 Seurat QC failed"
                    echo "log: ${OUTDIR}/${filename}.qc_mt_clusters_on_seurat.best_k.log"
                    echo "best_k: $BEST_K"
                    echo "cluster_file: $MAIN_CLUSTER_FILE"
                    echo "filtered_long_file: $FILTERED_LONG_FILE"
                } > "$SKIP_REASON_FILE"

                MT_DOWNSTREAM_OK=0
            fi
        fi

    else
        echo "Skipping STEP 15"
        if [[ "${MT_DOWNSTREAM_OK:-1}" != "1" ]]; then
            echo "Reason: publication-like auto-k clustering was not suitable for Seurat mt-cluster annotation"
        fi
    fi
    
    # ------------------------------------------------------------
    # STEP 15B. QC fuzzy mt clusters on Seurat for one auto-selected fuzzy best_k
    # ------------------------------------------------------------
    #
    # This step is the Seurat/QC counterpart of STEP 12B.
    # It does not replace STEP 15.
    #
    # STEP 15 uses the hclust result from STEP 12:
    #   ${filename}.clusters.${use_DPTH}.k<best_k>.tsv
    #
    # STEP 15B uses the fuzzy k-medoids result from STEP 12B:
    #   ${filename}.clusters_fuzzy.${use_DPTH}.k<fuzzy_best_k>.tsv
    #
    # Important:
    #   The fuzzy cluster file is assigned-only. Cells that were not assigned
    #   with high confidence in STEP 12B are not present in this cluster file.
    #   In the Seurat output they should therefore remain unassigned/NA for
    #   the fuzzy mt-cluster metadata column.
    #
    # Required R script:
    #   ${SCRIPT_NEW_DIR}/qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R
    #
    # ------------------------------------------------------------

    if run_step 151 && [[ "${MT_FUZZY_DOWNSTREAM_OK:-0}" == "1" ]]; then
        echo "STEP 15B: qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R for fuzzy k-medoids best_k"

        FUZZY_VALID_K_FILE="${OUTDIR}/${filename}.valid_k_fuzzy.${use_DPTH}.txt"
        FUZZY_BEST_K_FILE="${OUTDIR}/${filename}.best_k_fuzzy.${use_DPTH}.txt"
        STEP15B_CLUSTER_LIST="${OUTDIR}/${filename}.step15B_fuzzy_cluster_files.${use_DPTH}.txt"
        FILTERED_LONG_FILE="${OUTDIR}/${filename}.filtered_variants_by_cell.long.tsv.gz"
        STEP15_R_SCRIPT="${SCRIPT_NEW_DIR}/qc_mt_clusters_on_seurat_safe_terminal_v7_multik.R"

        if [[ ! -s "$FUZZY_BEST_K_FILE" ]]; then
            echo "WARNING: fuzzy best_k file does not exist or is empty:"
            echo "$FUZZY_BEST_K_FILE"
            echo "Skipping STEP 15B for $filename."

        elif [[ ! -s "$FILTERED_LONG_FILE" ]]; then
            echo "WARNING: filtered variants long file does not exist or is empty:"
            echo "$FILTERED_LONG_FILE"
            echo "Skipping STEP 15B for $filename."

        elif [[ ! -s "$STEP15_R_SCRIPT" ]]; then
            echo "ERROR: STEP 15 R script does not exist or is empty:"
            echo "$STEP15_R_SCRIPT"
            exit 1

        else
            FUZZY_BEST_K=$(cat "$FUZZY_BEST_K_FILE")

            if [[ -z "$FUZZY_BEST_K" || ! "$FUZZY_BEST_K" =~ ^[0-9]+$ ]]; then
                echo "ERROR: FUZZY_BEST_K is empty or not an integer:"
                echo "$FUZZY_BEST_K"
                exit 1
            fi

            FUZZY_MAIN_CLUSTER_FILE="${OUTDIR}/${filename}.clusters_fuzzy.${use_DPTH}.k${FUZZY_BEST_K}.tsv"
            FUZZY_ASSIGNMENT_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy.fuzzy_kmedoids.assignments.k${FUZZY_BEST_K}.tsv"
            FUZZY_MEMBERSHIP_FILE="${OUTDIR}/${filename}.publication_mito.${use_DPTH}.fuzzy.fuzzy_kmedoids.membership.k${FUZZY_BEST_K}.tsv"

            if [[ ! -s "$FUZZY_MAIN_CLUSTER_FILE" ]]; then
                echo "ERROR: fuzzy best-k assigned-only cluster file does not exist or is empty:"
                echo "$FUZZY_MAIN_CLUSTER_FILE"
                exit 1
            fi

            echo "$FUZZY_BEST_K" > "$FUZZY_VALID_K_FILE"

            echo -n "" > "$STEP15B_CLUSTER_LIST"
            printf "%s\t%s\n" "$FUZZY_BEST_K" "$FUZZY_MAIN_CLUSTER_FILE" >> "$STEP15B_CLUSTER_LIST"

            n_fuzzy_assigned_cells=$(awk 'NR > 1 { n++ } END { print n + 0 }' "$FUZZY_MAIN_CLUSTER_FILE")

            echo "Fuzzy best k selected for STEP 15B:"
            echo "$FUZZY_BEST_K"
            echo "Fuzzy assigned cells in cluster file: $n_fuzzy_assigned_cells"
            echo "Fuzzy assigned-only cluster file selected for STEP 15B:"
            cat "$STEP15B_CLUSTER_LIST"
            echo "Fuzzy full assignment file:"
            echo "$FUZZY_ASSIGNMENT_FILE"
            echo "Fuzzy membership file:"
            echo "$FUZZY_MEMBERSHIP_FILE"

            if Rscript "$STEP15_R_SCRIPT" \
                "$SEURAT_RDS" \
                "${OUTDIR}/${filename}.uniq_mt_reads_per_cell.tsv" \
                "$FILTERED_LONG_FILE" \
                "$STEP15B_CLUSTER_LIST" \
                "${OUTDIR}/${filename}.mt.fuzzy_best_k" \
                "${filename}_" \
                "${filename}.mt_cluster_fuzzy" \
                sample \
                "$filename" \
                "$PUBLICATION_COVERAGE_THRESHOLD" \
                > "${OUTDIR}/${filename}.qc_mt_clusters_on_seurat.fuzzy_best_k.log" 2>&1
            then
                echo "Finished STEP 15B for fuzzy auto-selected best_k=$FUZZY_BEST_K for $filename"
            else
                echo "WARNING: STEP 15B failed for $filename."
                echo "See log:"
                echo "${OUTDIR}/${filename}.qc_mt_clusters_on_seurat.fuzzy_best_k.log"

                {
                    echo "sample: $filename"
                    echo "reason: STEP 15B Seurat QC failed"
                    echo "log: ${OUTDIR}/${filename}.qc_mt_clusters_on_seurat.fuzzy_best_k.log"
                    echo "fuzzy_best_k: $FUZZY_BEST_K"
                    echo "fuzzy_cluster_file: $FUZZY_MAIN_CLUSTER_FILE"
                    echo "filtered_long_file: $FILTERED_LONG_FILE"
                } > "$SKIP_REASON_FILE"
            
                MT_FUZZY_DOWNSTREAM_OK=0
            fi
        fi

    else
        echo "Skipping STEP 15B"
        if [[ "${MT_FUZZY_DOWNSTREAM_OK:-0}" != "1" ]]; then
            echo "Reason: fuzzy k-medoids STEP 12B was not suitable for Seurat mt-cluster annotation or was not run"
        fi
    fi
    
    echo "Finished: $filename"
    echo

    # ------------------------------------------------------------
    # Move completed sample results from fast SSD to slower storage
    # ------------------------------------------------------------
    #
    # Do not move results after partial runs such as:
    #   START_STEP=12 STOP_STEP=13 ...
    #
    # ------------------------------------------------------------

    if (( SAMPLE_START_STEP > 1 || STOP_STEP < 15 )); then
        echo "This was a partial run."
        echo "START_STEP=$SAMPLE_START_STEP, STOP_STEP=$STOP_STEP"
        echo "Not moving results to final storage."
        echo
        continue
    fi

    echo "Preparing to move results to final storage:"
    echo "  from: $OUTDIR"
    echo "  to:   $FINAL_OUTDIR"

    if [[ -e "$FINAL_OUTDIR" ]]; then
        echo "ERROR: Destination directory already exists:"
        echo "$FINAL_OUTDIR"
        echo "Not moving results to avoid overwriting previous analysis."
        exit 1
    fi

    mv "$OUTDIR" "$RESULTS_BASE/"

    echo "Moved results for $filename to:"
    echo "$FINAL_OUTDIR"
    echo

done

echo "All selected/default samples finished."
