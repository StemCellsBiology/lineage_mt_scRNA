#!/usr/bin/env python3

import sys
import gzip
import csv
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter per-cell per-allele mtDNA count table by empirical base-quality "
            "threshold, with the 1/500 low-impact allele exception. "
            "This should be run before call_candidate_alleles.py."
        )
    )

    parser.add_argument(
        "per_cell_counts",
        help=(
            "Input table from per_cell_per_allele_counts1.py "
            "(barcode, pos, A, C, G, T, depth, A_mean_bq, C_mean_bq, G_mean_bq, T_mean_bq)"
        )
    )

    parser.add_argument(
        "bq_threshold",
        type=float,
        help="Empirical mean base-quality threshold, e.g. 23.8"
    )

    parser.add_argument(
        "out_prefix",
        help="Output prefix"
    )

    parser.add_argument(
        "--low-impact-af",
        type=float,
        default=1 / 500,
        help="Low-impact allele fraction threshold. Default: 1/500 = 0.002"
    )

    return parser.parse_args()


args = parse_args()

in_path = args.per_cell_counts
bq_threshold = args.bq_threshold
out_prefix = args.out_prefix
low_impact_af = args.low_impact_af

bases = ("A", "C", "G", "T")

out_filtered = out_prefix + ".BQfiltered.per_cell_alleles.tsv.gz"
out_qc = out_prefix + ".position_allele_BQ_QC.tsv.gz"
out_bad_positions = out_prefix + ".bad_positions.txt"
out_good_positions = out_prefix + ".good_positions.txt"
out_summary = out_prefix + ".BQ_filter_summary.txt"


# ------------------------------------------------------------
# PASS 1:
# Aggregate per-position, per-allele counts and BQ.
#
# Input table has per-cell values:
# barcode pos A C G T depth A_mean_bq C_mean_bq G_mean_bq T_mean_bq
#
# For each position+allele:
# total_count = sum allele counts across all cells
# total_bq_sum = sum(count_in_cell * mean_bq_in_cell)
# mean_bq = total_bq_sum / total_count
#
# Also compute total depth per position:
# position_depth = sum A+C+G+T across all cells
#
# allele_fraction = total_count / position_depth
# ------------------------------------------------------------

allele_agg = defaultdict(lambda: {"count": 0, "bq_sum": 0.0})
position_depth = defaultdict(int)
all_positions = set()

with gzip.open(in_path, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")

    required = {
        "barcode", "pos", "A", "C", "G", "T", "depth",
        "A_mean_bq", "C_mean_bq", "G_mean_bq", "T_mean_bq"
    }

    missing = required - set(reader.fieldnames)
    if missing:
        sys.stderr.write(
            "Missing required columns: "
            + ", ".join(sorted(missing))
            + "\n"
        )
        sys.exit(1)

    for row in reader:
        pos = int(row["pos"])
        all_positions.add(pos)

        # Use depth column as the per-cell total depth at this position.
        # This should equal A+C+G+T from the same row.
        row_depth = int(row["depth"])
        position_depth[pos] += row_depth

        for b in bases:
            count = int(row[b])
            if count == 0:
                continue

            mean_bq_str = row[f"{b}_mean_bq"]
            if mean_bq_str == "NA" or mean_bq_str == "":
                continue

            mean_bq = float(mean_bq_str)

            allele_agg[(pos, b)]["count"] += count
            allele_agg[(pos, b)]["bq_sum"] += count * mean_bq


# ------------------------------------------------------------
# Determine problematic alleles and bad positions.
#
# Paper-like logic:
#
# A position is removed if it contains one or more alleles with:
#   mean_bq < threshold
# unless that allele has a non-significant effect on heteroplasmy:
#   allele_fraction < 1/500
#
# Therefore:
# problematic allele =
#   mean_bq < threshold AND allele_fraction >= low_impact_af
# ------------------------------------------------------------

qc_rows = []
bad_positions = set()

for pos in sorted(all_positions):
    pos_depth = position_depth[pos]

    for allele in bases:
        key = (pos, allele)

        if key not in allele_agg:
            # Allele not observed at this position.
            continue

        total_count = allele_agg[key]["count"]
        bq_sum = allele_agg[key]["bq_sum"]

        if total_count <= 0 or pos_depth <= 0:
            continue

        mean_bq = bq_sum / total_count
        allele_fraction = total_count / pos_depth

        low_bq = mean_bq < bq_threshold
        low_impact = allele_fraction < low_impact_af
        problematic = low_bq and not low_impact

        if problematic:
            bad_positions.add(pos)

        qc_rows.append({
            "pos": pos,
            "allele": allele,
            "position_depth": pos_depth,
            "allele_count": total_count,
            "allele_fraction": allele_fraction,
            "mean_bq": mean_bq,
            "bq_threshold": bq_threshold,
            "low_impact_af": low_impact_af,
            "low_bq": low_bq,
            "low_impact_below_1over500": low_impact,
            "problematic_for_position": problematic
        })


good_positions = all_positions - bad_positions


# ------------------------------------------------------------
# Save QC table
# ------------------------------------------------------------

with gzip.open(out_qc, "wt") as out:
    out.write(
        "pos\tallele\tposition_depth\tallele_count\tallele_fraction\t"
        "mean_bq\tbq_threshold\tlow_impact_af\t"
        "low_bq\tlow_impact_below_1over500\tproblematic_for_position\n"
    )

    for r in qc_rows:
        out.write(
            f"{r['pos']}\t"
            f"{r['allele']}\t"
            f"{r['position_depth']}\t"
            f"{r['allele_count']}\t"
            f"{r['allele_fraction']:.8f}\t"
            f"{r['mean_bq']:.6f}\t"
            f"{r['bq_threshold']:.6f}\t"
            f"{r['low_impact_af']:.8f}\t"
            f"{r['low_bq']}\t"
            f"{r['low_impact_below_1over500']}\t"
            f"{r['problematic_for_position']}\n"
        )


with open(out_bad_positions, "w") as out:
    for pos in sorted(bad_positions):
        out.write(f"{pos}\n")


with open(out_good_positions, "w") as out:
    for pos in sorted(good_positions):
        out.write(f"{pos}\n")


# ------------------------------------------------------------
# PASS 2:
# Filter original per-cell table.
#
# Conservative position-level filtering:
# remove all rows for positions in bad_positions.
# ------------------------------------------------------------

kept_rows = 0
removed_rows = 0

with gzip.open(in_path, "rt") as f, gzip.open(out_filtered, "wt") as out:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames

    out.write("\t".join(fieldnames) + "\n")

    for row in reader:
        pos = int(row["pos"])

        if pos in bad_positions:
            removed_rows += 1
            continue

        out.write("\t".join(row[col] for col in fieldnames) + "\n")
        kept_rows += 1


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

n_low_bq = sum(1 for r in qc_rows if r["low_bq"])
n_low_impact = sum(1 for r in qc_rows if r["low_impact_below_1over500"])
n_problematic = sum(1 for r in qc_rows if r["problematic_for_position"])

with open(out_summary, "w") as out:
    out.write(f"input_table\t{in_path}\n")
    out.write(f"bq_threshold\t{bq_threshold}\n")
    out.write(f"low_impact_af\t{low_impact_af}\n")
    out.write(f"total_positions\t{len(all_positions)}\n")
    out.write(f"good_positions\t{len(good_positions)}\n")
    out.write(f"bad_positions\t{len(bad_positions)}\n")
    out.write(f"observed_position_alleles\t{len(qc_rows)}\n")
    out.write(f"low_bq_position_alleles\t{n_low_bq}\n")
    out.write(f"low_impact_position_alleles\t{n_low_impact}\n")
    out.write(f"problematic_position_alleles\t{n_problematic}\n")
    out.write(f"rows_kept\t{kept_rows}\n")
    out.write(f"rows_removed\t{removed_rows}\n")
    out.write(f"filtered_table\t{out_filtered}\n")
    out.write(f"qc_table\t{out_qc}\n")
    out.write(f"bad_positions_file\t{out_bad_positions}\n")
    out.write(f"good_positions_file\t{out_good_positions}\n")


sys.stderr.write(f"Input table: {in_path}\n")
sys.stderr.write(f"BQ threshold: {bq_threshold}\n")
sys.stderr.write(f"Low-impact AF threshold: {low_impact_af}\n")
sys.stderr.write(f"Filtered table: {out_filtered}\n")
sys.stderr.write(f"QC table: {out_qc}\n")
sys.stderr.write(f"Bad positions: {out_bad_positions}\n")
sys.stderr.write(f"Good positions: {out_good_positions}\n")
sys.stderr.write(f"Summary: {out_summary}\n")
sys.stderr.write(f"Total positions: {len(all_positions)}\n")
sys.stderr.write(f"Good positions: {len(good_positions)}\n")
sys.stderr.write(f"Bad positions: {len(bad_positions)}\n")
sys.stderr.write(f"Observed position-alleles: {len(qc_rows)}\n")
sys.stderr.write(f"Low-BQ position-alleles: {n_low_bq}\n")
sys.stderr.write(f"Low-impact position-alleles: {n_low_impact}\n")
sys.stderr.write(f"Problematic position-alleles: {n_problematic}\n")
sys.stderr.write(f"Rows kept: {kept_rows}\n")
sys.stderr.write(f"Rows removed: {removed_rows}\n")

