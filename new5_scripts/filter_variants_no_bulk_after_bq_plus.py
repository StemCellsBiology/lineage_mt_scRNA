#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter candidate mtDNA variants after prior BQ filtering. "
            "This script assumes that base-quality filtering was already done "
            "upstream on the per-cell allele count table."
        )
    )

    parser.add_argument(
        "candidate_alleles",
        help="Input candidate_alleles.tsv from call_candidate_alleles.py"
    )

    parser.add_argument(
        "long_by_cell",
        help="Input long_by_cell.tsv.gz from call_candidate_alleles.py"
    )

    parser.add_argument(
        "filtered_candidates",
        help="Output filtered candidate variants TSV"
    )

    parser.add_argument(
        "filtered_long",
        help="Output filtered long-by-cell TSV.GZ"
    )

    parser.add_argument(
        "--min-heteroplasmy",
        type=float,
        default=0.005,
        help="Minimum population heteroplasmy / population VAF. Default: 0.005"
    )

    parser.add_argument(
        "--min-cells-alt-positive",
        type=int,
        default=1,
        help="Minimum number of cells with alt_count > 0. Default: 1"
    )

    parser.add_argument(
        "--min-total-alt-count",
        type=int,
        default=1,
        help="Minimum total ALT read count across all cells. Default: 1"
    )

    parser.add_argument(
        "--min-total-depth",
        type=int,
        default=1,
        help="Minimum total depth at the position across all cells. Default: 1"
    )

    return parser.parse_args()


args = parse_args()

cand = pd.read_csv(args.candidate_alleles, sep="\t")
long = pd.read_csv(args.long_by_cell, sep="\t")


required_cand_cols = {
    "variant",
    "mean_heteroplasmy_population",
    "n_cells_alt_positive",
    "total_alt_count",
    "total_depth"
}

missing = required_cand_cols - set(cand.columns)
if missing:
    sys.stderr.write(
        "Missing required columns in candidate table: "
        + ", ".join(sorted(missing))
        + "\n"
    )
    sys.exit(1)

if "variant" not in long.columns:
    sys.stderr.write("Missing required column in long table: variant\n")
    sys.exit(1)


# Ensure numeric columns
numeric_cols = [
    "mean_heteroplasmy_population",
    "n_cells_alt_positive",
    "total_alt_count",
    "total_depth"
]

for col in numeric_cols:
    cand[col] = pd.to_numeric(cand[col], errors="coerce")


pass_heteroplasmy = (
    cand["mean_heteroplasmy_population"] >= args.min_heteroplasmy
)

pass_cells = (
    cand["n_cells_alt_positive"] >= args.min_cells_alt_positive
)

pass_alt_count = (
    cand["total_alt_count"] >= args.min_total_alt_count
)

pass_depth = (
    cand["total_depth"] >= args.min_total_depth
)

final_pass = (
    pass_heteroplasmy &
    pass_cells &
    pass_alt_count &
    pass_depth
)

filtered = cand[final_pass].copy()

# Diagnostic QC table
qc = cand.copy()
qc["pass_min_heteroplasmy"] = pass_heteroplasmy
qc["pass_min_cells_alt_positive"] = pass_cells
qc["pass_min_total_alt_count"] = pass_alt_count
qc["pass_min_total_depth"] = pass_depth
qc["pass_final"] = final_pass

qc_out = args.filtered_candidates.replace(".tsv", ".QC_all_candidates.tsv")
qc.to_csv(qc_out, sep="\t", index=False)

keep = set(filtered["variant"])
long_f = long[long["variant"].isin(keep)].copy()

filtered.to_csv(args.filtered_candidates, sep="\t", index=False)
long_f.to_csv(args.filtered_long, sep="\t", index=False, compression="gzip")


sys.stderr.write(f"candidate_file={args.candidate_alleles}\n")
sys.stderr.write(f"long_file={args.long_by_cell}\n")
sys.stderr.write(f"min_heteroplasmy={args.min_heteroplasmy}\n")
sys.stderr.write(f"min_cells_alt_positive={args.min_cells_alt_positive}\n")
sys.stderr.write(f"min_total_alt_count={args.min_total_alt_count}\n")
sys.stderr.write(f"min_total_depth={args.min_total_depth}\n")
sys.stderr.write(f"variants_before={cand.shape[0]}\n")
sys.stderr.write(f"variants_after={filtered.shape[0]}\n")
sys.stderr.write(f"long_rows_before={long.shape[0]}\n")
sys.stderr.write(f"long_rows_after={long_f.shape[0]}\n")
sys.stderr.write(f"filtered_candidates={args.filtered_candidates}\n")
sys.stderr.write(f"filtered_long={args.filtered_long}\n")
sys.stderr.write(f"qc_all_candidates={qc_out}\n")

