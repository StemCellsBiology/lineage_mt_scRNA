#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: call_candidate_alleles.py <per_cell_counts.tsv.gz> <candidate_alleles.tsv> <long_by_cell.tsv.gz>\n"
    )
    sys.exit(1)

counts_file = sys.argv[1]
candidate_out = sys.argv[2]
long_out = sys.argv[3]

df = pd.read_csv(counts_file, sep="\t")
bases = ["A", "C", "G", "T"]

candidate_rows = []
long_rows = []

pooled = df.groupby("pos")[bases].sum().reset_index()

for _, row in pooled.iterrows():
    pos = int(row["pos"])
    total_depth = int(sum(row[b] for b in bases))

    if total_depth == 0:
        continue

    sorted_bases = sorted(
        [(b, int(row[b])) for b in bases],
        key=lambda x: x[1],
        reverse=True
    )

    ref_like = sorted_bases[0][0]
    pos_df = df[df["pos"] == pos].copy()
    pos_total_depth = pos_df["depth"].sum()

    for alt, total_alt_count in sorted_bases[1:]:
        if total_alt_count == 0:
            continue

        alt_bq_col = f"{alt}_mean_bq"

        tmp = pos_df[["barcode", "pos", "depth", alt, alt_bq_col]].copy()
        tmp = tmp.rename(columns={
            alt: "alt_count",
            alt_bq_col: "alt_mean_bq"
        })

        tmp["vaf"] = tmp["alt_count"] / tmp["depth"]
        tmp["variant"] = f"{pos}_{ref_like}>{alt}"
        tmp["ref_like"] = ref_like
        tmp["alt"] = alt

        mean_population_af = (
            tmp["alt_count"].sum() / pos_total_depth
            if pos_total_depth > 0 else 0
        )

        observed_bq = pd.to_numeric(
            tmp.loc[tmp["alt_count"] > 0, "alt_mean_bq"],
            errors="coerce"
        )

        candidate_rows.append({
            "variant": f"{pos}_{ref_like}>{alt}",
            "pos": pos,
            "ref_like": ref_like,
            "alt": alt,
            "total_depth": int(pos_total_depth),
            "total_alt_count": int(tmp["alt_count"].sum()),
            "mean_heteroplasmy_population": mean_population_af,
            "n_cells_alt_positive": int((tmp["alt_count"] > 0).sum()),
            "mean_alt_bq_observed": (
                float(observed_bq.mean())
                if len(observed_bq.dropna()) > 0 else np.nan
            )
        })

        long_rows.append(
            tmp[[
                "barcode", "variant", "pos", "ref_like", "alt",
                "depth", "alt_count", "vaf", "alt_mean_bq"
            ]]
        )

candidates = pd.DataFrame(candidate_rows)

if not candidates.empty:
    candidates = candidates.sort_values(
        ["mean_heteroplasmy_population", "n_cells_alt_positive"],
        ascending=[False, False]
    )

candidates.to_csv(candidate_out, sep="\t", index=False)

if long_rows:
    long = pd.concat(long_rows, ignore_index=True)
else:
    long = pd.DataFrame(columns=[
        "barcode", "variant", "pos", "ref_like", "alt",
        "depth", "alt_count", "vaf", "alt_mean_bq"
    ])

long.to_csv(long_out, sep="\t", index=False, compression="gzip")
