#!/usr/bin/env python3

import sys
import gzip
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import norm
from sklearn.mixture import GaussianMixture


if len(sys.argv) != 3:
    sys.stderr.write(
        "Usage: plot_bq_gmm_from_percell_table.py <per_cell_table.tsv.gz> <out_prefix>\n"
    )
    sys.exit(1)

infile = sys.argv[1]
out_prefix = sys.argv[2]

plot_path = out_prefix + ".GMM_meanBQ_density.png"
qc_path = out_prefix + ".GMM_meanBQ_QC.tsv.gz"


bases = ("A", "C", "G", "T")

# Aggregation:
# key = (pos, allele)
# value = {
#   total_count: total reads supporting this allele at this position across all cells
#   total_bq_sum: total sum of base qualities reconstructed as count * mean_bq for each cell
# }
agg = defaultdict(lambda: {"total_count": 0, "total_bq_sum": 0.0})

# Read per-cell table and aggregate into per-position, per-allele mean BQ
with gzip.open(infile, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")

    required_cols = {
        "barcode", "pos", "A", "C", "G", "T", "depth",
        "A_mean_bq", "C_mean_bq", "G_mean_bq", "T_mean_bq"
    }
    missing = required_cols - set(reader.fieldnames)
    if missing:
        sys.stderr.write(f"Missing required columns: {', '.join(sorted(missing))}\n")
        sys.exit(1)

    for row in reader:
        pos = int(row["pos"])

        for b in bases:
            cnt = int(row[b])
            if cnt == 0:
                continue

            mean_bq_str = row[f"{b}_mean_bq"]
            if mean_bq_str == "NA" or mean_bq_str == "":
                continue

            mean_bq = float(mean_bq_str)

            agg[(pos, b)]["total_count"] += cnt
            agg[(pos, b)]["total_bq_sum"] += cnt * mean_bq


# Build table of per-position, per-allele mean BQ values
records = []
for (pos, allele), d in agg.items():
    if d["total_count"] <= 0:
        continue
    mean_bq = d["total_bq_sum"] / d["total_count"]
    records.append((pos, allele, d["total_count"], mean_bq))

if len(records) < 10:
    sys.stderr.write(
        "Too few per-position/per-allele records to fit a 3-component GMM.\n"
    )
    sys.exit(1)

values = np.array([r[3] for r in records], dtype=float).reshape(-1, 1)

# Fit 3-component Gaussian mixture model
gmm = GaussianMixture(
    n_components=3,
    covariance_type="full",
    random_state=0,
    n_init=20,
    reg_covar=1e-4
)
gmm.fit(values)

means = gmm.means_.flatten()
weights = gmm.weights_.flatten()
covars = gmm.covariances_.reshape(-1)
sds = np.sqrt(covars)

# Order components by mean
order = np.argsort(means)
low_comp = order[0]
mid_comp = order[1]
high_comp = order[2]   # highest mean = high-confidence component

# Posterior probability that each observed value belongs to high-confidence component
obs_probs = gmm.predict_proba(values)[:, high_comp]

# Save QC table
with gzip.open(qc_path, "wt") as out:
    out.write("pos\tallele\tn_reads\tmean_bq\tp_high_conf\tpass_99pct\n")
    for rec, p in zip(records, obs_probs):
        pos, allele, n_reads, mean_bq = rec
        passed = "YES" if p > 0.99 else "NO"
        out.write(f"{pos}\t{allele}\t{n_reads}\t{mean_bq:.6f}\t{p:.6f}\t{passed}\n")

# Build density curves
x_min = max(0.0, float(np.floor(values.min())) - 2.0)
x_max = float(np.ceil(values.max())) + 2.0
x = np.linspace(x_min, x_max, 2000).reshape(-1, 1)

total_density = np.exp(gmm.score_samples(x))
grid_probs = gmm.predict_proba(x)[:, high_comp]

# Threshold = smallest x with P(high-confidence) > 0.99
threshold = None
idx = np.where(grid_probs > 0.99)[0]
if len(idx) > 0:
    threshold = float(x[idx[0], 0])

component_densities = {}
for k in range(3):
    component_densities[k] = weights[k] * norm.pdf(
        x[:, 0], loc=means[k], scale=sds[k]
    )

# Plot
plt.figure(figsize=(10, 6))

# Histogram as density
plt.hist(
    values[:, 0],
    bins=40,
    density=True,
    alpha=0.35,
    edgecolor="black",
    label="Observed distribution"
)

component_colors = {
    low_comp: "red",
    mid_comp: "green",
    high_comp: "blue"
}

component_labels = {
    low_comp: f"Low-BQ component (mean={means[low_comp]:.2f})",
    mid_comp: f"Intermediate-BQ component (mean={means[mid_comp]:.2f})",
    high_comp: f"High-confidence component (mean={means[high_comp]:.2f})"
}

for k in order:
    plt.plot(
        x[:, 0],
        component_densities[k],
        color=component_colors[k],
        linewidth=2,
        label=component_labels[k]
    )

plt.plot(
    x[:, 0],
    total_density,
    color="black",
    linewidth=2,
    linestyle="-",
    label="Total mixture density"
)

if threshold is not None:
    plt.axvline(
        threshold,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label=f"99% threshold = {threshold:.2f}"
    )

plt.xlabel("Per-position, per-allele mean base quality")
plt.ylabel("Density")
plt.title("Gaussian mixture model fit over per-position, per-allele mean base qualities")
plt.legend(fontsize=8)
plt.tight_layout()

# Optional visual cap, because you asked for Y from 0 to 1.
# This is only a display choice; density is not a probability.
plt.ylim(0, 1)

plt.savefig(plot_path, dpi=300)
plt.close()

sys.stderr.write(f"Saved plot: {plot_path}\n")
sys.stderr.write(f"Saved QC table: {qc_path}\n")

if threshold is not None:
    sys.stderr.write(f"Estimated 99% high-confidence threshold: {threshold:.4f}\n")
else:
    sys.stderr.write("Could not determine a 99% high-confidence threshold.\n")

