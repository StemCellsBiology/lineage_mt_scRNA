#!/usr/bin/env python3

import sys
import gzip
import pysam
import matplotlib

matplotlib.use("Agg")  # ważne na serwerze bez GUI

import matplotlib.pyplot as plt
from collections import defaultdict, Counter


if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: per_cell_per_allele_counts.py <unique_mt.bam> <barcodes.tsv> <out.tsv.gz>\n"
    )
    sys.exit(1)


bam_path = sys.argv[1]
barcodes_path = sys.argv[2]
out_path = sys.argv[3]

# Automatyczna nazwa wykresu
if out_path.endswith(".tsv.gz"):
    plot_path = out_path.replace(".tsv.gz", ".base_quality_distribution.png")
elif out_path.endswith(".gz"):
    plot_path = out_path.replace(".gz", ".base_quality_distribution.png")
else:
    plot_path = out_path + ".base_quality_distribution.png"


valid = set()
with open(barcodes_path) as f:
    for line in f:
        bc = line.strip()
        if bc:
            valid.add(bc)


bases = ("A", "C", "G", "T")

# Globalny rozkład jakości baz
# Zliczamy wszystkie jakości baz, które przechodzą filtry:
# mapped, primary alignment, CB tag, valid barcode, A/C/G/T
bq_hist = Counter()

bam = pysam.AlignmentFile(bam_path, "rb")


with gzip.open(out_path, "wt") as out:
    out.write(
        "barcode\tpos\tA\tC\tG\tT\tdepth\t"
        "A_mean_bq\tC_mean_bq\tG_mean_bq\tT_mean_bq\n"
    )

    for pileupcolumn in bam.pileup(
        truncate=True,
        stepper="samtools",
        min_base_quality=0,      # bardzo ważne: nie filtrujemy po jakości bazy
        ignore_overlaps=False,
        ignore_orphans=False
    ):
        pos = pileupcolumn.reference_pos + 1

        per_cell = defaultdict(lambda: {
            "A": 0, "C": 0, "G": 0, "T": 0,
            "A_bq": 0, "C_bq": 0, "G_bq": 0, "T_bq": 0
        })

        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue

            read = pileupread.alignment

            if read.is_unmapped:
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            if not read.has_tag("CB"):
                continue

            cb = read.get_tag("CB")
            if cb not in valid:
                continue

            qpos = pileupread.query_position
            if qpos is None:
                continue

            base = read.query_sequence[qpos].upper()
            if base not in bases:
                continue

            bq = read.query_qualities[qpos]

            # Nie filtrujemy po base quality.
            # Zapisujemy każdą jakość do globalnego histogramu.
            bq_hist[bq] += 1

            per_cell[cb][base] += 1
            per_cell[cb][f"{base}_bq"] += bq

        for cb, cnt in per_cell.items():
            depth = cnt["A"] + cnt["C"] + cnt["G"] + cnt["T"]
            if depth == 0:
                continue

            mean_bq = {}
            for b in bases:
                mean_bq[b] = cnt[f"{b}_bq"] / cnt[b] if cnt[b] > 0 else "NA"

            out.write(
                f"{cb}\t{pos}\t"
                f"{cnt['A']}\t{cnt['C']}\t{cnt['G']}\t{cnt['T']}\t{depth}\t"
                f"{mean_bq['A']}\t{mean_bq['C']}\t{mean_bq['G']}\t{mean_bq['T']}\n"
            )


bam.close()


# Wykres rozkładu jakości baz
if len(bq_hist) > 0:
    qualities = sorted(bq_hist.keys())
    counts = [bq_hist[q] for q in qualities]

    plt.figure(figsize=(10, 6))
    plt.bar(qualities, counts)
    plt.xlabel("Base quality score")
    plt.ylabel("Number of observed bases")
    plt.title("Distribution of base quality scores")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()

    sys.stderr.write(f"Saved base quality distribution plot: {plot_path}\n")
else:
    sys.stderr.write("No base qualities collected. Plot was not generated.\n")

