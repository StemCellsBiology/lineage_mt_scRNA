#!/usr/bin/env python3
import sys
import pandas as pd

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: select_barcodes_by_mt_reads.py <reads_per_cell.tsv> <min_reads> <out_barcodes.tsv>\n"
    )
    sys.exit(1)

in_file = sys.argv[1]
min_reads = int(sys.argv[2])
out_file = sys.argv[3]

df = pd.read_csv(in_file, sep="\t")

keep = df[df["unique_mt_reads"] >= min_reads]

with open(out_file, "w") as out:
    for bc in keep["barcode"]:
        out.write(f"{bc}\n")

sys.stderr.write(f"barcodes_before={df.shape[0]}\n")
sys.stderr.write(f"barcodes_after={keep.shape[0]}\n")
