#!/usr/bin/env python3
import sys
import pysam
import pandas as pd
from collections import Counter

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: count_unique_mt_reads_per_cell.py <unique_mt.bam> <barcodes.tsv> <out.tsv>\n"
    )
    sys.exit(1)

bam_path = sys.argv[1]
barcodes_path = sys.argv[2]
out_path = sys.argv[3]

valid = set()
with open(barcodes_path) as f:
    for line in f:
        bc = line.strip()
        if bc:
            valid.add(bc)

counts = Counter()

bam = pysam.AlignmentFile(bam_path, "rb")

for read in bam.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    if read.is_secondary or read.is_supplementary:
        continue
    if not read.has_tag("CB"):
        continue

    cb = read.get_tag("CB")
    if cb in valid:
        counts[cb] += 1

bam.close()

rows = []
for bc in sorted(valid):
    rows.append({
        "barcode": bc,
        "unique_mt_reads": counts.get(bc, 0)
    })

pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
