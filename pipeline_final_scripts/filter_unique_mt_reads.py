#!/usr/bin/env python3
import sys
import pysam

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: filter_unique_mt_reads.py <input.bam> <mt_chr> <output.bam> <mapq_min>\n"
    )
    sys.exit(1)

in_bam = sys.argv[1]
mt_chr = sys.argv[2]
out_bam = sys.argv[3]
mapq_min = int(sys.argv[4])

bam = pysam.AlignmentFile(in_bam, "rb")
out = pysam.AlignmentFile(out_bam, "wb", template=bam)

n_total = 0
n_pass = 0

for read in bam.fetch(until_eof=True):
    n_total += 1

    if read.is_unmapped:
        continue
    if read.is_secondary or read.is_supplementary:
        continue
    if read.reference_name != mt_chr:
        continue
    if read.mapping_quality < mapq_min:
        continue

    # strict: must be uniquely mapped according to aligner
    if not read.has_tag("NH"):
        continue
    if read.get_tag("NH") != 1:
        continue

    out.write(read)
    n_pass += 1

bam.close()
out.close()

sys.stderr.write(f"total_reads_seen={n_total}\n")
sys.stderr.write(f"unique_mt_reads_written={n_pass}\n")
