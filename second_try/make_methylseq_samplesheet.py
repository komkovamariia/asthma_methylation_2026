#!/usr/bin/env python3
"""Build a paired-end samplesheet for hla1 through hla66, retaining lane rows."""

import argparse
import csv
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq_dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    fastq_dir = args.fastq_dir.expanduser().resolve()
    if not fastq_dir.is_dir():
        raise NotADirectoryError(fastq_dir)

    rows = []
    for r1 in sorted(fastq_dir.rglob("*_R1_*.fastq.gz")):
        sample = re.sub(r"_S[0-9]+_L[0-9]+_R1_.*", "", r1.name)
        if re.fullmatch(r"hla([1-9]|[1-5][0-9]|6[0-6])", sample) is None:
            continue
        r2 = r1.with_name(r1.name.replace("_R1_", "_R2_", 1))
        for fastq in (r1, r2):
            if not fastq.is_file() or fastq.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty FASTQ: {fastq}")
        rows.append((sample, str(r1), str(r2)))
    if not rows:
        raise ValueError(
            f"No hla1 through hla66 paired-end FASTQs found in {fastq_dir}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("sample", "fastq_1", "fastq_2"))
        writer.writerows(rows)
    print(
        f"Wrote {len(rows)} lane rows for {len({row[0] for row in rows})} samples: {args.out}"
    )


if __name__ == "__main__":
    main()
