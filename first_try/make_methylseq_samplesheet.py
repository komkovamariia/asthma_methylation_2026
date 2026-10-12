#!/usr/bin/env python3
"""Build the historical four-column hg38 samplesheet from flat Illumina FASTQs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

IGNORE_PREFIXES = {"undetermined", "methylom1", "methylom2", "methylom3"}


def _extract_sample_name(path: Path) -> str | None:
    token = path.name.split("_", 1)[0].strip()
    if not token:
        return None
    if token.lower() in IGNORE_PREFIXES:
        return None
    return token


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--fastq_dir",
        required=True,
        type=Path,
        help="Directory containing gzipped FASTQ files (Illumina).",
    )
    ap.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output CSV path (samplesheet.csv).",
    )
    ap.add_argument(
        "--genome",
        default="",
        help=(
            "Value to put into optional 'genome' column. "
            "Use a short reference key (no slashes). "
            "and only works if configured via nf-core configs / iGenomes. "
            "If unsure, leave empty (default) and provide --fasta/--bismark_index when running the pipeline."
        ),
    )
    ap.add_argument(
        "--sample_regex",
        default=r"^(?P<prefix>hla)(?P<num>\d+)$",
        help="Regex to identify sample IDs to include (default: ^(hla)(\\d+)$).",
    )
    ap.add_argument(
        "--min_sample_num",
        type=int,
        default=1,
        help="Minimum sample number to include (default: 1).",
    )
    ap.add_argument(
        "--max_sample_num",
        type=int,
        default=66,
        help="Maximum sample number to include (default: 66).",
    )

    args = ap.parse_args()

    fastq_dir: Path = args.fastq_dir.expanduser().resolve()
    out: Path = args.out.expanduser().resolve()

    if not fastq_dir.exists():
        raise FileNotFoundError(f"--fastq_dir does not exist: {fastq_dir}")

    rx = re.compile(args.sample_regex, flags=re.IGNORECASE)

    r1_files = sorted(fastq_dir.glob("*_R1_001.fastq.gz"))
    rows: list[dict[str, str]] = []

    for r1 in r1_files:
        sample = _extract_sample_name(r1)
        if sample is None:
            continue

        m = rx.match(sample)
        if not m:
            continue

        num = int(m.group("num"))
        if not (args.min_sample_num <= num <= args.max_sample_num):
            continue

        r2 = Path(str(r1).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.exists():
            raise FileNotFoundError(
                f"Missing R2 pair for {r1.name}: expected {r2.name}"
            )

        rows.append(
            {
                "sample": sample,
                "fastq_1": str(r1),
                "fastq_2": str(r2),
                "genome": args.genome or "",
            }
        )

    if not rows:
        raise RuntimeError(
            f"No FASTQs matched in {fastq_dir}. "
            "Expected pattern '*_R1_001.fastq.gz' with sample like 'hla1'."
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["sample", "fastq_1", "fastq_2", "genome"]
        )
        writer.writeheader()
        writer.writerows(rows)

    sample_counts = Counter(r["sample"] for r in rows)
    sample_counts_sorted = dict(sorted(sample_counts.items(), key=lambda kv: kv[0]))
    print(f"Wrote samplesheet: {out}")
    print(f"Rows: {len(rows)}  |  Unique samples: {len(sample_counts_sorted)}")
    print("Example sample lane counts:")
    for i, (sid, cnt) in enumerate(sample_counts_sorted.items()):
        if i >= 10:
            break
        print(f"{sid}\t{cnt}")


if __name__ == "__main__":
    main()
