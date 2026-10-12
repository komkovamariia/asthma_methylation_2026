#!/usr/bin/env python3
"""Prepare lane pairs for the historical read-merging trial.

Writes merge_manifest.tsv and a header-only samplesheet placeholder.
This script does not merge reads or populate the final samplesheet.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

R1_PATTERNS = [
    re.compile(r"(.*)_R1(?:_001)?\.f(ast)?q\.gz$", re.IGNORECASE),
    re.compile(r"(.*)_1\.f(ast)?q\.gz$", re.IGNORECASE),
]
R2_REPLACE = [("_R1", "_R2"), ("_R1_001", "_R2_001"), ("_1.", "_2.")]
LANE_RE = re.compile(r"_L00(\d{1,2})_")
SAMPLE_RE = re.compile(
    r"^(?P<sample>.+?)(_S\d+_L00\d+_|_L00\d+_|_R1|_R2|_1\.|_2\.)", re.IGNORECASE
)


def infer_lane(name: str) -> str:
    m = LANE_RE.search(name)
    if m:
        return f"L{int(m.group(1)):03d}"
    return "L001"


def infer_sample(stem: str) -> str:
    m = SAMPLE_RE.search(stem)
    if m:
        return m.group("sample")
    return re.sub(r"(_R[12].*)$", "", stem, flags=re.IGNORECASE)


def find_r2_for_r1(r1_path: Path) -> Path | None:
    s = r1_path.name
    candidates = set()
    for a, b in R2_REPLACE:
        if a in s:
            candidates.add(s.replace(a, b, 1))
    for cand in candidates:
        p = r1_path.with_name(cand)
        if p.exists():
            return p
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq_dir", required=True)
    ap.add_argument("--merged_dir", required=True)
    ap.add_argument("--out_manifest", required=True)
    ap.add_argument("--out_samplesheet", required=True)
    args = ap.parse_args()

    fastq_dir = Path(args.fastq_dir).resolve()
    merged_dir = Path(args.merged_dir).resolve()
    merged_dir.mkdir(parents=True, exist_ok=True)

    r1_files = []
    for p in fastq_dir.rglob("*.f*q.gz"):
        if re.search(
            r"_R1(?:_001)?\.f(ast)?q\.gz$", p.name, re.IGNORECASE
        ) or re.search(r"_1\.f(ast)?q\.gz$", p.name, re.IGNORECASE):
            r1_files.append(p)
    if not r1_files:
        print(f"ERROR: No R1 fastq.gz files found under: {fastq_dir}", file=sys.stderr)
        sys.exit(2)

    rows = []
    missing = 0
    for r1 in sorted(r1_files):
        r2 = find_r2_for_r1(r1)
        if not r2:
            missing += 1
            continue
        lane = infer_lane(r1.name)
        stem = r1.name
        sample = infer_sample(stem)
        merged_lane = merged_dir / "lanes" / sample / f"{sample}_{lane}.merged.fq.gz"
        merged_lane.parent.mkdir(parents=True, exist_ok=True)
        rows.append((sample, lane, str(r1), str(r2), str(merged_lane)))

    if missing:
        print(
            f"WARNING: {missing} R1 files had no matching R2; they were skipped.",
            file=sys.stderr,
        )

    if not rows:
        print("ERROR: Could not pair any R1/R2 files.", file=sys.stderr)
        sys.exit(3)

    out_manifest = Path(args.out_manifest)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with out_manifest.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "lane", "r1", "r2", "merged_lane_fastq"])
        w.writerows(rows)

    out_samplesheet = Path(args.out_samplesheet)
    out_samplesheet.parent.mkdir(parents=True, exist_ok=True)
    with out_samplesheet.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "fastq_1", "fastq_2"])

    samples = sorted({r[0] for r in rows})
    print(
        f"OK: manifest={out_manifest}  samplesheet={out_samplesheet} (will be populated after concat)  samples={len(samples)}  lane_pairs={len(rows)}"
    )


if __name__ == "__main__":
    main()
