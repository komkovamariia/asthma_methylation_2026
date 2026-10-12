#!/usr/bin/env python3
"""Map hla1 through hla66 to patient, cell type and timepoint.

Each patient has CD4 before, CD8 before, CD4 after and CD8 after in that order.
Samples 61 through 66 are CD4/CD8 pairs from three controls.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

CONTROL_MAPPING: dict[int, dict[str, str]] = {
    61: {"patient_id": "control_1", "cell_type": "CD4+"},
    62: {"patient_id": "control_1", "cell_type": "CD8+"},
    63: {"patient_id": "control_2", "cell_type": "CD4+"},
    64: {"patient_id": "control_2", "cell_type": "CD8+"},
    65: {"patient_id": "control_3", "cell_type": "CD4+"},
    66: {"patient_id": "control_3", "cell_type": "CD8+"},
}


def parse_sample_number(sample_id: str) -> int:
    """Extract the numeric hla sample identifier."""
    m = re.search(r"(?i)\bhla(?P<num>\d+)\b", sample_id)
    if not m:
        raise ValueError(
            f"Cannot extract number from sample_id='{sample_id}' (expected 'hlaN')"
        )
    return int(m.group("num"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--samplesheet_csv",
        required=True,
        type=Path,
        help="methylseq samplesheet.csv (with 'sample' column).",
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
        help="Output TSV path (sample_sheet.tsv).",
    )
    args = ap.parse_args()

    samplesheet_csv: Path = args.samplesheet_csv.expanduser().resolve()
    out_tsv: Path = args.out_tsv.expanduser().resolve()

    df = pd.read_csv(samplesheet_csv)
    if "sample" not in df.columns:
        raise ValueError("samplesheet_csv must contain 'sample' column")

    unique_samples = sorted(set(df["sample"].astype(str)))
    rows: list[dict[str, object]] = []

    for sample_id in unique_samples:
        n = parse_sample_number(sample_id)
        if not (1 <= n <= 66):
            continue

        if 1 <= n <= 60:
            patient_num = ((n - 1) // 4) + 1
            within_patient = ((n - 1) % 4) + 1
            if within_patient == 1:
                cell_type, timepoint = "CD4+", "before"
            elif within_patient == 2:
                cell_type, timepoint = "CD8+", "before"
            elif within_patient == 3:
                cell_type, timepoint = "CD4+", "after"
            else:
                cell_type, timepoint = "CD8+", "after"

            rows.append(
                {
                    "sample_id": sample_id,
                    "number": n,
                    "patient_id": f"patient_{patient_num}",
                    "cell_type": cell_type,
                    "timepoint": timepoint,
                    "group": "patient",
                }
            )
        else:
            info = CONTROL_MAPPING.get(n)
            if not info:
                continue
            rows.append(
                {
                    "sample_id": sample_id,
                    "number": n,
                    "patient_id": info["patient_id"],
                    "cell_type": info["cell_type"],
                    "timepoint": "",
                    "group": "control",
                }
            )

    if not rows:
        raise ValueError("No hla1 through hla66 samples found in samplesheet")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(rows).sort_values("number").reset_index(drop=True)
    out_df.to_csv(out_tsv, sep="\t", index=False)

    print(f"Wrote metadata TSV: {out_tsv} (rows={len(out_df)})")
    print(out_df.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
