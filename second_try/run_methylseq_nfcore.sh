#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR to the second_try directory}"
: "${RESULTS_DIR:?Set RESULTS_DIR}"
: "${WORK_DIR:?Set WORK_DIR}"
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA}"
: "${BISMARK_INDEX_DIR:?Set BISMARK_INDEX_DIR}"
: "${NON_DIRECTIONAL:=1}"
: "${TARGET_BED:=}"
: "${NFCOREFLOW:=nf-core/methylseq}"
: "${NFCORE_VERSION:=4.2.0}"
: "${PROFILE:=conda}"
: "${RESUME:=1}"
: "${SKIP_SAMPLE_IDS:=hla4,hla5}"

mkdir -p "${PROJECT_DIR}/logs" "${RESULTS_DIR}" "${WORK_DIR}"
SAMPLESHEET="${PROJECT_DIR}/samplesheet_methylseq.csv"
SAMPLESHEET_FILTERED="${PROJECT_DIR}/samplesheet_methylseq.filtered.csv"

# CSV parsing also supports commas in FASTQ paths.
python3 - "${SAMPLESHEET}" "${SAMPLESHEET_FILTERED}" "${SKIP_SAMPLE_IDS}" <<'CSV'
import csv
import sys
from pathlib import Path

source, destination, skip_ids = sys.argv[1:]
skip = {s.strip() for s in skip_ids.split(",") if s.strip()} | {"Undetermined"}
with open(source, newline="") as handle:
    reader = csv.DictReader(handle)
    if not {"sample", "fastq_1", "fastq_2"}.issubset(reader.fieldnames or []):
        raise ValueError("Samplesheet requires sample, fastq_1 and fastq_2 columns")
    fields = reader.fieldnames
    rows = [row for row in reader if row["sample"] and row["sample"] not in skip]
if not rows:
    raise ValueError("No samples remain after filtering")
for row in rows:
    for field in ("fastq_1", "fastq_2"):
        fastq = row[field]
        if field == "fastq_2" and not fastq:
            continue
        if not fastq or not Path(fastq).is_file() or Path(fastq).stat().st_size == 0:
            raise FileNotFoundError(f"{row['sample']}: missing or empty {field}: {fastq}")
with open(destination, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
print(f"Validated {len(rows)} lane rows; excluded: {', '.join(sorted(skip))}")
CSV

[[ -s "${REFERENCE_FASTA}" ]] || { echo "Missing reference: ${REFERENCE_FASTA}" >&2; exit 1; }
[[ -d "${BISMARK_INDEX_DIR}/Bisulfite_Genome" ]] || { echo "Missing Bismark index: ${BISMARK_INDEX_DIR}" >&2; exit 1; }
[[ -f "${PROJECT_DIR}/aldan3_methylseq.config" ]] || { echo "Missing ALDAN3 config" >&2; exit 1; }
if [[ -n "${TARGET_BED}" && ! -s "${TARGET_BED}" ]]; then
  echo "Missing target BED: ${TARGET_BED}" >&2
  exit 1
fi

cmd=(nextflow run "${NFCOREFLOW}" -r "${NFCORE_VERSION}"
  -profile "${PROFILE}"
  -work-dir "${WORK_DIR}"
  -c "${PROJECT_DIR}/aldan3_methylseq.config"
  --input "${SAMPLESHEET_FILTERED}"
  --outdir "${RESULTS_DIR}"
  --fasta "${REFERENCE_FASTA}"
  --bismark_index "${BISMARK_INDEX_DIR}"
  --skip_deduplication
  --save_reference
  --cytosine_report)

if [[ "${NON_DIRECTIONAL}" == "1" ]]; then
  cmd+=(--non_directional)
fi
if [[ -n "${TARGET_BED}" ]]; then
  cmd+=(--run_targeted_sequencing --target_regions_file "${TARGET_BED}"
    --run_qualimap --bamqc_regions_file "${TARGET_BED}")
fi
if [[ "${RESUME}" == "1" ]]; then
  cmd+=(-resume)
fi

printf 'Running:'
printf ' %q' "${cmd[@]}"
printf '\n'
# Preserve Nextflow's exit status for SLURM.
exec "${cmd[@]}"
