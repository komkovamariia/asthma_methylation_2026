#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR to the second_try_after directory}"
: "${INPUT_FASTQ_DIR:?Set INPUT_FASTQ_DIR}"
: "${REF_DIR:=/projects/mice_transplant_2025/reference}"
: "${REF_FASTA:=${REF_DIR}/hg38.fa}"
: "${ALIGN_DIR:=${PROJECT_DIR}/hg38_realign_out}"
: "${TMP_DIR:=${PROJECT_DIR}/hg38_realign_tmp}"
: "${ARRAY_CONCURRENCY:=10}"
: "${PARTITION:=long}"
: "${ONLY_SAMPLE_IDS:=}"
: "${ONLY_SAMPLES_FILE:=}"
: "${SKIP_SAMPLE_IDS:=Undetermined,hla4,hla5}"

[[ -s "${REF_FASTA}" && -d "${REF_DIR}/Bisulfite_Genome" ]] || {
  echo "Missing hg38 FASTA or Bismark index in ${REF_DIR}" >&2
  exit 2
}
mkdir -p "${ALIGN_DIR}" "${TMP_DIR}" "${PROJECT_DIR}/logs"
SAMPLESHEET_TSV="${PROJECT_DIR}/samplesheet_realign_hg38.tsv"

# A separate array row per lane avoids mixing input libraries.
python3 - "${INPUT_FASTQ_DIR}" "${SAMPLESHEET_TSV}" "${SKIP_SAMPLE_IDS}" \
  "${ONLY_SAMPLE_IDS}" "${ONLY_SAMPLES_FILE}" <<'TSV'
import csv
import re
import sys
from pathlib import Path

input_dir, output, skip_ids, only_ids, only_file = sys.argv[1:]
input_dir = Path(input_dir)
if not input_dir.is_dir():
    raise NotADirectoryError(input_dir)
skip = {s.strip() for s in skip_ids.split(",") if s.strip()}
only = {s.strip() for s in only_ids.split(",") if s.strip()}
if only_file:
    only.update(s.strip() for s in Path(only_file).read_text().splitlines() if s.strip())
r1_files = sorted(input_dir.rglob("*_R1_*.fastq.gz"))
rows = []
for r1 in r1_files or sorted(input_dir.glob("hla*.merged.fq.gz")):
    if r1_files:
        sample = re.sub(r"_S[0-9]+_L[0-9]+_R1_.*", "", r1.name)
        r2 = r1.with_name(r1.name.replace("_R1_", "_R2_", 1))
    else:
        sample = r1.name.split(".", 1)[0].split("_merged", 1)[0]
        r2 = None
    if sample in skip or (only and sample not in only):
        continue
    for fastq in (r1, r2):
        if fastq is not None and (not fastq.is_file() or fastq.stat().st_size == 0):
            raise FileNotFoundError(f"Missing or empty FASTQ: {fastq}")
    rows.append((sample, str(r1.resolve()), str(r2.resolve()) if r2 else ""))
if not rows:
    raise ValueError("No input rows remain after filtering")
with open(output, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(("sample", "fastq_1", "fastq_2"))
    writer.writerows(rows)
print(f"Wrote {len(rows)} input rows to {output}")
TSV

COUNT=$(($(wc -l < "${SAMPLESHEET_TSV}") - 1))
export SAMPLESHEET_TSV ALIGN_DIR TMP_DIR REF_DIR
cd "${PROJECT_DIR}"
sbatch --partition="${PARTITION}" --array="1-${COUNT}%${ARRAY_CONCURRENCY}" \
  "${PROJECT_DIR}/slurm_realign_hg38.sh"
