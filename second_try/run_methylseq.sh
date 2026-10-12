#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR to the second_try directory}"
: "${RAW_FASTQ_DIR:?Set RAW_FASTQ_DIR to the paired-end FASTQ directory}"

mkdir -p "${PROJECT_DIR}/logs"
python3 "${PROJECT_DIR}/make_methylseq_samplesheet.py" \
  --fastq_dir "${RAW_FASTQ_DIR}" \
  --out "${PROJECT_DIR}/samplesheet_methylseq.csv"
