#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-${SCRIPT_DIR}}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/../second_try/methylseq_out}"
METADATA_TSV="${METADATA_TSV:-${PROJECT_DIR}/sample_sheet.tsv}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/methylseq_postprocess}"
mkdir -p "${RESULTS_DIR}"

exec Rscript "${PROJECT_DIR}/postprocess_dss.R" \
  --methylseq_outdir "${OUTDIR}" \
  --metadata_tsv "${METADATA_TSV}" \
  --results_dir "${RESULTS_DIR}" \
  --min_cov "${MIN_COV:-5}" \
  --p_threshold "${Q_THRESHOLD:-0.05}" \
  --delta "${DELTA:-0.05}" \
  --top_n_heatmap "${TOP_N_HEATMAP:-500}"
