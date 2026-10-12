#!/usr/bin/env bash
#SBATCH --job-name=methylseq_nfcore
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=infinite

set -euo pipefail

export PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}}"
export RAW_FASTQ_DIR="${RAW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/2025_12_kofiadi_seq}"
export RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/methylseq_out}"
export WORK_DIR="${WORK_DIR:-${PROJECT_DIR}/methylseq_work}"
export REFERENCE_FASTA="${REFERENCE_FASTA:-${PROJECT_DIR}/reference/tz_amplicons/tz_targets.fasta}"
export BISMARK_INDEX_DIR="${BISMARK_INDEX_DIR:-${PROJECT_DIR}/reference/tz_amplicons}"
export NON_DIRECTIONAL="${NON_DIRECTIONAL:-1}"
export NFCOREFLOW="${NFCOREFLOW:-nf-core/methylseq}"
export NFCORE_VERSION="${NFCORE_VERSION:-4.2.0}"
export PROFILE="${PROFILE:-conda}"
export RESUME="${RESUME:-1}"
export NXF_HOME="${NXF_HOME:-${PROJECT_DIR}/.nextflow}"

cd "${PROJECT_DIR}"
mkdir -p logs
printf 'Job %s on %s\n' "${SLURM_JOB_ID:-local}" "${SLURM_NODELIST:-localhost}"

exec bash "${PROJECT_DIR}/run_methylseq_nfcore.sh"
