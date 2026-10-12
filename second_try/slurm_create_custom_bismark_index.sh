#!/usr/bin/env bash
#SBATCH --job-name=custom_bismark_index
#SBATCH --output=logs/custom_bismark_index_%j.out
#SBATCH --error=logs/custom_bismark_index_%j.err
#SBATCH --partition=infinite
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference/tz_amplicons}"
REFERENCE_FASTA="${REFERENCE_FASTA:-${REFERENCE_DIR}/tz_targets.fasta}"
BISMARK_INDEX_DIR="${BISMARK_INDEX_DIR:-${REFERENCE_DIR}}"
mkdir -p "${PROJECT_DIR}/logs"
[[ -s "${REFERENCE_FASTA}" ]] || { echo "Run create_tz_reference.py first: ${REFERENCE_FASTA}" >&2; exit 1; }

# Bismark writes its index beside the FASTA, so these directories must agree.
REF_FASTA_DIR="$(cd -- "$(dirname -- "${REFERENCE_FASTA}")" && pwd)"
INDEX_DIR="$(cd -- "${BISMARK_INDEX_DIR}" && pwd)"
if [[ "${REF_FASTA_DIR}" != "${INDEX_DIR}" ]]; then
  echo "BISMARK_INDEX_DIR must be the directory containing REFERENCE_FASTA" >&2
  exit 1
fi

index_complete() {
  local conversion count
  for conversion in CT GA; do
    count=$(find "${INDEX_DIR}/Bisulfite_Genome/${conversion}_conversion" \
      -maxdepth 1 -type f -name '*.bt2*' -size +0c 2>/dev/null | wc -l) || return 1
    [[ "${count}" -ge 6 ]] || return 1
  done
}

if index_complete; then
  echo "Bismark index already exists: ${INDEX_DIR}/Bisulfite_Genome"
  exit 0
fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
bismark_genome_preparation --verbose "${REF_FASTA_DIR}"
index_complete || { echo "Bismark index is incomplete" >&2; exit 1; }
echo "Bismark index ready: ${INDEX_DIR}/Bisulfite_Genome"
