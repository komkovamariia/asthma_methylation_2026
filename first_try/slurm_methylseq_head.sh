#!/bin/bash
#SBATCH --job-name=nf_methylseq
#SBATCH --output=nf_methylseq_%j.out
#SBATCH --error=nf_methylseq_%j.err
#SBATCH --partition=infinite
#SBATCH --time=20-00:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G

set -euo pipefail
export PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}}"
export RAW_FASTQ_DIR="${RAW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/2025_12_kofiadi_seq}"
export NXF_HOME="${NXF_HOME:-$PROJECT_DIR/.nextflow}"
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-asthma_methylseq}"
cd "$PROJECT_DIR"
exec bash "$PROJECT_DIR/run_methylseq.sh"
