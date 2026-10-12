#!/usr/bin/env bash
set -euo pipefail

RAW_FASTQ_DIR="${RAW_FASTQ_DIR:-/projects/mice_transplant_2025/raw/2025_12_kofiadi_seq}"

PROJECT_DIR="${PROJECT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

OUTDIR="${OUTDIR:-$PROJECT_DIR/methylseq_out}"

WORKDIR="${WORKDIR:-$PROJECT_DIR/methylseq_work}"

GENOME="${GENOME:-}"
FASTA="${FASTA:-/projects/mice_transplant_2025/reference/hg38.fa}"
BISMARK_INDEX="${BISMARK_INDEX:-/projects/mice_transplant_2025/reference}"

PROFILE="${PROFILE:-conda}"

RESUME="${RESUME:-true}"  # true/false

SKIP_DEDUPLICATION="${SKIP_DEDUPLICATION:-true}"  # true/false

RUN_TARGETED_SEQUENCING="${RUN_TARGETED_SEQUENCING:-false}"  # true/false
TARGET_BED="${TARGET_BED:-$PROJECT_DIR/targets/targets.bed}"  # used only if RUN_TARGETED_SEQUENCING=true

SAMPLESHEET="$PROJECT_DIR/samplesheet_methylseq.csv"
METADATA_TSV="$PROJECT_DIR/sample_sheet.tsv"
ALDAN3_CONFIG="$PROJECT_DIR/aldan3_methylseq.config"

mkdir -p "$PROJECT_DIR" "$OUTDIR" "$WORKDIR"

python3 "$PROJECT_DIR/make_methylseq_samplesheet.py" \
  --fastq_dir "$RAW_FASTQ_DIR" \
  --out "$SAMPLESHEET"

python3 "$PROJECT_DIR/make_tz_metadata.py" \
  --samplesheet_csv "$SAMPLESHEET" \
  --out_tsv "$METADATA_TSV"

EXTRA_ARGS=()
if [[ "$RESUME" == "true" || "$RESUME" == "1" ]]; then
  EXTRA_ARGS+=( -resume )
fi
if [[ "$SKIP_DEDUPLICATION" == "true" ]]; then
  EXTRA_ARGS+=( --skip_deduplication )
fi

if [[ "$RUN_TARGETED_SEQUENCING" == "true" ]]; then
  EXTRA_ARGS+=( --run_targeted_sequencing --target_regions_file "$TARGET_BED" )
fi

REF_ARGS=()
if [[ -n "${GENOME}" ]]; then
  REF_ARGS+=( --genome "$GENOME" )
else
  if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA not found: $FASTA" >&2
    echo "Set FASTA=/path/to/hg38.fa" >&2
    exit 2
  fi
  if [[ ! -d "$BISMARK_INDEX/Bisulfite_Genome" ]]; then
    echo "ERROR: Bismark index not found under: $BISMARK_INDEX/Bisulfite_Genome" >&2
    echo "Set BISMARK_INDEX=/path/to/reference_dir (that contains Bisulfite_Genome/)" >&2
    exit 2
  fi
  REF_ARGS+=( --fasta "$FASTA" --bismark_index "$BISMARK_INDEX" )
fi

set -x
nextflow run nf-core/methylseq \
  -r 4.2.0 \
  -profile "$PROFILE" \
  -c "$ALDAN3_CONFIG" \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  -work-dir "$WORKDIR" \
  "${REF_ARGS[@]}" \
  "${EXTRA_ARGS[@]}"

