#!/usr/bin/env bash
#SBATCH --job-name=hg38_realign
#SBATCH --output=logs/%x_%j_%a.out
#SBATCH --error=logs/%x_%j_%a.err
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

: "${SAMPLESHEET_TSV:?Set SAMPLESHEET_TSV}"
: "${ALIGN_DIR:?Set ALIGN_DIR}"
: "${TMP_DIR:?Set TMP_DIR}"
: "${REF_DIR:?Set REF_DIR}"
: "${SLURM_ARRAY_TASK_ID:?This script requires a SLURM array task}"
: "${SLURM_JOB_ID:?This script requires a SLURM job}"
command -v bismark >/dev/null
command -v samtools >/dev/null

TOTAL=$(($(wc -l < "${SAMPLESHEET_TSV}") - 1))
IDX="${SLURM_ARRAY_TASK_ID}"
if [[ "${IDX}" -lt 1 || "${IDX}" -gt "${TOTAL}" ]]; then
  echo "Invalid array index ${IDX} for ${TOTAL} input rows" >&2
  exit 2
fi
line=$(awk -v n=$((IDX + 1)) 'NR == n {print; exit}' "${SAMPLESHEET_TSV}")
line="${line%$'\r'}"
sample=$(printf '%s\n' "${line}" | cut -f1)
fq1=$(printf '%s\n' "${line}" | cut -f2)
fq2=$(printf '%s\n' "${line}" | cut -f3)
for fastq in "${fq1}" "${fq2}"; do
  [[ -z "${fastq}" ]] && continue
  [[ -r "${fastq}" && -s "${fastq}" ]] || { echo "Missing or empty FASTQ: ${fastq}" >&2; exit 2; }
done
[[ -d "${REF_DIR}/Bisulfite_Genome" ]] || { echo "Missing Bismark index: ${REF_DIR}" >&2; exit 2; }

# Lane tasks for the same sample need separate outputs and temporary directories.
TASK_ALIGN_DIR="${ALIGN_DIR}/${sample}/task_${IDX}"
temp="${TMP_DIR}/tmp_${sample}_${SLURM_JOB_ID}_${IDX}"
mkdir -p "${TASK_ALIGN_DIR}" "${temp}"
trap 'rm -rf -- "${temp}"' EXIT
THREADS="${SLURM_CPUS_PER_TASK:-16}"
# Each Bismark instance also runs Bowtie2 and compression workers.
BISMARK_MULTICORE="${BISMARK_MULTICORE:-$((THREADS / 4))}"
[[ "${BISMARK_MULTICORE}" -ge 1 ]] || BISMARK_MULTICORE=1
export OMP_NUM_THREADS="${THREADS}"
align_args=(--genome "${REF_DIR}" --output_dir "${TASK_ALIGN_DIR}" --temp_dir "${temp}"
  --basename "${sample}" --multicore "${BISMARK_MULTICORE}"
  --bowtie2 --score_min "L,0,-0.6" --minins 0 --maxins 2000 --non_directional)

if [[ -n "${fq2}" ]]; then
  bismark "${align_args[@]}" -1 "${fq1}" -2 "${fq2}"
  raw_bam="${TASK_ALIGN_DIR}/${sample}_pe.bam"
  report="${TASK_ALIGN_DIR}/${sample}_PE_report.txt"
else
  bismark "${align_args[@]}" "${fq1}"
  raw_bam="${TASK_ALIGN_DIR}/${sample}.bam"
  report="${TASK_ALIGN_DIR}/${sample}_SE_report.txt"
fi
[[ -s "${raw_bam}" ]] || { echo "Expected BAM not found: ${raw_bam}" >&2; exit 3; }
sorted_bam="${TASK_ALIGN_DIR}/${sample}.hg38.sorted.bam"
samtools sort -@ "${THREADS}" -o "${sorted_bam}" "${raw_bam}"
samtools index -@ "${THREADS}" "${sorted_bam}"
samtools flagstat "${sorted_bam}" > "${TASK_ALIGN_DIR}/${sample}.flagstat.txt"
samtools idxstats "${sorted_bam}" > "${TASK_ALIGN_DIR}/${sample}.idxstats.txt"
if [[ -f "${report}" ]]; then
  grep -E 'Mapping efficiency|Total number|Number of alignments|CT/GA conversion|C methylated in CpG context' \
    "${report}" > "${TASK_ALIGN_DIR}/${sample}.bismark_key_metrics.txt" || true
fi
printf 'Finished %s, input row %s\n' "${sample}" "${IDX}"
