# Asthma methylation

19-amplicon bisulfite panel. nf-core/methylseq 4.2.0, Bismark and DSS on ALDAN3.

15 patients: CD4+ and CD8+, before and after treatment. Three controls: CD4+ and CD8+. GIPP: patients 13, 14 and 15. ASIT: patients 1 through 12. Design: `postanalysis/sample_sheet.tsv`. `hla4` and `hla5` are excluded by default: 64 of 66 samples.

## Workflow

nf-core checked read counts, lengths and base quality, then trimmed adapters and low-quality ends. Bismark aligned the reads to the custom 19-amplicon reference from the assignment. Coordinates are local to each amplicon.

Non-directional alignment searched four bisulfite strands: OT, original top; OB, original bottom; CTOT, complementary to original top; CTOB, complementary to original bottom. This allowed reads from any library orientation. Deduplication was disabled for the amplicon library.

Bismark reported mapped and unmapped reads and counted methylated and unmethylated observations at each covered position. MultiQC combined read QC, trimming, alignment and methylation-context reports. Run notes recorded low mapping and poor CHH/CHG metrics: biological interpretation remains limited.

The `.cov.gz` files were imported into bsseq as position-by-sample count matrices. Columns: target contig, 1-based start, end, methylation percentage, methylated count and unmethylated count. Coverage is the sum of both counts. The current `--CX` setting includes all cytosine contexts; these files have no context column. CpG-only summaries require context-validated input.

Per-target methylation was calculated as total methylated counts divided by total coverage, using sites with coverage ≥1. Metadata supplied patient ID, cell type, timepoint and patient/control status. Each target had four plots: CD4 ASIT, CD8 ASIT, CD4 GIPP and CD8 GIPP. Patients needed both before and after values with coverage for that target.

Plot change: after minus before. Absolute change <0.01 was `no_change`; larger changes were `increase` or `decrease`. Controls appeared as purple reference lines. Missing pairs or gene values were recorded in `gene_q2_exclusions.txt`; skipped plots were reported in the run log.

## Run

From the repository root:

```bash
conda env create -f second_try/environment.yml
conda env create -f postanalysis/postprocess.yml
conda activate asthma_methylseq
cd second_try
mkdir -p logs
export RAW_FASTQ_DIR=/projects/mice_transplant_2025/raw/2025_12_kofiadi_seq
python create_tz_reference.py
export PROJECT_DIR="$PWD"
bash run_methylseq.sh
python make_tz_metadata.py --samplesheet_csv samplesheet_methylseq.csv \
  --out_tsv ../postanalysis/sample_sheet.tsv
sbatch slurm_create_custom_bismark_index.sh
# Wait for the index job to finish successfully.
sbatch slurm_methylseq_nfcore.sh
```

After methylseq finishes:

```bash
cd ../postanalysis
conda activate asthma_methylseq_post
unset PROJECT_DIR
bash run_postprocess.sh
```

Preprocessing overrides: `SKIP_SAMPLE_IDS`, `RESUME`, `REFERENCE_FASTA`, `BISMARK_INDEX_DIR`, `RESULTS_DIR`, `WORK_DIR`. Optional `TARGET_BED` must match the reference. Postanalysis overrides: `OUTDIR`, `METADATA_TSV`, `RESULTS_DIR`, `MIN_COV`, `Q_THRESHOLD`, `DELTA`, `TOP_N_HEATMAP`.

## Statistics

| Contrast | Model | Table effect |
| --- | --- | --- |
| Q1: CD4 versus CD8 | DSS `~ pair_id + cell_type`, complete patient/timepoint pairs | CD4 minus CD8 |
| Q2: before versus after | DSS `~ patient_id + timepoint`, complete pairs within each cell type | Before minus after |
| Q3: patients versus controls | Unpaired DSS, separately for CD4 and CD8 | Patients minus controls |

DML filter: coverage ≥5 in at least one sample. Significance: BH-adjusted `padj <= 0.05` and `abs(mean_diff) >= 0.05`, per contrast. This is separate from the 0.01 plot-color cutoff. Q1 and Q2 use unsmoothed multifactor fits; Q3 uses smoothing and pools patient timepoints. Repeated patient samples and three controls limit inference.

Outputs: DML tables, volcano plots, paired plots, gene plots and QC. `gene_mean_methylation.tsv` contains target summaries; `coverage_file_map.tsv` records inputs.

## Other routes

| Path | Use |
| --- | --- |
| `first_try/` | Historical directional hg38 run |
| `second_try/prep_merge_manifest.py` | Historical lane manifest; creates a header-only samplesheet and performs no read merging |
| `second_try_after/` | hg38 alignment diagnostics, separate output per lane |
| `postanalysis/postprocess_dss_all_genes_archive.R` | Historical paired-delta t-tests for Q1/Q2, DSS for Q3 and weighted global plots |

For the first hg38 run, activate `asthma_methylseq`, set `FASTA` and `BISMARK_INDEX`, then run `bash run_methylseq.sh` from `first_try`. For diagnostics, set `PROJECT_DIR`, `INPUT_FASTQ_DIR` and `REF_DIR`, then run `bash run_realign_hg38.sh` from `second_try_after`.

Keep historical postanalysis output separate. Its Q1 effect is CD8 minus CD4; Q2 is after minus before. Constant paired differences can make its t-test fail.

```bash
Rscript postanalysis/postprocess_dss_all_genes_archive.R \
  --methylseq_outdir second_try/methylseq_out \
  --metadata_tsv postanalysis/sample_sheet.tsv \
  --results_dir postanalysis/methylseq_postprocess_archive
```

## Checks

```bash
python -m unittest discover -s tests -v
Rscript tests/test_postprocess.R
```

Tests use synthetic inputs. FASTQs, indexes, coverage files, logs and generated results stay outside Git.

[nf-core/methylseq](https://nf-co.re/methylseq/4.2.0) · [Bismark](https://github.com/FelixKrueger/Bismark) · [DSS](https://bioconductor.org/packages/DSS/) · [bsseq](https://bioconductor.org/packages/bsseq/).
