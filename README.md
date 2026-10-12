# Targeted methylation analysis in asthma

Bisulfite sequencing of a 19-amplicon panel in CD4+ and CD8+ cells, before and after treatment. The project contains paired-end preprocessing with nf-core/methylseq 4.2.0, an ALDAN3 SLURM configuration and methylation analysis with DSS and bsseq.

## Analysis routes

| Path | Purpose | Status |
| --- | --- | --- |
| `first_try/` | Original directional PE nf-core run against hg38 | Historical baseline |
| `second_try/` | Non-directional PE nf-core run against the supplied 19-amplicon reference | Current preprocessing |
| `second_try/prep_merge_manifest.py` | Lane manifest for the attempted read-merging route | Historical preparation only; no merging implementation in this archive |
| `second_try_after/` | Independent hg38 alignment of PE lanes or available merged SE reads | Mapping diagnostics |
| `postanalysis/postprocess_dss.R` | Paired DSS models, Q3 and per-gene ASIT/GIPP plots | Current postanalysis |
| `postanalysis/postprocess_dss_all_genes_archive.R` | One-sample t-tests of paired methylation differences for Q1/Q2, DSS for Q3 and coverage-weighted global plots | Historical alternative method |

All distinct executable routes from the supplied archive are retained. Generated samplesheets and duplicate metadata are regenerated from the shared design. `down` and `git` were empty. The screenshots also show `all_contigs.txt`, which was absent from the archive; neither R script reads it.

## Samples and reference

The study has 15 patients with four samples each and three controls with two samples each. Patient samples follow this order: CD4+ before, CD8+ before, CD4+ after, CD8+ after. Samples `hla61` through `hla66` are the control CD4+/CD8+ pairs. `postanalysis/sample_sheet.tsv` records this design. GIPP includes patients 1, 7 and 10; the other patients receive ASIT. The preprocessing launcher excludes `hla4` and `hla5` by default, leaving 64 samples. `SKIP_SAMPLE_IDS` sets the exclusion list.

`second_try/create_tz_reference.py` contains the supplied bisulfite-converted sequences. Their names and sequences are preserved. Positions in downstream tables are local amplicon coordinates. A genomic coordinate map and the assignment's highlighted CpG subset are not included in this source archive.

## Run on ALDAN3

Run these commands from this directory. FASTQ files, reference indexes and analysis outputs are kept outside Git.

```bash
conda env create -f second_try/environment.yml
conda env create -f postanalysis/postprocess.yml
conda activate asthma_methylseq
cd second_try
mkdir -p logs
export RAW_FASTQ_DIR=/projects/mice_transplant_2025/raw/2025_12_kofiadi_seq
python create_tz_reference.py
```

For direct script execution, set `PROJECT_DIR` to the absolute path of `second_try` first:

```bash
export PROJECT_DIR="$PWD"
bash run_methylseq.sh
python make_tz_metadata.py --samplesheet_csv samplesheet_methylseq.csv \
  --out_tsv ../postanalysis/sample_sheet.tsv
```

For SLURM, submit from `second_try`. The submit directory becomes the default project directory. Create `logs` before submission because SLURM opens its log files before executing the script.

```bash
sbatch slurm_create_custom_bismark_index.sh
# Wait for a successful index build, then start methylseq.
sbatch slurm_methylseq_nfcore.sh
```

`slurm_methylseq_head.sh` can submit samplesheet preparation as a separate job. It scans FASTQs recursively, retains all lane rows and fails on missing or empty mates. Repeated sample IDs let nf-core concatenate lanes. No fastp read merging is used. The pipeline runs non-directional Bismark, skips deduplication for the amplicon library and requests cytosine reports. `RESUME=1` retains Nextflow resume behavior; set it to `0` for a fresh run. `REFERENCE_FASTA`, `BISMARK_INDEX_DIR`, `RESULTS_DIR` and `WORK_DIR` can be overridden. `TARGET_BED` enables the optional targeted-region and Qualimap steps; its coordinates must match the selected reference.

The head jobs request 8 CPUs and 64 GB. The uploaded process-specific settings, including the `infinite` partition and extended alignment/extraction time limits, remain in `aldan3_methylseq.config`. Review them for a different cluster.

After methylseq finishes:

```bash
cd ../postanalysis
conda activate asthma_methylseq_post
# PROJECT_DIR from preprocessing refers to a different directory.
unset PROJECT_DIR
bash run_postprocess.sh
```

`OUTDIR`, `METADATA_TSV` and `RESULTS_DIR` override the postanalysis paths. Defaults are the sibling `second_try/methylseq_out`, the supplied metadata and `postanalysis/methylseq_postprocess`.

## Statistical analysis and outputs

| Comparison | Design | Reported mean difference |
| --- | --- | --- |
| Q1: CD4+ versus CD8+ | Complete patient/timepoint pairs, DSS model `~ pair_id + cell_type` | CD4+ minus CD8+ |
| Q2: before versus after, within each cell type | Complete patient pairs, DSS model `~ patient_id + timepoint` | Before minus after |
| Q3: patients versus controls, within each cell type | Unpaired DSS two-group comparison | Patients minus controls |

Significant rows satisfy BH-adjusted `padj <= 0.05` and an absolute methylation difference of at least `0.05`, meaning five percentage points. The existing CLI name `--p_threshold` now refers to the BH-adjusted threshold. Runner overrides are `Q_THRESHOLD`, `DELTA`, `MIN_COV` and `TOP_N_HEATMAP`. BH correction is applied separately within each contrast. The minimum coverage of 5 retains a site when at least one sample passes it; it does not require coverage in every sample. Paired models use DSS's unsmoothed multifactor fit, while Q3 retains the uploaded `smoothing=TRUE` setting.

Gene-level plots use coverage-weighted methylation: the sum of methylated counts divided by the sum of total counts across covered sites in an amplicon. They retain sites with coverage at least 1. Each target has CD4+ and CD8+ plots for ASIT and GIPP when complete pairs exist. Lines join each patient's before/after values and purple lines show the controls. Plot colors use `after minus before`; changes below one percentage point are marked `no_change`. This descriptive threshold is separate from statistical significance.

Outputs include full and significant-site TSVs, adjusted-p-value volcano plots, paired mean-methylation plots, gene-level plots and coverage/PCA/heatmap QC. `gene_mean_methylation.tsv` exposes the weighted summaries; `coverage_file_map.tsv` records the exact input mapping. Missing timepoints or gene values are recorded in `gene_q2_exclusions.txt` when a gene plot has available input. Historical QC filenames containing `CpGs` are retained for compatibility.

The historical extractor setting `--CX` exports all cytosine contexts. Coverage files have no context column, so the postanalysis cannot verify CpG specificity from these files alone. CpG-specific inference requires context-validated coverage input. The cleaned configuration retains `--CX` inside its process selector and preserves nf-core's overlap and read-end options. Q3 pools available before/after patient samples, so repeated samples from a patient are not independent replicates; treat Q3 as exploratory. Q1 blocks patient/timepoint pairs, with two timepoints from some patients. Further inference about treatment or disease should consider these repeated measurements and the three-control sample size.

## Optional hg38 diagnostic alignment

`second_try_after` contains alignment-only diagnostics for poorly mapped libraries. It uses an existing hg38 Bismark index. Results are stored separately for each input lane, avoiding output collisions.

```bash
cd ../second_try_after
mkdir -p logs
export PROJECT_DIR="$PWD"
export INPUT_FASTQ_DIR=/projects/mice_transplant_2025/raw/2025_12_kofiadi_seq
export REF_DIR=/projects/mice_transplant_2025/reference
bash run_realign_hg38.sh
```

`ONLY_SAMPLE_IDS` or `ONLY_SAMPLES_FILE` restricts the inputs. The default excludes `Undetermined`, `hla4` and `hla5`. `PARTITION` sets the submitted partition and `ARRAY_CONCURRENCY` defaults to 10. Each task requests 16 CPUs and 128 GB; the default Bismark instance count is one quarter of the allocation because each instance launches additional workers. These BAMs are diagnostic outputs and do not replace the panel coverage files.

## Checks and changes from the uploaded archive

```bash
python -m unittest discover -s tests -v
Rscript tests/test_postprocess.R
```

The tests use temporary synthetic files. They check the panel sequence hash, all 66 metadata assignments, lane discovery, missing mates, sample exclusions and Nextflow options/exit status. R checks exercise paired sample selection and both postanalysis methods on synthetic coverage. No biological results are shipped with the source code.

Cleanup removes Jupyter checkpoints, generated absolute-path samplesheets and redundant launch notes. Historical analysis scripts remain in their original directories. Functional corrections in the current workflow restore GIPP patients 1/7/10, use BH-adjusted significance with the agreed 0.05 effect cutoff, filter CD4 and CD8 pairs correctly, name BSseq samples explicitly and restrict pairing to available coverage inputs. The hg38 baseline now resolves its own script directory and honors `RESUME`; its model and directional setting are retained. These corrections can change results relative to the uploaded scripts. Real-data reruns are required to establish biological conclusions.

The historical R alternative preserves its original statistics and global weighting. Its Q1 effect is CD8 minus CD4 and Q2 is after minus before, opposite to the current tables. It does not produce the therapy-stratified gene plots. Sites with constant paired differences can make the original t-test fail. Run it explicitly into a separate results directory:

```bash
Rscript postanalysis/postprocess_dss_all_genes_archive.R \
  --methylseq_outdir second_try/methylseq_out \
  --metadata_tsv postanalysis/sample_sheet.tsv \
  --results_dir postanalysis/methylseq_postprocess_archive
```

To reproduce the first hg38 route, activate `asthma_methylseq`, change into `first_try` and run `bash run_methylseq.sh` or `sbatch slurm_methylseq_head.sh`. Set `FASTA` and `BISMARK_INDEX` to the hg38 reference. This generator scans a flat FASTQ directory and includes all 66 samples. Its head job retains 6 CPUs, 32 GB and 20 days on `infinite`.

Methods: [nf-core/methylseq 4.2.0](https://nf-co.re/methylseq/4.2.0), [Bismark](https://github.com/FelixKrueger/Bismark), [DSS](https://bioconductor.org/packages/DSS/) and [bsseq](https://bioconductor.org/packages/bsseq/).
