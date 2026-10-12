#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
  library(bsseq)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

option_list <- list(
  make_option(c("--methylseq_outdir"), type="character", help="Path to methylseq outdir (e.g. methylseq_out)"),
  make_option(c("--metadata_tsv"), type="character", help="Path to sample metadata TSV"),
  make_option(c("--results_dir"), type="character", default="methylseq_postprocess", help="Directory for results"),
  make_option(c("--min_cov"), type="integer", default=5, help="Minimum coverage in at least one sample to keep a site (default 5)"),
  make_option(c("--p_threshold"), type="double", default=0.05, help="BH-adjusted p-value threshold (default 0.05)"),
  make_option(c("--delta"), type="double", default=0.05, help="Absolute methylation difference threshold (default 0.05)"),
  make_option(c("--top_n_heatmap"), type="integer", default=500, help="Top N coverage sites by variance for heatmap and PCA (default 500)")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$methylseq_outdir) || is.null(opt$metadata_tsv)) {
  stop("Both --methylseq_outdir and --metadata_tsv are required")
}
if (!is.finite(opt$p_threshold) || opt$p_threshold <= 0 || opt$p_threshold > 1 ||
    !is.finite(opt$delta) || opt$delta < 0 || opt$delta > 1 ||
    is.na(opt$min_cov) || opt$min_cov < 1 ||
    is.na(opt$top_n_heatmap) || opt$top_n_heatmap < 1) {
  stop("Invalid coverage, adjusted p-value, delta or heatmap threshold")
}

outdir <- normalizePath(opt$methylseq_outdir, mustWork=TRUE)
metadata_path <- normalizePath(opt$metadata_tsv, mustWork=TRUE)
results_dir <- opt$results_dir
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

message("methylseq outdir: ", outdir)
message("metadata:        ", metadata_path)
message("results_dir:     ", normalizePath(results_dir))

exclusion_log_file <- file.path(results_dir, "gene_q2_exclusions.txt")
if (file.exists(exclusion_log_file)) file.remove(exclusion_log_file)

meta <- fread(metadata_path)
stopifnot(all(c("sample_id", "patient_id", "cell_type", "timepoint", "group") %in% colnames(meta)))
if (anyDuplicated(meta$sample_id)) stop("Duplicate sample IDs in metadata")
if (anyNA(meta[, .(sample_id, patient_id, cell_type, group)])) stop("Missing required metadata values")
if (!all(meta$cell_type %in% c("CD4+", "CD8+"))) stop("Unknown cell_type in metadata")
if (!all(meta$group %in% c("patient", "control"))) stop("Unknown group in metadata")
if (!all(meta[group == "patient", timepoint] %in% c("before", "after"))) stop("Unknown patient timepoint")
if (!all(grepl("^patient_[0-9]+$", meta[group == "patient", patient_id]))) stop("Invalid patient_id")
if (anyDuplicated(meta[, .(patient_id, cell_type, timepoint)])) stop("Duplicate patient/cell/timepoint rows")

meta[group == "patient", patient_num := as.integer(sub("patient_", "", patient_id))]
meta[group == "patient", therapy := ifelse(patient_num %in% c(13, 14, 15), "GIPP", "ASIT")]
meta[group == "control", therapy := "control"]

cov_files <- list.files(outdir, pattern="\\.cov\\.gz$", recursive=TRUE, full.names=TRUE)
if (length(cov_files) == 0) stop("No .cov.gz files found under: ", outdir)

sample_ids <- unique(meta$sample_id)
pick_sample <- function(path, ids) {
  filename <- basename(path)
  hits <- ids[vapply(ids, function(id) {
    grepl(paste0("(^|[^[:alnum:]])", id, "([^[:alnum:]]|$)"), filename)
  }, logical(1))]
  if (length(hits) > 1) stop("Ambiguous sample ID in coverage file: ", path)
  if (length(hits) == 1) hits[[1]] else NA_character_
}

file_map <- data.table(
  file = cov_files,
  sample_id = vapply(cov_files, pick_sample, character(1), ids=sample_ids)
)

file_map <- file_map[sample_id %in% sample_ids]
file_map <- unique(file_map, by=c("file","sample_id"))

if (nrow(file_map) == 0) stop("Could not map any .cov.gz files to sample_ids from metadata.")

message("Found cov files mapped to metadata samples: ", nrow(file_map))

if (anyDuplicated(file_map$sample_id)) {
  stop("Multiple coverage files mapped to one sample. Use an outdir with one final file per sample.")
}
fwrite(file_map, file.path(results_dir, "coverage_file_map.tsv"), sep = "\t")
missing_samples <- setdiff(meta$sample_id, file_map$sample_id)
message("Samples without coverage files: ", paste(missing_samples, collapse = ", "))
sample_names <- file_map$sample_id
bs <- bsseq::read.bismark(files = file_map$file)
colnames(bs) <- sample_names
meta <- meta[match(sample_names, sample_id)]
stopifnot(identical(meta$sample_id, colnames(bs)))

# Gene summaries retain sites with coverage >= 1; DML uses a separate filter.
bs_all <- bs
cov_mat_all <- getCoverage(bs_all, type = "Cov")

# Retain sites with coverage >= min_cov in at least one sample.
keep <- rowSums(cov_mat_all >= opt$min_cov) >= 1
bs_dml <- bs_all[keep,]
if (nrow(bs_dml) == 0) stop("No coverage sites pass --min_cov")
cov_mat <- cov_mat_all[keep, , drop = FALSE]
message("Sites after coverage filter (DML): ", nrow(bs_dml))

meth_mat <- getMeth(bs_dml, type = "raw")  # proportion

gr_all <- rowRanges(bs_all)
chr_all <- as.character(GenomicRanges::seqnames(gr_all))

# Reference contig names differ from the gene labels used in reports.
gene_map <- list(
  "IL4"             = "IL-4",
  "IL5"             = "IL-5",
  "IL13"            = "IL-13",
  "CYSLTR1"         = "CYSLTR1",
  "JAK3"            = "JAK3",
  "GATA3(Promoter)" = "GATA3_Promoter",
  "GATA3(Exon1)"    = "GATA3_Exon1",
  "PRG2"            = "PRG2",
  "EPX"             = "EPX",
  "RNASE2"          = "RNASE2",
  "IL13RA1"         = "IL13RA1",
  "STAT6"           = "STAT6",
  "PTGDR"           = "PTGDR",
  "ALOX15"          = "ALOX15",
  "MMP9"            = "MMP9",
  "HDC"             = "HDC",
  "FCER1G"          = "FCER1G",
  "MMP2"            = "MMP2",
  "MRGPRX2"         = "MRGPRX2"
)

calc_gene_means <- function(gene_label, contig, min_cov_geneplot = 1L) {
  idx <- which(chr_all == contig)
  if (!length(idx)) return(NULL)

  cov_all <- getCoverage(bs_all, type = "Cov")[idx, , drop = FALSE]
  met_all <- getCoverage(bs_all, type = "M")[idx, , drop = FALSE]

  mask <- cov_all >= min_cov_geneplot
  den <- colSums(cov_all * mask, na.rm = TRUE)
  num <- colSums(met_all * mask, na.rm = TRUE)

  data.table(
    sample_id = colnames(cov_all),
    mean_meth = ifelse(den > 0, num / den, NA_real_),
    gene      = gene_label
  )
}

gene_means <- rbindlist(
  lapply(names(gene_map), function(glab) {
    calc_gene_means(glab, gene_map[[glab]], min_cov_geneplot = 1L)
  }),
  use.names = TRUE, fill = TRUE
)
if (!nrow(gene_means)) stop("No panel contigs found in coverage input; check the reference")
gene_means <- merge(gene_means, meta, by = "sample_id", all.x = TRUE)
fwrite(gene_means, file.path(results_dir, "gene_mean_methylation.tsv"), sep = "\t")

top_n <- min(opt$top_n_heatmap, nrow(meth_mat))
var_vec <- apply(meth_mat, 1, function(x) var(x, na.rm = TRUE))
top_idx <- order(var_vec, decreasing = TRUE)[seq_len(top_n)]
meth_top <- meth_mat[top_idx, , drop = FALSE]

row_ok <- apply(meth_top, 1, function(x) any(is.finite(x) & !is.na(x)))
col_ok <- apply(meth_top, 2, function(x) any(is.finite(x) & !is.na(x)))
meth_top <- meth_top[row_ok, col_ok, drop = FALSE]

annot <- as.data.frame(meta[match(colnames(meth_top), meta$sample_id), .(group, patient_id, cell_type, timepoint)])
rownames(annot) <- colnames(meth_top)

if (nrow(meth_top) >= 2 && ncol(meth_top) >= 2) {
  mt <- meth_top
  mt[!is.finite(mt)] <- NA
  row_ok2 <- apply(mt, 1, function(x) all(!is.na(x)))
  col_ok2 <- apply(mt, 2, function(x) all(!is.na(x)))
  mt <- mt[row_ok2, col_ok2, drop = FALSE]

  if (nrow(mt) >= 2 && ncol(mt) >= 2) {
    pca <- prcomp(t(mt), center = TRUE, scale. = FALSE)
    pca_df <- data.frame(
      sample_id = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2]
    )
    annot_pca <- annot[colnames(mt), , drop = FALSE]
    pca_df <- merge(pca_df, cbind(sample_id = rownames(annot_pca), annot_pca), by = "sample_id", all.x = TRUE)

    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type, shape = group)) +
      geom_point(size = 2.5, alpha = 0.9) +
      theme_bw() +
      labs(title = "PCA (top variable coverage sites)", subtitle = paste0("Top ", nrow(mt), " coverage sites by variance"))
    ggsave(file.path(results_dir, "QC_PCA_topCpGs.pdf"), p_pca, width = 7, height = 5)

    pdf(filepath <- file.path(results_dir, "QC_heatmap_topCpGs.pdf"), width = 9, height = 7)
    pheatmap(mt, show_rownames = FALSE, annotation_col = annot_pca, fontsize_col = 8)
    dev.off()
  } else {
    message("Skip PCA/heatmap after NA cleanup: ", nrow(mt), " CpGs, ", ncol(mt), " samples.")
  }
} else {
  message("Skip PCA/heatmap: not enough non-NA CpGs or samples (", nrow(meth_top), " CpGs, ", ncol(meth_top), " samples).")
}

cov_summary <- data.table(sample_id=colnames(cov_mat), mean_cov=colMeans(cov_mat), median_cov=apply(cov_mat,2,median))
cov_summary <- merge(cov_summary, meta, by="sample_id", all.x=TRUE)
p_cov <- ggplot(cov_summary, aes(x=cell_type, y=mean_cov, fill=group)) +
  geom_boxplot(outlier.size=0.6) +
  geom_jitter(width=0.15, size=1, alpha=0.7) +
  theme_bw() +
  labs(title="Mean coverage per sample", y="Mean coverage")
ggsave(file.path(results_dir, "QC_mean_coverage_boxplot.pdf"), p_cov, width=7, height=5)

get_pcol <- function(dt) {
  for (c in c("pval","pvals","p.value","P.Value")) {
    if (c %in% colnames(dt)) return(c)
  }
  stop("Cannot find p-value column in DSS result. Columns: ", paste(colnames(dt), collapse=", "))
}

mean_meth_diff <- function(bs_obj, group1_ids, group2_ids) {
  meth <- getMeth(bs_obj, type="raw")
  g1 <- intersect(group1_ids, colnames(meth))
  g2 <- intersect(group2_ids, colnames(meth))
  if (length(g1) == 0 || length(g2) == 0) stop("Empty groups for mean_meth_diff()")
  rowMeans(meth[, g1, drop=FALSE], na.rm=TRUE) - rowMeans(meth[, g2, drop=FALSE], na.rm=TRUE)
}

contrast_groups <- function(label) {
  if (startsWith(label, "Q1")) return(c("CD4+", "CD8+"))
  if (startsWith(label, "Q2")) return(c("before", "after"))
  c("patients", "controls")
}

contrast_text <- function(label) {
  paste(contrast_groups(label), collapse = " minus ")
}

plot_paired_mean <- function(bs_obj, design_df, group1_ids, group2_ids, label) {
  meth <- getMeth(bs_obj, type = "raw")
  means <- data.table(sample_id = colnames(meth), mean_meth = colMeans(meth, na.rm = TRUE))
  design <- as.data.table(design_df)
  block_column <- if ("pair_id" %in% names(design)) "pair_id" else "patient_id"
  design[, plot_pair := as.character(get(block_column))]
  means <- merge(means, design[, .(sample_id, plot_pair)], by = "sample_id")
  groups <- contrast_groups(label)
  means[, comparison := factor(ifelse(sample_id %in% group1_ids, groups[[1]], groups[[2]]),
                                levels = groups)]
  p <- ggplot(means, aes(x = comparison, y = mean_meth, group = plot_pair)) +
    geom_line(color = "grey65", alpha = 0.7) +
    geom_point(aes(color = comparison), size = 2) +
    theme_bw() +
    labs(title = paste0("Paired mean methylation: ", label), x = "", y = "Mean methylation")
  ggsave(file.path(results_dir, paste0(label, "_paired_mean_methylation.pdf")), p,
         width = 7, height = 5)
}

plot_volcano <- function(dt, label, diff_col, p_col) {
  volc <- data.frame(
    diff = dt[[diff_col]],
    neglog10p = -log10(pmax(dt$padj, 1e-300)),
    sig = (dt$padj <= opt$p_threshold) & (abs(dt[[diff_col]]) >= opt$delta)
  )
  p_volc <- ggplot(volc, aes(x=diff, y=neglog10p, color=sig)) +
    geom_point(size=0.7, alpha=0.7) +
    theme_bw() +
    labs(title=paste0("Volcano: ", label), x=paste0("Mean methylation difference (", contrast_text(label), ")"), y="-log10(BH-adjusted p-value)") +
    scale_color_manual(values=c("FALSE"="grey70","TRUE"="#d62728"))
  ggsave(file.path(results_dir, paste0(label, "_volcano.pdf")), p_volc, width=7, height=5)
}

plot_mean_meth_box <- function(bs_obj, meta_dt, group1_ids, group2_ids, label) {
  meth <- getMeth(bs_obj, type="raw")
  means <- data.frame(
    sample_id = colnames(meth),
    mean_meth = colMeans(meth, na.rm=TRUE)
  )
  means <- merge(means, meta_dt, by="sample_id", all.x=TRUE)
  group_labels <- contrast_groups(label)
  means$grp <- ifelse(means$sample_id %in% group1_ids, group_labels[[1]],
                      ifelse(means$sample_id %in% group2_ids, group_labels[[2]], NA))
  means$grp <- factor(means$grp, levels = group_labels)
  means <- means[!is.na(means$grp),]
  p_box <- ggplot(means, aes(x=grp, y=mean_meth)) +
    geom_boxplot(outlier.size=0.6) +
    geom_jitter(width=0.15, size=1, alpha=0.7) +
    theme_bw() +
    labs(title=paste0("Mean methylation per sample: ", label), x="", y="Mean methylation")
  ggsave(file.path(results_dir, paste0(label, "_mean_methylation_boxplot.pdf")), p_box, width=6, height=5)
}

make_q2_gene_plot <- function(gene_label, therapy_val, cell_val, delta_cut = 0.01) {
  dt <- gene_means[
    gene == gene_label &
      group == "patient" &
      therapy == therapy_val &
      cell_type == cell_val
  ]

  if (!nrow(dt)) {
    message("SKIP gene plot: ", gene_label, " ", therapy_val, " ", cell_val, " (no data)")
    return(NULL)
  }

  log_exclusion <- function(patients_vec, reason) {
    if (!length(patients_vec)) return(invisible(NULL))
    log_dt <- data.table(
      gene       = gene_label,
      therapy    = therapy_val,
      cell_type  = cell_val,
      patient_id = patients_vec,
      reason     = reason
    )
    data.table::fwrite(
      log_dt,
      file      = exclusion_log_file,
      sep       = "\t",
      append    = file.exists(exclusion_log_file),
      col.names = !file.exists(exclusion_log_file)
    )
  }

  all_patients <- unique(dt$patient_id)
  have_before <- unique(dt[timepoint == "before", patient_id])
  have_after  <- unique(dt[timepoint == "after",  patient_id])
  keep_patients <- intersect(have_before, have_after)

  log_exclusion(setdiff(all_patients, have_before), "no_before_timepoint")
  log_exclusion(setdiff(all_patients, have_after),  "no_after_timepoint")

  dt <- dt[patient_id %in% keep_patients]
  if (!nrow(dt)) {
    message("SKIP gene plot (no complete pairs): ", gene_label, " ", therapy_val, " ", cell_val)
    return(NULL)
  }

  wide <- dcast(dt, patient_id ~ timepoint, value.var = "mean_meth")
  wide_all <- copy(wide)
  wide <- wide[!is.na(before) & !is.na(after)]

  log_exclusion(setdiff(wide_all$patient_id, wide$patient_id), "NA_mean_meth_before_or_after_for_gene")

  if (!nrow(wide)) {
    message("SKIP gene plot (no before/after values): ", gene_label, " ", therapy_val, " ", cell_val)
    return(NULL)
  }

  wide[, delta := after - before]
  wide[, change := fifelse(abs(delta) < delta_cut, "no_change",
                           fifelse(delta > 0, "increase", "decrease"))]

  long <- melt(
    wide,
    id.vars = c("patient_id", "change"),
    measure.vars = c("before", "after"),
    variable.name = "timepoint",
    value.name = "mean_meth"
  )

  ctr <- gene_means[
    gene == gene_label &
      group == "control" &
      cell_type == cell_val
  ]
  if (nrow(ctr)) {
    ctr[, change := "control"]
  }

  lab_before <- wide[, .(patient_id, timepoint = "before", y = before)]
  lab_after  <- wide[, .(patient_id, timepoint = "after",  y = after)]

  long[, timepoint := factor(timepoint, levels = c("before", "after"))]

  p <- ggplot() +
    geom_line(
      data = long,
      aes(x = timepoint, y = mean_meth, group = patient_id, color = change),
      alpha = 0.7
    ) +
    geom_point(
      data = long,
      aes(x = timepoint, y = mean_meth, color = change),
      size = 1.6, alpha = 0.9
    ) +
    { if (nrow(ctr)) geom_hline(
        data = ctr,
        aes(yintercept = mean_meth, color = change),
        alpha = 0.8, linewidth = 0.7
      ) else NULL } +
    geom_text(
      data = lab_before,
      aes(x = timepoint, y = y, label = patient_id),
      hjust = 1.05, vjust = 0.5, size = 2.3, color = "grey30"
    ) +
    geom_text(
      data = lab_after,
      aes(x = timepoint, y = y, label = patient_id),
      hjust = -0.05, vjust = 0.5, size = 2.3, color = "grey30"
    ) +
    scale_color_manual(
      values = c(
        increase  = "red3",
        decrease  = "steelblue3",
        no_change = "grey70",
        control   = "purple"
      )
    ) +
    theme_bw() +
    labs(
      title    = paste0("gene: ", gene_label,
                        ", cell_type: ", cell_val,
                        ", therapy: ", therapy_val),
      subtitle = paste0("before vs after the therapy, |diff| < ", delta_cut, " means no_change"),
      x = "",
      y = "Coverage-weighted mean methylation"
    )

  fname <- sprintf("Q2_gene_%s_%s_%s.pdf",
                   gsub("[^A-Za-z0-9]+","_", gene_label),
                   gsub("[+ ]","", cell_val),
                   therapy_val)
  ggsave(file.path(results_dir, fname), p, width = 5.5, height = 4.5)
  invisible(p)
}

run_dml_two_group <- function(bs_obj, meta_dt, group1_ids, group2_ids, label) {
  g1 <- intersect(group1_ids, colnames(bs_obj))
  g2 <- intersect(group2_ids, colnames(bs_obj))
  if (length(g1) < 2 || length(g2) < 2) {
    message("SKIP ", label, ": need >=2 samples per group after intersection with BSseq (got ", length(g1), " vs ", length(g2), ")")
    return(NULL)
  }
  message("Running DSS DMLtest (unpaired): ", label, " (", length(g1), " vs ", length(g2), ")")
  dml <- DMLtest(bs_obj, group1 = g1, group2 = g2, smoothing = TRUE)
  dml_dt <- as.data.table(dml)
  p_col <- get_pcol(dml_dt)
  dml_dt[, padj := p.adjust(get(p_col), method = "BH")]
  dml_dt[, mean_diff := mean_meth_diff(bs_obj, g1, g2)]
  sig <- dml_dt[padj <= opt$p_threshold & abs(mean_diff) >= opt$delta]

  fwrite(dml_dt, file.path(results_dir, paste0(label, "_DML_all.tsv")), sep="\t")
  fwrite(sig,    file.path(results_dir, paste0(label, "_DML_sig.tsv")), sep="\t")
  plot_volcano(dml_dt, label, "mean_diff", p_col)
  plot_mean_meth_box(bs_obj, meta_dt, g1, g2, label)
  return(list(all=dml_dt, sig=sig))
}

complete_pairs <- function(meta_dt, cell_type=NULL) {
  m <- copy(meta_dt)
  if (!is.null(cell_type)) {
    selected_cell <- cell_type
    m <- m[get("cell_type") == selected_cell]
  }
  have_before <- unique(m[timepoint == "before", patient_id])
  have_after  <- unique(m[timepoint == "after",  patient_id])
  keep_patients <- intersect(have_before, have_after)
  m[patient_id %in% keep_patients]
}

run_dml_paired_multifactor <- function(bs_obj, meta_dt, label, group1_ids, group2_ids, design_df, formula_str, coef_name) {
  sample_ids <- unique(c(group1_ids, group2_ids))
  sample_ids <- intersect(sample_ids, colnames(bs_obj))
  if (length(sample_ids) < 6) {
    message("SKIP ", label, ": too few samples after intersection with BSseq")
    return(NULL)
  }

  bs_sub <- bs_obj[, sample_ids]
  design_df <- as.data.frame(design_df)
  rownames(design_df) <- design_df$sample_id
  design_df <- design_df[colnames(bs_sub), , drop = FALSE]
  design_df <- droplevels(design_df)
  design_matrix <- model.matrix(as.formula(formula_str), design_df)
  if (qr(design_matrix)$rank < ncol(design_matrix) ||
      nrow(design_matrix) <= ncol(design_matrix) + 1) {
    message("SKIP ", label, ": insufficient degrees of freedom for the paired DSS model")
    return(NULL)
  }

  message("Running DSS multifactor (paired): ", label, " n=", ncol(bs_sub), " formula=", formula_str, " coef=", coef_name)
  fit <- DMLfit.multiFactor(bs_sub, design=design_df, formula=as.formula(formula_str))
  res <- DMLtest.multiFactor(fit, coef=coef_name)
  dt <- as.data.table(res)
  p_col <- get_pcol(dt)
  dt[, padj := p.adjust(get(p_col), method="BH")]
  dt[, mean_diff := mean_meth_diff(bs_sub, group1_ids, group2_ids)]
  sig <- dt[padj <= opt$p_threshold & abs(mean_diff) >= opt$delta]

  fwrite(dt,  file.path(results_dir, paste0(label, "_DML_all.tsv")), sep="\t")
  fwrite(sig, file.path(results_dir, paste0(label, "_DML_sig.tsv")), sep="\t")
  plot_volcano(dt, label, "mean_diff", p_col)
  plot_paired_mean(bs_sub, design_df, group1_ids, group2_ids, label)
  return(list(all=dt, sig=sig))
}

patients <- meta[group == "patient"]
controls <- meta[group == "control"]

patients_cd4 <- complete_pairs(patients, cell_type="CD4+")
q2_cd4_g1 <- patients_cd4[timepoint == "before", sample_id]
q2_cd4_g2 <- patients_cd4[timepoint == "after",  sample_id]
design_cd4 <- patients_cd4[, .(sample_id, patient_id=factor(patient_id), timepoint=factor(timepoint, levels=c("before","after")))]
run_dml_paired_multifactor(bs_dml, meta, "Q2_paired_CD4_before_vs_after", q2_cd4_g1, q2_cd4_g2, design_cd4, "~ patient_id + timepoint", "timepointafter")

patients_cd8 <- complete_pairs(patients, cell_type="CD8+")
q2_cd8_g1 <- patients_cd8[timepoint == "before", sample_id]
q2_cd8_g2 <- patients_cd8[timepoint == "after",  sample_id]
design_cd8 <- patients_cd8[, .(sample_id, patient_id=factor(patient_id), timepoint=factor(timepoint, levels=c("before","after")))]
run_dml_paired_multifactor(bs_dml, meta, "Q2_paired_CD8_before_vs_after", q2_cd8_g1, q2_cd8_g2, design_cd8, "~ patient_id + timepoint", "timepointafter")

paired_cells <- patients[, .(complete = all(c("CD4+", "CD8+") %in% cell_type)),
                            by = .(patient_id, timepoint)][complete == TRUE]
patients_q1 <- merge(patients, paired_cells[, .(patient_id, timepoint)],
                     by = c("patient_id", "timepoint"))
patients_q1[, pair_id := interaction(patient_id, timepoint, drop = TRUE)]
q1_g1 <- patients_q1[cell_type == "CD4+", sample_id]
q1_g2 <- patients_q1[cell_type == "CD8+", sample_id]
design_q1 <- patients_q1[, .(
  sample_id,
  pair_id = factor(pair_id),
  cell_type = factor(cell_type, levels = c("CD4+", "CD8+"))
)]
run_dml_paired_multifactor(bs_dml, meta, "Q1_paired_patients_CD4_vs_CD8_adj_timepoint", q1_g1, q1_g2, design_q1, "~ pair_id + cell_type", "cell_typeCD8+")

# Controls are separate individuals; Q3 has no cross-group pairing.
q3_cd4_g1 <- patients[cell_type == "CD4+", sample_id]
q3_cd4_g2 <- controls[cell_type == "CD4+", sample_id]
run_dml_two_group(bs_dml, meta, q3_cd4_g1, q3_cd4_g2, "Q3_unpaired_CD4_patients_vs_controls")

q3_cd8_g1 <- patients[cell_type == "CD8+", sample_id]
q3_cd8_g2 <- controls[cell_type == "CD8+", sample_id]
run_dml_two_group(bs_dml, meta, q3_cd8_g1, q3_cd8_g2, "Q3_unpaired_CD8_patients_vs_controls")

for (g in names(gene_map)) {
  for (therapy_val in c("ASIT","GIPP")) {
    for (cell_val in c("CD4+","CD8+")) {
      make_q2_gene_plot(gene_label = g, therapy_val = therapy_val, cell_val = cell_val, delta_cut = 0.01)
    }
  }
}

message("Post-processing complete. Outputs in: ", normalizePath(results_dir))
