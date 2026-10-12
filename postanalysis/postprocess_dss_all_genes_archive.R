#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
  library(bsseq)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

option_list <- list(
  make_option(c("--methylseq_outdir"), type="character",
              help="Path to methylseq outdir, for example methylseq_out"),
  make_option(c("--metadata_tsv"), type="character",
              help="Path to sample_sheet.tsv with the study design"),
  make_option(c("--results_dir"), type="character", default="methylseq_postprocess",
              help="Directory for results"),
  make_option(c("--min_cov"), type="integer", default=5,
              help="Min coverage per CpG per sample to keep site, default 5"),
  make_option(c("--p_threshold"), type="double", default=NA,
              help="DEPRECATED: use --q_threshold, BH adjusted p value threshold"),
  make_option(c("--q_threshold"), type="double", default=0.05,
              help="BH adjusted p value threshold FDR, default 0.05"),
  make_option(c("--delta"), type="double", default=0.05,
              help="Min absolute methylation difference to call significant, default 0.05"),
  make_option(c("--top_n_heatmap"), type="integer", default=500,
              help="Top N CpGs by variance for heatmap and PCA, default 500"),
  make_option(c("--min_pairs"), type="integer", default=2,
              help="Min number of pairs required per CpG to run paired test, default 2")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (!is.na(opt$p_threshold)) opt$q_threshold <- opt$p_threshold

outdir <- normalizePath(opt$methylseq_outdir, mustWork=TRUE)
metadata_path <- normalizePath(opt$metadata_tsv, mustWork=TRUE)
results_dir <- opt$results_dir
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

message("methylseq outdir: ", outdir)
message("metadata:        ", metadata_path)
message("results_dir:     ", normalizePath(results_dir))

meta <- fread(metadata_path)
stopifnot(all(c("sample_id","patient_id","cell_type","timepoint","group") %in% colnames(meta)))

meta[, sample_id := trimws(sample_id)]
meta[, patient_id := trimws(patient_id)]
meta[, cell_type := trimws(cell_type)]
meta[, timepoint := trimws(timepoint)]
meta[, group := trimws(group)]

cov_files <- list.files(outdir, pattern="\\.cov\\.gz$", recursive=TRUE, full.names=TRUE)
if (length(cov_files) == 0) stop("No .cov.gz files found under: ", outdir)

sample_ids <- unique(meta$sample_id)

pick_sample <- function(path, ids) {
  bn <- basename(path)
  hits <- ids[sapply(ids, function(s) grepl(s, bn, fixed=TRUE))]
  if (length(hits) == 1) return(hits[[1]])
  if (length(hits) > 1) {
    hits2 <- hits[order(nchar(hits), decreasing=TRUE)]
    return(hits2[[1]])
  }
  return(sub("\\..*$","", bn))
}

file_map <- data.table(
  file = cov_files,
  sample_id = vapply(cov_files, pick_sample, character(1), ids=sample_ids)
)

file_map <- file_map[sample_id %in% sample_ids]
file_map <- unique(file_map, by=c("file","sample_id"))
if (nrow(file_map) == 0) stop("Could not map any .cov.gz files to sample_ids from metadata.")
message("Found cov files mapped to metadata samples: ", nrow(file_map))

sample_names <- file_map$sample_id
bs <- bsseq::read.bismark(files = file_map$file)

if (ncol(bs) != length(sample_names)) {
  stop("ncol bs differs from length sample_names, check file_map and file to sample mapping.")
}
colnames(bs) <- sample_names

meta <- meta[sample_id %in% colnames(bs)]
if (nrow(meta) == 0) stop("After intersecting metadata with BSseq sample names, meta is empty.")

meta <- meta[match(colnames(bs), sample_id)]
stopifnot(all(meta$sample_id == colnames(bs)))

message("BSseq sample names first 10: ", paste(head(colnames(bs), 10), collapse=", "))

bs_raw <- bs
cov_mat <- getCoverage(bs_raw, type="Cov")
keep <- rowSums(cov_mat >= opt$min_cov, na.rm=TRUE) >= 1
bs <- bs_raw[keep,]
message("CpGs after coverage filter: ", nrow(bs))

meth_mat <- getMeth(bs, type="raw") # proportions in [0,1], may contain NA

theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face="bold"),
      legend.title = element_text(face="bold"),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(filename, plot_obj, w=7, h=5) {
  ggsave(file.path(results_dir, filename), plot_obj, width=w, height=h)
}

plot_volcano <- function(dml_dt, title, xlab) {
  volc <- data.frame(
    diff = dml_dt$diff,
    neglog10q = -log10(pmax(dml_dt$padj, 1e-300))
  )
  volc$class <- "ns"
  volc$class[dml_dt$padj <= opt$q_threshold & dml_dt$diff >=  opt$delta] <- "hypermethylated"
  volc$class[dml_dt$padj <= opt$q_threshold & dml_dt$diff <= -opt$delta] <- "hypomethylated"

  ggplot(volc, aes(x = diff, y = neglog10q, color = class)) +
    geom_point(size = 1.0, alpha = 0.75) +
    theme_pub() +
    labs(title = title, x = xlab, y = "Minus log10 BH adjusted p") +
    scale_color_manual(values = c(
      "ns" = "grey75",
      "hypomethylated" = "#1f77b4",
      "hypermethylated" = "#d62728"
    )) +
    guides(color = guide_legend(title = "CpG class"))
}

top_n <- min(opt$top_n_heatmap, nrow(meth_mat))
if (top_n < 2) {
  message("SKIP QC PCA/heatmap: too few CpGs after filtering.")
} else {

  var_vec <- apply(meth_mat, 1, function(x) var(x, na.rm=TRUE))
  top_idx <- order(var_vec, decreasing=TRUE)[seq_len(top_n)]
  meth_top <- meth_mat[top_idx, , drop=FALSE]

  annot <- as.data.frame(meta[, .(group, patient_id, cell_type, timepoint)])
  rownames(annot) <- meta$sample_id

  meth_top2 <- meth_top[apply(meth_top, 1, function(r) sum(is.finite(r)) >= 2), , drop=FALSE]

  if (nrow(meth_top2) < 2) {
    message("SKIP QC PCA/heatmap: too few CpGs with finite values.")
  } else {
    row_means <- rowMeans(meth_top2, na.rm=TRUE)
    bad <- which(!is.finite(meth_top2), arr.ind=TRUE)
    if (nrow(bad) > 0) meth_top2[bad] <- row_means[bad[,1]]

    pca <- prcomp(t(meth_top2), center=TRUE, scale.=FALSE)
    pca_df <- data.frame(
      sample_id = rownames(pca$x),
      PC1 = pca$x[,1],
      PC2 = pca$x[,2]
    )
    pca_df <- merge(pca_df, cbind(sample_id=rownames(annot), annot), by="sample_id", all.x=TRUE)

    p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=cell_type, shape=group)) +
      geom_point(size=2.7, alpha=0.9) +
      theme_pub() +
      labs(title="PCA on top variable CpGs",
           subtitle=paste0("Top ", nrow(meth_top2), " CpGs by variance"),
           color="Cell type", shape="Group") +
      scale_color_brewer(palette="Set1")
    save_plot("QC_PCA_topCpGs.pdf", p_pca, w=7, h=5)

    hm_cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    pdf(file.path(results_dir, "QC_heatmap_topCpGs.pdf"), width=9, height=7)
    pheatmap(meth_top2,
             show_rownames=FALSE,
             annotation_col=annot[colnames(meth_top2), , drop=FALSE],
             fontsize_col=8,
             color=hm_cols,
             border_color=NA)
    dev.off()
  }
}

cov_summary <- data.table(sample_id=colnames(cov_mat),
                          mean_cov=colMeans(cov_mat, na.rm=TRUE),
                          median_cov=apply(cov_mat,2,median, na.rm=TRUE))
cov_summary <- merge(cov_summary, meta, by="sample_id", all.x=TRUE)

p_cov <- ggplot(cov_summary, aes(x=cell_type, y=mean_cov, fill=group)) +
  geom_boxplot(outlier.size=0.8, alpha=0.85) +
  geom_jitter(width=0.15, size=1.3, alpha=0.7) +
  theme_pub() +
  labs(title="Mean CpG coverage per sample",
       x="Cell type", y="Mean coverage", fill="Group") +
  scale_fill_brewer(palette="Set2")
save_plot("QC_mean_coverage_boxplot.pdf", p_cov, w=7, h=5)

# Historical route: paired deltas and one-sample t-tests for Q1 and Q2.
# Q3 uses DSS. Keep results separate from the current multifactor analysis.
has_multifactor <- FALSE

call_DMLfit_multiFactor <- function(bs_obj, design_mat, pdata_dt, formula_obj) {
  f <- DSS::DMLfit.multiFactor
  fmls <- names(formals(f))

  pdata <- as.data.frame(pdata_dt)
  if ("sample_id" %in% colnames(pdata)) rownames(pdata) <- pdata$sample_id

  attempts <- list()

  if ("formula" %in% fmls) {
    data_keys <- intersect(c("pData", "pdata", "phenoData", "data"), fmls)
    if (length(data_keys) == 0) data_keys <- NA_character_

    for (dk in data_keys) {
      args1 <- list(bs_obj, formula = formula_obj)
      if (!is.na(dk)) args1[[dk]] <- pdata
      if ("design" %in% fmls) args1[["design"]] <- design_mat
      if ("X" %in% fmls)      args1[["X"]]      <- design_mat
      attempts[[length(attempts) + 1]] <- args1

      args2 <- list(bs_obj, formula = formula_obj)
      if (!is.na(dk)) args2[[dk]] <- pdata
      attempts[[length(attempts) + 1]] <- args2
    }
  }

  if ("design" %in% fmls) attempts[[length(attempts) + 1]] <- list(bs_obj, design = design_mat)
  if ("X" %in% fmls)      attempts[[length(attempts) + 1]] <- list(bs_obj, X = design_mat)

  last_err <- NULL
  for (a in attempts) {
    res <- tryCatch(do.call(f, a), error = function(e) { last_err <<- e; NULL })
    if (!is.null(res)) return(res)
  }

  stop("Failed DSS DMLfit multiFactor calls. formals = ", paste(fmls, collapse = ", "),
       ". Last error: ", if (!is.null(last_err)) last_err$message else "unknown")
}

call_DMLtest_multiFactor <- function(fit_obj, coef_name) {
  g <- DSS::DMLtest.multiFactor
  gmls <- names(formals(g))
  if ("coef" %in% gmls) return(g(fit_obj, coef = coef_name))
  if ("coef.name" %in% gmls) return(g(fit_obj, coef.name = coef_name))
  stop("Unsupported DSS DMLtest multiFactor signature. formals = ", paste(gmls, collapse=", "))
}

rr <- rowRanges(bs)
chr_vec <- as.character(GenomeInfoDb::seqnames(rr))
pos_vec <- GenomicRanges::start(rr)

paired_delta_test <- function(delta_mat) {
  pvals <- apply(delta_mat, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < opt$min_pairs) return(NA_real_)
    suppressWarnings(t.test(x, mu=0)$p.value)
  })
  diff <- rowMeans(delta_mat, na.rm=TRUE)
  dt <- data.table(chr = chr_vec, pos = pos_vec, diff = diff, pval = pvals)
  dt <- dt[!is.na(pval)]
  dt[, padj := p.adjust(pval, method="BH")]
  dt
}

# Sum methylated counts over total coverage on sites with coverage >= 1.
global_mean_methylation <- function(bs_obj, sample_ids, min_cov_mean = 1) {
  cov <- getCoverage(bs_obj, type = "Cov")[, sample_ids, drop = FALSE]
  meth <- getMeth(bs_obj, type = "raw")[, sample_ids, drop = FALSE]
  ok <- (cov >= min_cov_mean) & is.finite(meth)
  num <- colSums(meth * cov * ok, na.rm = TRUE)
  den <- colSums(cov * ok, na.rm = TRUE)
  mean_meth <- num / den
  mean_meth[!is.finite(mean_meth) | den == 0] <- NA_real_
  data.table(sample_id = sample_ids, mean_meth = as.numeric(mean_meth))
}

run_Q1_paired_CD4_CD8 <- function(bs_obj, meta_dt, label) {
  dt <- as.data.table(meta_dt)[group == "patient"]
  dt[, cell := fifelse(cell_type == "CD4+", "CD4", "CD8")]
  dt[, tp := factor(timepoint, levels=c("before","after"))]
  dt[, pair_id := interaction(patient_id, tp, drop=TRUE)]

  wide <- dcast(dt, pair_id + patient_id + tp ~ cell, value.var="sample_id")
  wide <- wide[!is.na(CD4) & !is.na(CD8)]
  if (nrow(wide) < 2) {
    message("SKIP ", label, ": too few complete patient by timepoint pairs.")
    return(NULL)
  }

  samp <- c(wide$CD4, wide$CD8)
  bs_sub <- bs_obj[, samp]

  dml_dt <- NULL
  if (has_multifactor) {
    dml_dt <- tryCatch({
      dt2 <- dt[sample_id %in% samp]
      dt2[, cell := factor(cell, levels=c("CD4","CD8"))]
      dt2[, pair_id := factor(pair_id)]
      design <- model.matrix(~ cell + pair_id, data=dt2)
      fit <- call_DMLfit_multiFactor(bs_sub, design, dt2, ~ cell + pair_id)
      test <- call_DMLtest_multiFactor(fit, "cellCD8")  # CD8 - CD4
      dml_dt <- as.data.table(test)
      dml_dt[, padj := p.adjust(pval, method="BH")]
      if (!("chr" %in% names(dml_dt))) dml_dt[, chr := chr_vec]
      if (!("pos" %in% names(dml_dt))) dml_dt[, pos := pos_vec]
      dml_dt
    }, error = function(e) {
      message("DSS multiFactor failed for ", label, ": ", e$message,
              " | falling back to paired deltas")
      NULL
    })
  }
  if (is.null(dml_dt)) {
    meth <- getMeth(bs_sub, type="raw")
    delta_mat <- sapply(seq_len(nrow(wide)), function(i) {
      cd8 <- meth[, wide$CD8[i]]
      cd4 <- meth[, wide$CD4[i]]
      cd8 - cd4
    })
    if (is.null(dim(delta_mat))) delta_mat <- matrix(delta_mat, ncol=1)
    dml_dt <- paired_delta_test(delta_mat)
  }

  sig <- dml_dt[padj <= opt$q_threshold & abs(diff) >= opt$delta]

  fwrite(dml_dt, file.path(results_dir, paste0(label, "_DML_all.tsv")), sep="\t")
  fwrite(sig,    file.path(results_dir, paste0(label, "_DML_sig.tsv")), sep="\t")

  p_volc <- plot_volcano(dml_dt,
                         title = paste0("Volcano: ", label),
                         xlab  = "Delta methylation CD8 minus CD4")
  save_plot(paste0(label, "_volcano.pdf"), p_volc, w=7, h=5)

  means <- global_mean_methylation(bs_raw, samp, min_cov_mean = 1)

  long <- rbind(
    data.table(pair_id = wide$pair_id, patient_id = wide$patient_id, tp = wide$tp,
               cell = "CD4", sample_id = wide$CD4),
    data.table(pair_id = wide$pair_id, patient_id = wide$patient_id, tp = wide$tp,
               cell = "CD8", sample_id = wide$CD8)
  )
  long <- merge(long, means, by="sample_id", all.x=TRUE)

  message(label, " complete pairs: ", nrow(wide), " expected: ",
          uniqueN(dt$patient_id) * uniqueN(dt$tp))
  if (any(!is.finite(long$mean_meth))) {
    bad <- unique(long[!is.finite(mean_meth), sample_id])
    message(label, " samples with missing mean methylation: ",
            paste(head(bad, 10), collapse=", "),
            if (length(bad) > 10) " ..." else "")
  }

  long[, patient_id := factor(patient_id)]
  long[, tp := factor(tp, levels=c("before","after"))]
  long[, cell := factor(cell, levels=c("CD4","CD8"))]

  p_pair <- ggplot(long, aes(x=cell, y=mean_meth,
                             group=interaction(patient_id, tp),
                             color=patient_id, linetype=tp)) +
    geom_line(alpha=0.6, linewidth=0.85) +
    geom_point(size=2.4, alpha=0.95) +
    theme_pub() +
    labs(title="Paired comparison patients CD4 vs CD8",
         subtitle="Color shows patient, line type shows timepoint",
         x="Cell type", y="Mean methylation across CpGs",
         color="Patient", linetype="Timepoint") +
    scale_linetype_manual(values=c(before="solid", after="dashed")) +
    scale_color_hue()
  save_plot(paste0(label, "_paired_mean_methylation.pdf"), p_pair, w=9, h=6)

  invisible(list(all=dml_dt, sig=sig))
}

run_Q2_paired_before_after <- function(bs_obj, meta_dt, cell_label, label) {
  dt <- as.data.table(meta_dt)[group == "patient" & cell_type == cell_label]
  dt[, tp := factor(timepoint, levels=c("before","after"))]

  wide <- dcast(dt, patient_id ~ tp, value.var="sample_id")
  wide <- wide[!is.na(before) & !is.na(after)]
  if (nrow(wide) < 2) {
    message("SKIP ", label, ": too few complete before/after patient pairs.")
    return(NULL)
  }

  samp <- c(wide$before, wide$after)
  bs_sub <- bs_obj[, samp]

  dml_dt <- NULL
  if (has_multifactor) {
    dml_dt <- tryCatch({
      dt2 <- dt[sample_id %in% samp]
      dt2[, tp := factor(tp, levels=c("before","after"))]
      dt2[, patient_id := factor(patient_id)]
      design <- model.matrix(~ tp + patient_id, data=dt2)
      fit <- call_DMLfit_multiFactor(bs_sub, design, dt2, ~ tp + patient_id)
      test <- call_DMLtest_multiFactor(fit, "tpafter") # after - before
      dml_dt <- as.data.table(test)
      dml_dt[, padj := p.adjust(pval, method="BH")]
      if (!("chr" %in% names(dml_dt))) dml_dt[, chr := chr_vec]
      if (!("pos" %in% names(dml_dt))) dml_dt[, pos := pos_vec]
      dml_dt
    }, error = function(e) {
      message("DSS multiFactor failed for ", label, ": ", e$message,
              " | falling back to paired deltas")
      NULL
    })
  }
  if (is.null(dml_dt)) {
    meth <- getMeth(bs_sub, type="raw")
    delta_mat <- sapply(seq_len(nrow(wide)), function(i) {
      aft <- meth[, wide$after[i]]
      bef <- meth[, wide$before[i]]
      aft - bef
    })
    if (is.null(dim(delta_mat))) delta_mat <- matrix(delta_mat, ncol=1)
    dml_dt <- paired_delta_test(delta_mat)
  }

  sig <- dml_dt[padj <= opt$q_threshold & abs(diff) >= opt$delta]

  fwrite(dml_dt, file.path(results_dir, paste0(label, "_DML_all.tsv")), sep="\t")
  fwrite(sig,    file.path(results_dir, paste0(label, "_DML_sig.tsv")), sep="\t")

  p_volc <- plot_volcano(dml_dt,
                         title = paste0("Volcano: ", label),
                         xlab  = "Delta methylation after minus before")
  save_plot(paste0(label, "_volcano.pdf"), p_volc, w=7, h=5)

  samp <- c(wide$before, wide$after)
  means <- global_mean_methylation(bs_raw, samp, min_cov_mean = 1)
  setkey(means, sample_id)

  wide[, mean_before := means[before, mean_meth]]
  wide[, mean_after  := means[after,  mean_meth]]
  wide[, delta_mean  := mean_after - mean_before]
  wide[, direction := fifelse(delta_mean >= opt$delta, "increase",
                       fifelse(delta_mean <= -opt$delta, "decrease", "no_change"))]

  long <- rbind(
    data.table(patient_id = wide$patient_id, tp = "before", mean_meth = wide$mean_before, direction = wide$direction),
    data.table(patient_id = wide$patient_id, tp = "after",  mean_meth = wide$mean_after,  direction = wide$direction)
  )
  long[, tp := factor(tp, levels=c("before","after"))]
  long[, direction := factor(direction, levels=c("decrease","no_change","increase"))]

  message(label, " complete patient pairs: ", nrow(wide), " expected: ", uniqueN(dt$patient_id))
  if (any(!is.finite(long$mean_meth))) {
    badp <- unique(long[!is.finite(mean_meth), patient_id])
    message(label, " patients with missing mean methylation: ",
            paste(head(badp, 10), collapse=", "),
            if (length(badp) > 10) " ..." else "")
  }

  lab_before <- long[tp == "before"]
  lab_after  <- long[tp == "after"]

  p_pair <- ggplot(long, aes(x=tp, y=mean_meth, group=patient_id, color=direction)) +
    geom_line(alpha=0.55, linewidth=0.85) +
    geom_point(size=2.6, alpha=0.95) +
    ggrepel::geom_text_repel(data=lab_before, aes(label=patient_id),
                             nudge_x=-0.25, direction="y",
                             size=3, max.overlaps=Inf, min.segment.length=0) +
    ggrepel::geom_text_repel(data=lab_after, aes(label=patient_id),
                             nudge_x=0.25, direction="y",
                             size=3, max.overlaps=Inf, min.segment.length=0) +
    theme_pub() +
    theme(plot.margin = margin(12, 70, 12, 70)) +
    scale_y_continuous(expand=expansion(mult=c(0.06, 0.25))) +
    labs(title=paste0("Paired comparison patients ", cell_label, " before vs after"),
         subtitle=paste0("Direction is based on delta ", opt$delta),
         x="Timepoint", y="Mean methylation across CpGs", color="Change") +
    scale_color_manual(values=c(decrease="#d62728", no_change="#7f7f7f", increase="#2ca02c"))
  save_plot(paste0(label, "_paired_mean_methylation.pdf"), p_pair, w=10, h=7)

  invisible(list(all=dml_dt, sig=sig))
}

run_Q3_unpaired_patients_controls <- function(bs_obj, meta_dt, cell_label, label) {
  dt_pat <- as.data.table(meta_dt)[group == "patient" & cell_type == cell_label]
  dt_ctl <- as.data.table(meta_dt)[group == "control" & cell_type == cell_label]

  ids_pat <- dt_pat$sample_id
  ids_ctl <- dt_ctl$sample_id

  if (length(ids_pat) < 2 || length(ids_ctl) < 2) {
    message("SKIP ", label, ": need at least 2 samples per group, patients=", length(ids_pat),
            ", controls=", length(ids_ctl))
    return(NULL)
  }

  message("Running unpaired DMLtest: ", label, " patients n=", length(ids_pat),
          " vs controls n=", length(ids_ctl))

  dml <- DSS::DMLtest(bs_obj, ids_pat, ids_ctl, smoothing = TRUE)
  dml_dt <- as.data.table(dml)
  dml_dt[, padj := p.adjust(pval, method="BH")]

  sig <- dml_dt[padj <= opt$q_threshold & abs(diff) >= opt$delta]

  fwrite(dml_dt, file.path(results_dir, paste0(label, "_DML_all.tsv")), sep="\t")
  fwrite(sig,    file.path(results_dir, paste0(label, "_DML_sig.tsv")), sep="\t")

  p_volc <- plot_volcano(dml_dt,
                         title = paste0("Volcano: ", label),
                         xlab  = "Delta methylation patients minus controls")
  save_plot(paste0(label, "_volcano.pdf"), p_volc, w=7, h=5)

  samp <- c(ids_ctl, ids_pat)
  means <- global_mean_methylation(bs_raw, samp, min_cov_mean = 1)
  means <- merge(means, meta_dt[, .(sample_id, group, cell_type, timepoint)], by="sample_id", all.x=TRUE)
  means <- means[cell_type == cell_label & group %in% c("patient","control")]

  if (any(is.na(means$timepoint) & means$group == "control")) {
    means[group == "control" & is.na(timepoint), timepoint := "before"]
    message(label, " controls have missing timepoint, set to before for plotting")
  }

  means[, cohort := ifelse(group == "patient", "patients", "controls")]
  means[, cohort := factor(cohort, levels=c("controls","patients"))]
  means[, timepoint := factor(timepoint, levels=c("before","after"))]

  message(label, " sample counts by cohort and timepoint:")
  print(table(means$cohort, means$timepoint, useNA="ifany"))

  if (any(!is.finite(means$mean_meth))) {
    bad <- means[!is.finite(mean_meth), .(sample_id, cohort, timepoint)]
    message(label, " samples with missing mean methylation:")
    print(head(bad, 10))
  }

  p_box <- ggplot(means, aes(x=cohort, y=mean_meth, color=timepoint)) +
    geom_point(position=position_jitterdodge(jitter.width=0.12, dodge.width=0.5),
               size=2.4, alpha=0.9) +
    theme_pub() +
    labs(title=paste0("Patients vs controls ", cell_label, " mean methylation"),
         x="Cohort", y="Mean methylation across CpGs", color="Timepoint") +
    scale_color_manual(values=c(before="#1f77b4", after="#d62728"))
  save_plot(paste0(label, "_mean_methylation_boxplot.pdf"), p_box, w=7.5, h=5.5)

  invisible(list(all=dml_dt, sig=sig))
}

message("Using DSS multiFactor: ", ifelse(has_multifactor, "YES", "NO fallback uses paired deltas and one sample t test"))

run_Q1_paired_CD4_CD8(bs, meta, "Q1_patients_CD4_vs_CD8_paired")
run_Q2_paired_before_after(bs, meta, "CD4+", "Q2_patients_CD4_before_vs_after_paired")
run_Q2_paired_before_after(bs, meta, "CD8+", "Q2_patients_CD8_before_vs_after_paired")
run_Q3_unpaired_patients_controls(bs, meta, "CD4+", "Q3_CD4_patients_vs_controls")
run_Q3_unpaired_patients_controls(bs, meta, "CD8+", "Q3_CD8_patients_vs_controls")

message("Post-processing complete. Outputs in: ", normalizePath(results_dir))
