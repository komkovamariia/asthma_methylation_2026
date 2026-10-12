#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."))
script <- file.path(root, "postanalysis", "postprocess_dss.R")

# Load function definitions without running the analysis entry point.
helpers <- new.env(parent = globalenv())
for (expression in parse(script)) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      is.call(expression[[3]]) && identical(expression[[3]][[1]], as.name("function"))) {
    eval(expression, helpers)
  }
}
stopifnot(identical(helpers$pick_sample("hla10_bismark.cov.gz", c("hla1", "hla10")), "hla10"))
stopifnot(is.na(helpers$pick_sample("unrelated.cov.gz", c("hla1", "hla10"))))

meta <- fread(file.path(root, "postanalysis", "sample_sheet.tsv"))
available <- meta[!sample_id %in% c("hla4", "hla5") & group == "patient"]
cd4 <- helpers$complete_pairs(available, cell_type = "CD4+")
cd8 <- helpers$complete_pairs(available, cell_type = "CD8+")
stopifnot(all(cd4$cell_type == "CD4+"), all(cd8$cell_type == "CD8+"))
stopifnot(!"patient_2" %in% cd4$patient_id, !"patient_1" %in% cd8$patient_id)
stopifnot(uniqueN(cd4$patient_id) == 14, uniqueN(cd8$patient_id) == 14)

workspace <- tempfile("methylseq_smoke_")
dir.create(workspace)
coverage_dir <- file.path(workspace, "coverage")
results_dir <- file.path(workspace, "results")
dir.create(coverage_dir)
amplicons <- c("IL-4", "IL-5", "IL-13", "CYSLTR1", "JAK3", "GATA3_Promoter",
               "GATA3_Exon1", "PRG2", "EPX", "RNASE2", "IL13RA1", "STAT6", "PTGDR",
               "ALOX15", "MMP9", "HDC", "FCER1G", "MMP2", "MRGPRX2")
set.seed(1031)
for (i in which(!meta$sample_id %in% c("hla4", "hla5"))) {
  positions <- rep(c(10L, 30L, 50L), times = length(amplicons))
  total <- rep(100L, length(positions))
  shift <- ifelse(meta$cell_type[i] == "CD8+", 0.10, 0) +
    ifelse(is.na(meta$timepoint[i]), 0, ifelse(meta$timepoint[i] == "after", 0.10, 0))
  probability <- pmin(0.9, pmax(0.1, 0.4 + shift + rnorm(length(positions), 0, 0.04)))
  methylated <- rbinom(length(positions), total, probability)
  rows <- data.frame(rep(amplicons, each = 3), positions, positions,
                     100 * methylated / total, methylated, total - methylated)
  handle <- gzfile(file.path(coverage_dir, paste0(meta$sample_id[i], "_bismark.cov.gz")), "w")
  write.table(rows, handle, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(handle)
}
log_file <- file.path(workspace, "postprocess.log")
status <- system2(file.path(R.home("bin"), "Rscript"), args = c(
  shQuote(script), "--methylseq_outdir", shQuote(coverage_dir),
  "--metadata_tsv", shQuote(file.path(root, "postanalysis", "sample_sheet.tsv")),
  "--results_dir", shQuote(results_dir)
), stdout = log_file, stderr = log_file)
if (status != 0) {
  cat(readLines(log_file), sep = "\n")
  stop("Synthetic postanalysis failed; fixtures retained in ", workspace)
}

mapping <- fread(file.path(results_dir, "coverage_file_map.tsv"))
stopifnot(nrow(mapping) == 64, !anyDuplicated(mapping$sample_id))
genes <- fread(file.path(results_dir, "gene_mean_methylation.tsv"))
stopifnot(uniqueN(genes$gene) == 19)
stopifnot(setequal(unique(genes[therapy == "GIPP", patient_id]), c("patient_13", "patient_14", "patient_15")))
stopifnot(setequal(unique(genes[therapy == "control", patient_id]), paste0("control_", 1:3)))
all_tables <- list.files(results_dir, pattern = "_DML_all.tsv$", full.names = TRUE)
stopifnot(length(all_tables) == 5)
for (path in all_tables) {
  all <- fread(path)
  p_column <- intersect(c("pval", "pvals", "p.value", "P.Value"), names(all))[[1]]
  stopifnot(isTRUE(all.equal(all$padj, p.adjust(all[[p_column]], "BH"))))
  significant <- fread(sub("_DML_all.tsv$", "_DML_sig.tsv", path))
  expected <- all[padj <= 0.05 & abs(mean_diff) >= 0.05]
  stopifnot(isTRUE(all.equal(significant, expected)))
}
stopifnot(length(list.files(results_dir, pattern = "^Q2_gene_.*pdf$")) == 76)
stopifnot(length(list.files(results_dir, pattern = "_paired_mean_methylation.pdf$")) == 3)
stopifnot(length(list.files(results_dir, pattern = "_volcano.pdf$")) == 5)
archive_results <- file.path(workspace, "archive_results")
archive_script <- file.path(root, "postanalysis", "postprocess_dss_all_genes_archive.R")
archive_status <- system2(file.path(R.home("bin"), "Rscript"), args = c(
  shQuote(archive_script), "--methylseq_outdir", shQuote(coverage_dir),
  "--metadata_tsv", shQuote(file.path(root, "postanalysis", "sample_sheet.tsv")),
  "--results_dir", shQuote(archive_results)
), stdout = log_file, stderr = log_file)
if (archive_status != 0) {
  cat(readLines(log_file), sep = "\n")
  stop("Historical paired-delta analysis failed; fixtures retained in ", workspace)
}
archive_tables <- list.files(archive_results, pattern = "_DML_all.tsv$", full.names = TRUE)
stopifnot(length(archive_tables) == 5)
for (path in archive_tables) {
  all <- fread(path)
  stopifnot(isTRUE(all.equal(all$padj, p.adjust(all$pval, "BH"))))
  significant <- fread(sub("_DML_all.tsv$", "_DML_sig.tsv", path))
  stopifnot(isTRUE(all.equal(significant, all[padj <= 0.05 & abs(diff) >= 0.05])))
}
stopifnot(length(list.files(archive_results, pattern = "_paired_mean_methylation.pdf$")) == 3)
unlink(workspace, recursive = TRUE)
cat("Current and historical R postanalysis checks passed: pairing, mapping, therapy, BH significance and plots.\n")
