#!/usr/bin/env Rscript

# Export only LOX-family data and metadata from the official MouseThymusAgeing API.
# This script intentionally avoids raw sequencing/alignment files and full matrices.

project_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE)
out_dir <- file.path(project_root, "data", "external", "metadata", "emtab8560")
report_path <- file.path(project_root, "reports", "external_emtab8560_r_export_report.md")
feasibility_path <- file.path(project_root, "reports", "external_emtab8560_r_export_feasibility.md")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

required <- c(
  "MouseThymusAgeing",
  "SingleCellExperiment",
  "SummarizedExperiment"
)
optional <- c("scuttle", "scater", "dplyr", "readr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

write_feasibility <- function(lines) {
  writeLines(c(
    "# E-MTAB-8560 R Export Feasibility",
    "",
    lines,
    "",
    "Required installation example:",
    "",
    "```r",
    "install.packages('BiocManager')",
    "BiocManager::install(c('MouseThymusAgeing', 'SingleCellExperiment', 'SummarizedExperiment', 'scuttle'))",
    "```",
    "",
    "No E-MTAB-8560 per-mouse inference was performed."
  ), feasibility_path)
}

if (length(missing) > 0) {
  write_feasibility(c(
    paste0("Missing required R/Bioconductor packages: ", paste(missing, collapse = ", "), "."),
    "The official MouseThymusAgeing API could not be loaded.",
    "No Python/RDS fallback was used."
  ))
  stop(paste("Missing required packages:", paste(missing, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(MouseThymusAgeing)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

lox_genes <- c("Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4")

smart_days <- NULL
if (exists("SMARTseqMetadata")) {
  meta <- get("SMARTseqMetadata")
  if ("sample" %in% colnames(meta)) {
    smart_days <- unique(as.character(meta$sample))
  }
}
if (is.null(smart_days) || length(smart_days) == 0) {
  smart_days <- paste0("day", 1:5)
}

sce_list <- list()
for (day in smart_days) {
  sce <- MouseSMARTseqData(samples = day)
  sce$ExportSample <- day
  sce_list[[day]] <- sce
}
sce <- do.call(SingleCellExperiment::cbind, sce_list)

assay_inventory <- data.frame(
  assay_name = SummarizedExperiment::assayNames(sce),
  n_features = nrow(sce),
  n_cells = ncol(sce),
  stringsAsFactors = FALSE
)
write.table(assay_inventory, file.path(out_dir, "emtab8560_assay_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

coldata <- as.data.frame(SummarizedExperiment::colData(sce))
coldata$cell_barcode <- colnames(sce)
write.table(coldata, file.path(out_dir, "emtab8560_coldata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rowdata <- as.data.frame(SummarizedExperiment::rowData(sce))
rowdata$feature_id <- rownames(sce)
write.table(rowdata, file.path(out_dir, "emtab8560_rowdata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

gene_cols <- grep("symbol|gene", colnames(rowdata), ignore.case = TRUE, value = TRUE)
feature_symbol <- rep(NA_character_, nrow(rowdata))
for (col in gene_cols) {
  vals <- as.character(rowdata[[col]])
  hit <- vals %in% lox_genes
  feature_symbol[hit] <- vals[hit]
}
if (any(is.na(feature_symbol))) {
  feature_symbol[rownames(sce) %in% lox_genes] <- rownames(sce)[rownames(sce) %in% lox_genes]
}
lox_idx <- which(feature_symbol %in% lox_genes)

export_assay <- function(assay_name, out_name) {
  if (!(assay_name %in% SummarizedExperiment::assayNames(sce)) || length(lox_idx) == 0) {
    return(FALSE)
  }
  mat <- SummarizedExperiment::assay(sce, assay_name)[lox_idx, , drop = FALSE]
  df <- as.data.frame(as.matrix(mat))
  df <- cbind(
    feature_id = rownames(mat),
    gene = feature_symbol[lox_idx],
    df
  )
  write.table(df, file.path(out_dir, out_name), sep = "\t", quote = FALSE, row.names = FALSE)
  TRUE
}

counts_exported <- export_assay("counts", "emtab8560_lox_counts.tsv")
logcounts_exported <- export_assay("logcounts", "emtab8560_lox_logcounts.tsv")

individual_fields <- grep("individual|mouse|animal", colnames(coldata), ignore.case = TRUE, value = TRUE)

writeLines(c(
  "# E-MTAB-8560 R Export Report",
  "",
  "MouseThymusAgeing was loaded through the official R API.",
  "",
  "## Assays",
  "",
  paste0("- ", assay_inventory$assay_name),
  "",
  "## Exports",
  "",
  "- `data/external/metadata/emtab8560/emtab8560_coldata.tsv`",
  "- `data/external/metadata/emtab8560/emtab8560_rowdata.tsv`",
  "- `data/external/metadata/emtab8560/emtab8560_assay_inventory.tsv`",
  if (counts_exported) "- `data/external/metadata/emtab8560/emtab8560_lox_counts.tsv`" else "- LOX raw counts not exported.",
  if (logcounts_exported) "- `data/external/metadata/emtab8560/emtab8560_lox_logcounts.tsv`" else "- LOX logcounts not exported.",
  "",
  "## Biological ID Fields",
  "",
  if (length(individual_fields) > 0) {
    paste0("Candidate biological mouse/individual fields in colData: ", paste(individual_fields, collapse = ", "), ".")
  } else {
    "No biological mouse/individual field was found in exported colData."
  },
  "",
  "No FASTQ, BAM, SRA, or raw alignment files were downloaded."
), report_path)

message("E-MTAB-8560 MouseThymusAgeing export completed.")
