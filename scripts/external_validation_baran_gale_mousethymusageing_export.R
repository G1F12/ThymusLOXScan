#!/usr/bin/env Rscript

# Export small Baran-Gale / MouseThymusAgeing metadata and LOX-family summaries.

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "scripts/external_validation_baran_gale_mousethymusageing_export.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
export_dir <- file.path(project_root, "data", "external", "BaranGale")
report_path <- file.path(project_root, "reports", "baran_gale_mousethymusageing_feasibility.md")
feasibility_tsv <- file.path(project_root, "results", "tables", "external_baran_gale_feasibility.tsv")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(feasibility_tsv), recursive = TRUE, showWarnings = FALSE)

write_feasibility <- function(reason, detail) {
  writeLines(
    c(
      "# Baran-Gale / MouseThymusAgeing feasibility",
      "",
      "Baran-Gale / MouseThymusAgeing remains a planned external TEC-focused analysis because the required expression/metadata could not be loaded in the current environment.",
      "",
      "## Reason",
      "",
      reason,
      "",
      "## Detail",
      "",
      detail
    ),
    report_path
  )
  write.table(
    data.frame(status = "metadata/data loading not possible", reason = reason, detail = detail),
    feasibility_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

required_packages <- c("MouseThymusAgeing", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "Matrix")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  write_feasibility(
    "Required R/Bioconductor packages are not installed.",
    paste(missing_packages, collapse = ", ")
  )
  quit(save = "no", status = 0)
}

safe_col <- function(df, names, default = NA_character_) {
  for (name in names) {
    if (name %in% colnames(df)) {
      return(as.character(df[[name]]))
    }
  }
  rep(default, nrow(df))
}

lox_genes <- c("Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4")
lox_ensembl <- c(
  Lox = "ENSMUSG00000030084",
  Loxl1 = "ENSMUSG00000032383",
  Loxl2 = "ENSMUSG00000034205",
  Loxl3 = "ENSMUSG00000026922",
  Loxl4 = "ENSMUSG00000029723"
)

tryCatch({
  sce <- MouseThymusAgeing::MouseSMARTseqData()
  if (is.list(sce)) {
    sce <- do.call(SingleCellExperiment::cbind, sce)
  }
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  rd <- as.data.frame(SummarizedExperiment::rowData(sce))

  gene_symbol <- safe_col(rd, c("Symbol", "symbol", "Gene", "gene", "GeneSymbol", "gene_symbol", "GeneName"))
  gene_id <- safe_col(rd, c("Geneid", "gene_id", "ensembl", "Ensembl", "ID", "id"), rownames(rd))
  row_lookup <- setNames(seq_len(nrow(rd)), gene_symbol)
  for (gene in names(lox_ensembl)) {
    if (!(gene %in% names(row_lookup))) {
      hit <- which(gene_id == lox_ensembl[[gene]])
      if (length(hit) > 0) {
        row_lookup[[gene]] <- hit[[1]]
      }
    }
  }
  missing_genes <- setdiff(lox_genes, names(row_lookup))
  if (length(missing_genes) > 0) {
    stop(paste("LOX-family genes absent from rowData:", paste(missing_genes, collapse = ", ")))
  }

  assay_names <- SummarizedExperiment::assayNames(sce)
  assay_name <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
  counts <- SummarizedExperiment::assay(sce, assay_name)
  library_size <- Matrix::colSums(counts)

  metadata <- data.frame(
    cell_id = safe_col(cd, c("CellID", "cell_id", "cell")),
    sample_id = safe_col(cd, c("sample", "Sample", "SortDay", "PlateID")),
    age = safe_col(cd, c("Age", "age", "age_group")),
    sort_day = safe_col(cd, c("SortDay", "sort_day")),
    sort_type = safe_col(cd, c("SortType", "sort_type")),
    subtype = safe_col(cd, c("SubType", "subtype", "cell_type")),
    size_factor = suppressWarnings(as.numeric(safe_col(cd, c("sizeFactor", "size_factor"), NA_character_))),
    raw_library_size = as.numeric(library_size),
    stringsAsFactors = FALSE
  )

  value_frames <- list()
  for (gene in lox_genes) {
    vals <- as.numeric(counts[row_lookup[[gene]], ])
    value_frames[[gene]] <- data.frame(
      cell_id = metadata$cell_id,
      sample_id = metadata$sample_id,
      age = metadata$age,
      sort_day = metadata$sort_day,
      sort_type = metadata$sort_type,
      subtype = metadata$subtype,
      raw_library_size = metadata$raw_library_size,
      gene = gene,
      raw_counts = vals,
      detected = vals > 0,
      stringsAsFactors = FALSE
    )
  }
  lox_values <- do.call(rbind, value_frames)

  metadata_summary <- data.frame(
    field = c("dataset", "technology", "n_cells", "age_labels", "sort_type_labels", "subtype_labels"),
    value = c(
      "MouseThymusAgeing MouseSMARTseqData",
      "SMART-seq2",
      as.character(nrow(metadata)),
      paste(sort(unique(metadata$age)), collapse = ";"),
      paste(sort(unique(metadata$sort_type)), collapse = ";"),
      paste(sort(unique(metadata$subtype)), collapse = ";")
    ),
    stringsAsFactors = FALSE
  )

  write.table(metadata_summary, file.path(export_dir, "baran_gale_metadata_summary_export.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(metadata, file.path(export_dir, "baran_gale_cell_metadata_export.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(lox_values, file.path(export_dir, "baran_gale_lox_cell_values_export.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved Baran-Gale MouseThymusAgeing export under ", export_dir)
}, error = function(e) {
  write_feasibility("MouseThymusAgeing export failed.", paste(class(e)[[1]], conditionMessage(e), sep = ": "))
  quit(save = "no", status = 0)
})
