#!/usr/bin/env Rscript

# Validate the curated Sha/REMoDLB dataset files without modifying source data.

dataset_dir <- "sha"
report_dir <- "reports"
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  expr_total = file.path(dataset_dir, "sha_total_expression.rds"),
  expr_rchop = file.path(dataset_dir, "sha_rchop_expression.rds"),
  expr_rbchop = file.path(dataset_dir, "sha_rbchop_expression.rds"),
  clin_total = file.path(dataset_dir, "clin_sha_tot.rds"),
  clin_rchop = file.path(dataset_dir, "clin_sha_rchop.rds"),
  clin_rbchop = file.path(dataset_dir, "clin_sha_rbchop.rds")
)

stop_if_missing <- function(files) {
  missing <- files[!file.exists(unlist(files))]
  if (length(missing) > 0) {
    stop("Missing expected file(s): ", paste(unlist(missing), collapse = ", "))
  }
}

object_summary <- function(path) {
  object <- readRDS(path)
  data.frame(
    file = path,
    class = paste(class(object), collapse = ", "),
    rows = if (!is.null(dim(object))) nrow(object) else NA_integer_,
    columns = if (!is.null(dim(object))) ncol(object) else NA_integer_,
    duplicated_rows = if (!is.null(rownames(object))) anyDuplicated(rownames(object)) else NA_integer_,
    duplicated_columns = if (!is.null(colnames(object))) anyDuplicated(colnames(object)) else NA_integer_,
    any_na = anyNA(object),
    stringsAsFactors = FALSE
  )
}

md_table <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data[] <- lapply(data, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  header <- paste(names(data), collapse = " | ")
  divider <- paste(rep("---", ncol(data)), collapse = " | ")
  rows <- apply(data, 1, paste, collapse = " | ")
  paste(c(header, divider, rows), collapse = "\n")
}

fmt_bool <- function(value) {
  if (isTRUE(value)) "TRUE" else "FALSE"
}

stop_if_missing(paths)

expr_total <- readRDS(paths$expr_total)
expr_rchop <- readRDS(paths$expr_rchop)
expr_rbchop <- readRDS(paths$expr_rbchop)
clin_total <- readRDS(paths$clin_total)
clin_rchop <- readRDS(paths$clin_rchop)
clin_rbchop <- readRDS(paths$clin_rbchop)

summaries <- do.call(rbind, lapply(paths, object_summary))

sample_checks <- data.frame(
  check = c(
    "total expression columns match total clinical Sample_ID",
    "R-CHOP expression columns match R-CHOP clinical Sample_ID",
    "RB-CHOP expression columns match RB-CHOP clinical Sample_ID",
    "total clinical Sample_ID has no duplicates",
    "total expression column names have no duplicates",
    "total expression row names have no duplicates",
    "total expression and clinical are already in identical order"
  ),
  result = c(
    all(colnames(expr_total) %in% clin_total$Sample_ID) && all(clin_total$Sample_ID %in% colnames(expr_total)),
    all(colnames(expr_rchop) %in% clin_rchop$Sample_ID) && all(clin_rchop$Sample_ID %in% colnames(expr_rchop)),
    all(colnames(expr_rbchop) %in% clin_rbchop$Sample_ID) && all(clin_rbchop$Sample_ID %in% colnames(expr_rbchop)),
    anyDuplicated(clin_total$Sample_ID) == 0,
    anyDuplicated(colnames(expr_total)) == 0,
    anyDuplicated(rownames(expr_total)) == 0,
    identical(colnames(expr_total), clin_total$Sample_ID)
  ),
  stringsAsFactors = FALSE
)

treatment_counts <- as.data.frame(table(clin_total$Treatment, useNA = "ifany"), stringsAsFactors = FALSE)
names(treatment_counts) <- c("Treatment", "n")

class_counts <- as.data.frame(table(clin_total$Class, useNA = "ifany"), stringsAsFactors = FALSE)
names(class_counts) <- c("Class", "n")

survival_missing <- data.frame(
  field = c("Mesi_OS", "Evento_OS", "Mesi_PFS", "Evento_PFS"),
  missing_n = colSums(is.na(clin_total[, c("Mesi_OS", "Evento_OS", "Mesi_PFS", "Evento_PFS")])),
  stringsAsFactors = FALSE
)

clinical_missing <- data.frame(
  field = names(clin_total),
  missing_n = colSums(is.na(clin_total)),
  missing_pct = round(100 * colSums(is.na(clin_total)) / nrow(clin_total), 2),
  stringsAsFactors = FALSE
)
clinical_missing <- clinical_missing[order(-clinical_missing$missing_n, clinical_missing$field), ]

expr_range <- range(as.matrix(expr_total), na.rm = TRUE)

supplementary_files <- list.files(
  file.path(dataset_dir, "sha_supplementary_tables"),
  all.files = FALSE,
  full.names = FALSE
)
supplementary_files <- supplementary_files[order(supplementary_files)]
lock_files <- supplementary_files[grepl("^~\\$", supplementary_files)]

report <- c(
  "# Sha/REMoDLB Curated Dataset Validation",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Provenance Anchor",
  "",
  "- Local files match the Sha et al. 2019 JCO / REMoDLB public dataset structure.",
  "- GEO accession: GSE117556.",
  "- PubMed ID: 30523719.",
  "- DOI: 10.1200/JCO.18.01314.",
  "- Platform: GPL14951, Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip.",
  "",
  "## Object Summaries",
  "",
  md_table(summaries),
  "",
  "## Internal Consistency Checks",
  "",
  md_table(transform(sample_checks, result = vapply(result, fmt_bool, character(1)))),
  "",
  "Note: expression and clinical samples match, but the total files are not in identical order. Analysis scripts should reorder clinical metadata to expression columns before modeling.",
  "",
  "## Expression Summary",
  "",
  paste("- Total expression range:", paste(round(expr_range, 4), collapse = " to ")),
  paste("- Total expression contains NA:", fmt_bool(anyNA(expr_total))),
  "",
  "## Treatment Counts",
  "",
  md_table(treatment_counts),
  "",
  "## Molecular Class Counts",
  "",
  md_table(class_counts),
  "",
  "## Survival Field Missingness",
  "",
  md_table(survival_missing),
  "",
  "## Clinical Missingness By Field",
  "",
  md_table(clinical_missing),
  "",
  "## Supplementary Files",
  "",
  paste("-", supplementary_files),
  "",
  if (length(lock_files) > 0) {
    c(
      "## Notes",
      "",
      paste("- Temporary Excel lock file detected:", paste(lock_files, collapse = ", ")),
      "- This file should be ignored by analysis code."
    )
  } else {
    character(0)
  }
)

writeLines(report, file.path(report_dir, "sha_dataset_validation.md"))

cat("Wrote validation report to reports/sha_dataset_validation.md\n")
