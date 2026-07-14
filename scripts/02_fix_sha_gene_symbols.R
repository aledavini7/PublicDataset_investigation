#!/usr/bin/env Rscript

# Create corrected copies of the curated Sha/REMoDLB RDS files with obvious
# legacy/space-containing gene names converted to current symbol-style names.

input_dir <- "sha"
output_dir <- file.path(input_dir, "curated_gene_symbols")
report_dir <- "reports"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

expression_files <- c(
  sha_total_expression = file.path(input_dir, "sha_total_expression.rds"),
  sha_rchop_expression = file.path(input_dir, "sha_rchop_expression.rds"),
  sha_rbchop_expression = file.path(input_dir, "sha_rbchop_expression.rds")
)

clinical_files <- c(
  clin_sha_tot = file.path(input_dir, "clin_sha_tot.rds"),
  clin_sha_rchop = file.path(input_dir, "clin_sha_rchop.rds"),
  clin_sha_rbchop = file.path(input_dir, "clin_sha_rbchop.rds")
)

gene_symbol_map <- data.frame(
  old_symbol = c(
    "March 1", "March 2", "March 3", "March 4", "March 5", "March 6",
    "March 7", "March 8", "March 9", "March 10", "March 11",
    "Selenoprotein 15",
    "Septin 1", "Septin 2", "Septin 3", "Septin 4", "Septin 5",
    "Septin 6", "Septin 7", "Septin 9", "Septin 10", "Septin 11",
    "Septin 12", "Septin 13", "Septin 14"
  ),
  new_symbol = c(
    "MARCHF1", "MARCHF2", "MARCHF3", "MARCHF4", "MARCHF5", "MARCHF6",
    "MARCHF7", "MARCHF8", "MARCHF9", "MARCHF10", "MARCHF11",
    "SELENOF",
    "SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN4", "SEPTIN5",
    "SEPTIN6", "SEPTIN7", "SEPTIN9", "SEPTIN10", "SEPTIN11",
    "SEPTIN12", "SEPTIN7P2", "SEPTIN14"
  ),
  reason = c(
    rep("Legacy MARCH family name converted to current MARCHF-style symbol", 11),
    "Legacy selenoprotein name converted to current symbol",
    rep("Legacy septin name converted to current SEPTIN-style symbol", 11),
    "Legacy SEPT13 alias maps to SEPTIN7P2 in local org.Hs.eg.db annotation",
    "Legacy septin name converted to current SEPTIN-style symbol"
  ),
  stringsAsFactors = FALSE
)

stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing expected file(s): ", paste(missing, collapse = ", "))
  }
}

fix_expression_rownames <- function(path, mapping) {
  expression <- readRDS(path)
  old_rownames <- rownames(expression)
  new_rownames <- old_rownames
  matched <- match(new_rownames, mapping$old_symbol)
  to_replace <- !is.na(matched)
  new_rownames[to_replace] <- mapping$new_symbol[matched[to_replace]]

  if (anyDuplicated(new_rownames) != 0) {
    duplicated_symbols <- unique(new_rownames[duplicated(new_rownames)])
    stop(
      "Correction would create duplicated row names in ", path, ": ",
      paste(duplicated_symbols, collapse = ", ")
    )
  }

  rownames(expression) <- new_rownames
  expression
}

stop_if_missing(c(expression_files, clinical_files))

correction_summary <- data.frame()

for (name in names(expression_files)) {
  source_path <- expression_files[[name]]
  corrected <- fix_expression_rownames(source_path, gene_symbol_map)
  output_path <- file.path(output_dir, paste0(name, ".rds"))
  saveRDS(corrected, output_path)

  correction_summary <- rbind(
    correction_summary,
    data.frame(
      file = basename(output_path),
      rows = nrow(corrected),
      columns = ncol(corrected),
      duplicated_rows = anyDuplicated(rownames(corrected)),
      space_containing_rows = sum(grepl("[[:space:]]", rownames(corrected))),
      stringsAsFactors = FALSE
    )
  )
}

for (name in names(clinical_files)) {
  clinical <- readRDS(clinical_files[[name]])
  saveRDS(clinical, file.path(output_dir, paste0(name, ".rds")))
}

write.csv(
  gene_symbol_map,
  file.path(report_dir, "sha_gene_symbol_corrections.csv"),
  row.names = FALSE,
  quote = TRUE
)

report <- c(
  "# Sha Gene Symbol Correction Report",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Output",
  "",
  paste("- Corrected RDS folder:", output_dir),
  "- Original RDS files were not modified.",
  "- Clinical RDS files were copied unchanged into the corrected folder for convenience.",
  "",
  "## Correction Summary",
  "",
  paste(capture.output(print(correction_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Mapping",
  "",
  paste(capture.output(print(gene_symbol_map, row.names = FALSE)), collapse = "\n"),
  "",
  "## Notes",
  "",
  "- These corrections target the 25 obvious space-containing legacy names found in the expression row names.",
  "- `Septin 13` is handled conservatively using the local `org.Hs.eg.db` alias mapping from `SEPT13` to `SEPTIN7P2`.",
  "- Broader alias modernization across all 20,816 rows should be treated as a separate step because many old symbols are biologically valid historical aliases and may require careful one-to-many handling."
)

writeLines(report, file.path(report_dir, "sha_gene_symbol_correction_report.md"))

cat("Wrote corrected RDS files to ", output_dir, "\n", sep = "")
cat("Wrote mapping to reports/sha_gene_symbol_corrections.csv\n")
cat("Wrote report to reports/sha_gene_symbol_correction_report.md\n")
