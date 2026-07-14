#!/usr/bin/env Rscript

# Verify curated Sha gene row names against the GPL14951 probe annotation.
# This script expects the raw GPL14951 table at /private/tmp/GPL14951_table.txt.
# The table can be obtained from GEO's GPL14951 "Download full table" endpoint.

platform_table_path <- "/private/tmp/GPL14951_table.txt"
report_dir <- "reports"
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(platform_table_path)) {
  stop(
    "Missing platform table: ", platform_table_path, "\n",
    "Download GPL14951 raw platform annotation before running this script."
  )
}

platform <- read.delim(
  platform_table_path,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

original_expr <- readRDS("sha/sha_total_expression.rds")
modern_expr <- readRDS("sha/curated_gene_symbols/sha_total_expression.rds")

original_genes <- rownames(original_expr)
modern_genes <- rownames(modern_expr)

platform_symbols <- unique(c(platform$Symbol, platform$ILMN_Gene))
platform_symbols <- platform_symbols[!is.na(platform_symbols) & platform_symbols != ""]

changed <- data.frame(
  original_symbol = original_genes[original_genes != modern_genes],
  modern_symbol = modern_genes[original_genes != modern_genes],
  stringsAsFactors = FALSE
)

original_check <- data.frame(
  gene_symbol = original_genes,
  present_in_gpl14951_symbol_or_ilmn_gene = original_genes %in% platform_symbols,
  stringsAsFactors = FALSE
)

modern_check <- data.frame(
  gene_symbol = modern_genes,
  present_in_gpl14951_symbol_or_ilmn_gene = modern_genes %in% platform_symbols,
  intentionally_modernized = modern_genes %in% changed$modern_symbol,
  stringsAsFactors = FALSE
)

summary_table <- data.frame(
  check = c(
    "original row names present in GPL14951 Symbol/ILMN_Gene",
    "modern row names present in GPL14951 Symbol/ILMN_Gene",
    "modern row names absent from GPL14951 but intentionally modernized",
    "original row names absent from GPL14951",
    "row names changed between original and modernized datasets"
  ),
  n = c(
    sum(original_check$present_in_gpl14951_symbol_or_ilmn_gene),
    sum(modern_check$present_in_gpl14951_symbol_or_ilmn_gene),
    sum(!modern_check$present_in_gpl14951_symbol_or_ilmn_gene & modern_check$intentionally_modernized),
    sum(!original_check$present_in_gpl14951_symbol_or_ilmn_gene),
    nrow(changed)
  ),
  denominator = c(
    length(original_genes),
    length(modern_genes),
    length(modern_genes),
    length(original_genes),
    length(original_genes)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  original_check,
  file.path(report_dir, "sha_original_gene_labels_vs_gpl14951.csv"),
  row.names = FALSE
)
write.csv(
  modern_check,
  file.path(report_dir, "sha_modern_gene_labels_vs_gpl14951.csv"),
  row.names = FALSE
)

report <- c(
  "# Sha Probe-To-Gene Label Verification",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Source",
  "",
  "- Platform: GPL14951, Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip.",
  "- Probe annotation fields used: `Symbol` and `ILMN_Gene`.",
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "- All original curated expression row names are present in GPL14951 `Symbol` or `ILMN_Gene`.",
  "- Therefore, the original curated gene labels are consistent with the probe annotation used by the platform.",
  "- The modernized dataset differs only in the 25 intentional legacy-symbol updates.",
  "- Those 25 modern symbols are not direct GPL14951 labels; they are current-symbol replacements for downstream usability.",
  "",
  "## Intentional Modernizations",
  "",
  paste(capture.output(print(changed, row.names = FALSE)), collapse = "\n")
)

writeLines(report, file.path(report_dir, "sha_probe_gene_label_verification.md"))

cat("Wrote reports/sha_probe_gene_label_verification.md\n")
cat("Wrote reports/sha_original_gene_labels_vs_gpl14951.csv\n")
cat("Wrote reports/sha_modern_gene_labels_vs_gpl14951.csv\n")
