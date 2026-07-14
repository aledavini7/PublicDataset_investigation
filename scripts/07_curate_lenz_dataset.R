#!/usr/bin/env Rscript

# Curate Lenz GSE10846 expression and clinical data for the multi-dataset app.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

raw_dir <- "external_datasets_raw/lenz"
output_dir <- "curated_datasets/lenz"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expression_raw <- readRDS(file.path(raw_dir, "lenz_total_expression.rds"))
clinical_raw <- readRDS(file.path(raw_dir, "lenz_sha_tot.rds"))

clinical_raw$geo_accession <- as.character(clinical_raw$geo_accession)

if (!identical(colnames(expression_raw), clinical_raw$geo_accession)) {
  if (!all(colnames(expression_raw) %in% clinical_raw$geo_accession)) {
    missing_samples <- setdiff(colnames(expression_raw), clinical_raw$geo_accession)
    stop("Clinical metadata missing sample(s): ", paste(missing_samples, collapse = ", "))
  }
  clinical_raw <- clinical_raw[match(colnames(expression_raw), clinical_raw$geo_accession), , drop = FALSE]
}

if (!identical(colnames(expression_raw), clinical_raw$geo_accession)) {
  stop("Could not align Lenz clinical metadata to expression columns.")
}

ambiguous_rows <- grepl("///", rownames(expression_raw), fixed = TRUE)
excluded_ambiguous_gene_rows <- tibble(
  row_label = rownames(expression_raw)[ambiguous_rows],
  first_symbol = trimws(sub("///.*", "", row_label)),
  n_symbols = vapply(strsplit(row_label, "///", fixed = TRUE), length, integer(1))
)

expression <- expression_raw[!ambiguous_rows, , drop = FALSE]

if (any(duplicated(rownames(expression)))) {
  duplicated_symbols <- unique(rownames(expression)[duplicated(rownames(expression))])
  stop("Duplicated gene symbols after filtering: ", paste(head(duplicated_symbols, 20), collapse = ", "))
}

clean_missing <- function(x) {
  ifelse(is.na(x) | x == "" | toupper(as.character(x)) == "NA", NA, as.character(x))
}

clinical <- clinical_raw %>%
  transmute(
    Sample_ID = geo_accession,
    GEO_accession = geo_accession,
    Title = title,
    Individual = clean_missing(`Individual:ch1`),
    Disease_state = clean_missing(`Disease state:ch1`),
    Tissue = clean_missing(`Tissue:ch1`),
    Treatment = Treatment,
    COO_class = case_when(
      DLBCL_subtype == "ABC DLBCL" ~ "ABC",
      DLBCL_subtype == "GCB DLBCL" ~ "GCB",
      DLBCL_subtype == "Unclassified DLBCL" ~ "UNC",
      TRUE ~ NA_character_
    ),
    DLBCL_subtype_original = DLBCL_subtype,
    OS_months = OS_y * 12,
    OS_event = as.integer(Status_OS),
    OS_years = OS_y,
    Stage = clean_missing(Stage),
    Age = suppressWarnings(as.numeric(Age)),
    Gender = clean_missing(Gender),
    Clinical_info = `Clinical info:ch1`
  )

if (!identical(colnames(expression), clinical$Sample_ID)) {
  stop("Curated expression and clinical sample IDs are not aligned.")
}

dataset_info <- list(
  dataset_id = "lenz_gse10846",
  display_name = "Lenz GSE10846",
  reference_label = "Lenz et al. GSE10846",
  geo_accession = "GSE10846",
  platform = "GPL570 Affymetrix Human Genome U133 Plus 2.0 Array",
  pubmed_ids = c("19038878", "21546504"),
  source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846",
  expression_type = "microarray_log2_mas5_signal",
  gene_id_type = "gene_symbol",
  curation_notes = c(
    "Original expression rows are gene-symbol-like labels.",
    "Rows containing ambiguous multi-symbol labels with /// were excluded from direct app search.",
    "Only overall survival is available in the local/public metadata inspected.",
    "OS_event is coded as 1 = dead/event, 0 = alive/censored."
  ),
  n_samples = ncol(expression),
  n_genes = nrow(expression),
  excluded_ambiguous_gene_rows = nrow(excluded_ambiguous_gene_rows),
  cohorts = list(
    total = list(label = "Total", filter_column = NULL, filter_value = NULL, n = nrow(clinical)),
    chop = list(label = "CHOP", filter_column = "Treatment", filter_value = "CHOP", n = sum(clinical$Treatment == "CHOP")),
    rchop = list(label = "R-CHOP", filter_column = "Treatment", filter_value = "R-CHOP", n = sum(clinical$Treatment == "R-CHOP"))
  ),
  outcomes = list(
    OS = list(time = "OS_years", event = "OS_event", label = "Overall survival", xlab = "Time (years)")
  ),
  stratifications = list(
    treatment = list(column = "Treatment", label = "Treatment", available_in = "total"),
    coo_class = list(column = "COO_class", label = "COO class"),
    stage = list(column = "Stage", label = "Stage"),
    gender = list(column = "Gender", label = "Gender")
  ),
  unavailable = c(
    "PFS",
    "MYC_rearrangement",
    "BCL2_rearrangement",
    "BCL6_rearrangement",
    "hit_status",
    "double_expressor",
    "IHC"
  )
)

saveRDS(expression, file.path(output_dir, "expression.rds"))
saveRDS(clinical, file.path(output_dir, "clinical.rds"))
saveRDS(dataset_info, file.path(output_dir, "dataset_info.rds"))
write_csv(excluded_ambiguous_gene_rows, file.path(output_dir, "excluded_ambiguous_gene_rows.csv"))

summary_table <- tibble(
  metric = c(
    "raw_expression_rows",
    "curated_expression_rows",
    "excluded_ambiguous_rows",
    "samples",
    "chop_samples",
    "rchop_samples",
    "os_events",
    "os_censored"
  ),
  value = c(
    nrow(expression_raw),
    nrow(expression),
    nrow(excluded_ambiguous_gene_rows),
    nrow(clinical),
    sum(clinical$Treatment == "CHOP"),
    sum(clinical$Treatment == "R-CHOP"),
    sum(clinical$OS_event == 1),
    sum(clinical$OS_event == 0)
  )
)
write_csv(summary_table, file.path(output_dir, "curation_summary.csv"))

cat("Wrote curated Lenz dataset to ", output_dir, "\n", sep = "")
cat("Curated expression: ", nrow(expression), " genes x ", ncol(expression), " samples\n", sep = "")
cat("Excluded ambiguous rows: ", nrow(excluded_ambiguous_gene_rows), "\n", sep = "")
