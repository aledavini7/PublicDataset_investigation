#!/usr/bin/env Rscript

# Curate Chapuy GSE98588 expression and clinical data for the multi-dataset app.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

raw_dir <- "external_datasets_raw/chapuy"
output_dir <- "curated_datasets/chapuy"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expression_raw <- readRDS(file.path(raw_dir, "gep_expression_geo_GSE98588.rds"))
clinical_raw <- readRDS(file.path(raw_dir, "clin_137_gep.rds"))

expression <- as.matrix(expression_raw)
mode(expression) <- "numeric"

clinical_raw$pair_id <- as.character(clinical_raw$pair_id)

if (!identical(colnames(expression), clinical_raw$pair_id)) {
  if (!all(colnames(expression) %in% clinical_raw$pair_id)) {
    missing_samples <- setdiff(colnames(expression), clinical_raw$pair_id)
    stop("Clinical metadata missing sample(s): ", paste(missing_samples, collapse = ", "))
  }
  clinical_raw <- clinical_raw[match(colnames(expression), clinical_raw$pair_id), , drop = FALSE]
}

if (!identical(colnames(expression), clinical_raw$pair_id)) {
  stop("Could not align Chapuy clinical metadata to expression columns.")
}

if (any(duplicated(rownames(expression)))) {
  duplicated_symbols <- unique(rownames(expression)[duplicated(rownames(expression))])
  stop("Duplicated gene symbols: ", paste(head(duplicated_symbols, 20), collapse = ", "))
}

if (any(grepl("///", rownames(expression), fixed = TRUE))) {
  stop("Unexpected ambiguous gene-symbol labels containing ///.")
}

clean_missing <- function(x) {
  x <- as.character(x)
  x[x == "" | tolower(x) %in% c("na", "n/a", "nan")] <- NA_character_
  x
}

as_number <- function(x) {
  suppressWarnings(as.numeric(clean_missing(x)))
}

as_event <- function(x) {
  suppressWarnings(as.integer(as_number(x)))
}

standard_high_low <- function(x) {
  x <- clean_missing(x)
  dplyr::case_when(
    toupper(x) == "HIGH" ~ "high",
    toupper(x) == "LOW" ~ "low",
    TRUE ~ NA_character_
  )
}

standard_yes_no <- function(x) {
  x <- clean_missing(x)
  dplyr::case_when(
    toupper(x) %in% c("YES", "Y") ~ "yes",
    toupper(x) %in% c("NO", "N") ~ "no",
    TRUE ~ NA_character_
  )
}

clinical <- clinical_raw %>%
  transmute(
    Sample_ID = pair_id,
    Individual_ID = clean_missing(individual_id),
    GEP_ID = clean_missing(GEPID1),
    Source_cohort = clean_missing(Cohort),
    COO_class = case_when(
      COO_byGEP == "ABC" ~ "ABC",
      COO_byGEP == "GCB" ~ "GCB",
      COO_byGEP == "Unclassified" ~ "UNC",
      TRUE ~ NA_character_
    ),
    COO_any = case_when(
      `Any_COO (GEP+nanostring)` == "ABC" ~ "ABC",
      `Any_COO (GEP+nanostring)` == "GCB" ~ "GCB",
      `Any_COO (GEP+nanostring)` == "Unclassified" ~ "UNC",
      TRUE ~ NA_character_
    ),
    CCC_class = clean_missing(CCC),
    Genomic_complexity = case_when(
      `complex-clean` == "complex" ~ "complex",
      `complex-clean` == "clean" ~ "clean",
      TRUE ~ NA_character_
    ),
    Chapuy_cluster = if_else(!is.na(as_number(CLUSTER)), paste0("C", as.integer(as_number(CLUSTER))), NA_character_),
    Gender = clean_missing(Gender),
    Age = as_number(`Age-at first diagnosis`),
    RCHOP_like = standard_yes_no(`R-CHOP-like Chemo`),
    PFS_months = as_number(PFS),
    PFS_event = as_event(PFS_STAT),
    OS_months = as_number(OS),
    OS_event = as_event(OS_STAT),
    IPI = as_number(IPI),
    IPI_age = as_event(IPI_AGE),
    IPI_LDH = as_event(IPI_LDH),
    IPI_ECOG = as_event(IPI_ECOG),
    IPI_stage = as_event(IPI_STAGE),
    IPI_extranodal = as_event(IPI_EXBM),
    CNS_involvement = as_event(`ob24z (CNS involvment)`),
    Testicular_involvement = as_event(`hodz (Testicular invovlement)`),
    Ploidy = as_number(ploidy_absolute_reviewed),
    Purity = as_number(purity_absolute_reviewed),
    Number_mutations = as_number(numberOfMutations),
    Number_CNA = as_number(numberOfCNAs),
    Number_rearrangements = as_number(numberOfChromosomalRearrangements),
    Number_driver_mutations = as_number(numberOfDriver_Mutations),
    Number_driver_SCNA = as_number(numberOfDriver_SCNAs),
    Number_driver_SV = as_number(numberOfDriver_SVs),
    TP53_status = case_when(
      TP53 == "Mutated" ~ "mutated",
      TP53 == "Not Mutated" ~ "not-mutated",
      TRUE ~ NA_character_
    ),
    six_q_del = case_when(
      `6qdel` == "DEL" ~ "deleted",
      `6qdel` == "norm" ~ "normal",
      TRUE ~ NA_character_
    ),
    MYC_RNA = standard_high_low(MYC),
    MYC_RNA_quartile = standard_high_low(MYC_q),
    BCL2_RNA = standard_high_low(BCL2),
    BCL2_RNA_quartile = standard_high_low(BCL2_q),
    MTOR_RNA = standard_high_low(MTOR),
    MTOR_RNA_quartile = standard_high_low(MTOR_q),
    PIK3CD_RNA = standard_high_low(PIK3CD),
    PIK3CD_RNA_quartile = standard_high_low(PIK3CD_q),
    NFKBIE_RNA = standard_high_low(NFKBIE),
    NFKBIE_RNA_quartile = standard_high_low(NFKBIE_q),
    MYC_rearrangement = case_when(
      MYC_rearr == "YES" ~ "rearranged",
      MYC_rearr == "NO" ~ "not-rearranged",
      TRUE ~ NA_character_
    ),
    MYC_expression_reported = as_number(MYC_exp),
    PARP1_expression_reported = as_number(PARP1_exp)
  )

if (!identical(colnames(expression), clinical$Sample_ID)) {
  stop("Curated expression and clinical sample IDs are not aligned.")
}

dataset_info <- list(
  dataset_id = "chapuy_gse98588",
  display_name = "Chapuy GSE98588",
  reference_label = "Chapuy et al. Cancer Cell 2018",
  geo_accession = "GSE98588",
  platform = "GPL23432 Affymetrix Human Genome U133 Plus 2.0 Array, Brainarray ENSG v18 CDF",
  pubmed_ids = c("29713087"),
  source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98588",
  expression_type = "microarray_log2_brainarray_gene_signal",
  gene_id_type = "gene_symbol",
  curation_notes = c(
    "Expression columns match clinical pair_id exactly.",
    "Expression rows are already gene symbols with no duplicates and no /// ambiguous labels.",
    "OS and PFS are represented in months in the provided clinical file.",
    "PFS_event and OS_event are coded as 1 = event, 0 = censored."
  ),
  n_samples = ncol(expression),
  n_genes = nrow(expression),
  cohorts = list(
    total = list(label = "Total", filter_column = NULL, filter_value = NULL, n = nrow(clinical)),
    rchop_like = list(label = "R-CHOP-like", filter_column = "RCHOP_like", filter_value = "yes", n = sum(clinical$RCHOP_like == "yes", na.rm = TRUE)),
    non_rchop_like = list(label = "Non-R-CHOP-like", filter_column = "RCHOP_like", filter_value = "no", n = sum(clinical$RCHOP_like == "no", na.rm = TRUE))
  ),
  outcomes = list(
    OS = list(time = "OS_months", event = "OS_event", label = "Overall survival", xlab = "Time (months)"),
    PFS = list(time = "PFS_months", event = "PFS_event", label = "Progression-free survival", xlab = "Time (months)")
  )
)

saveRDS(expression, file.path(output_dir, "expression.rds"))
saveRDS(clinical, file.path(output_dir, "clinical.rds"))
saveRDS(dataset_info, file.path(output_dir, "dataset_info.rds"))

summary_table <- tibble(
  metric = c(
    "expression_rows",
    "samples",
    "cohort_total",
    "cohort_rchop_like",
    "cohort_non_rchop_like",
    "pfs_complete",
    "os_complete",
    "pfs_events",
    "os_events",
    "coo_abc",
    "coo_gcb",
    "coo_unc",
    "myc_rearranged"
  ),
  value = c(
    nrow(expression),
    nrow(clinical),
    nrow(clinical),
    sum(clinical$RCHOP_like == "yes", na.rm = TRUE),
    sum(clinical$RCHOP_like == "no", na.rm = TRUE),
    sum(!is.na(clinical$PFS_months) & !is.na(clinical$PFS_event)),
    sum(!is.na(clinical$OS_months) & !is.na(clinical$OS_event)),
    sum(clinical$PFS_event == 1, na.rm = TRUE),
    sum(clinical$OS_event == 1, na.rm = TRUE),
    sum(clinical$COO_class == "ABC", na.rm = TRUE),
    sum(clinical$COO_class == "GCB", na.rm = TRUE),
    sum(clinical$COO_class == "UNC", na.rm = TRUE),
    sum(clinical$MYC_rearrangement == "rearranged", na.rm = TRUE)
  )
)
write_csv(summary_table, file.path(output_dir, "curation_summary.csv"))

cat("Wrote curated Chapuy dataset to ", output_dir, "\n", sep = "")
cat("Curated expression: ", nrow(expression), " genes x ", ncol(expression), " samples\n", sep = "")
cat("Complete PFS samples: ", sum(!is.na(clinical$PFS_months) & !is.na(clinical$PFS_event)), "\n", sep = "")
cat("Complete OS samples: ", sum(!is.na(clinical$OS_months) & !is.na(clinical$OS_event)), "\n", sep = "")
