#!/usr/bin/env Rscript

# Curate Barrans GSE32918 expression and clinical data for the multi-dataset app.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

raw_dir <- "external_datasets_raw/barrans"
output_dir <- "curated_datasets/barrans"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

matrix_raw <- read.delim(
  file.path(raw_dir, "barrans_matrix.tsv"),
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata_raw <- read.delim(
  file.path(raw_dir, "barrans_metadata.tsv"),
  header = FALSE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation_raw <- read.delim(
  file.path(raw_dir, "GPL8432-11703.tsv"),
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

probe_ids <- matrix_raw$ID_REF
expression_probe <- as.matrix(matrix_raw[, -1, drop = FALSE])
mode(expression_probe) <- "numeric"
rownames(expression_probe) <- probe_ids

metadata_rows <- metadata_raw[[1]]
metadata_values <- metadata_raw[, -1, drop = FALSE]

get_metadata_row <- function(row_name) {
  idx <- which(metadata_rows == row_name)[1]
  if (is.na(idx)) {
    stop("Metadata row not found: ", row_name)
  }
  as.character(unlist(metadata_values[idx, ], use.names = FALSE))
}

get_characteristic <- function(prefix) {
  idx <- which(metadata_rows == "!Sample_characteristics_ch1" & grepl(paste0("^", prefix, ":"), metadata_raw[[2]]))[1]
  if (is.na(idx)) {
    stop("Metadata characteristic not found: ", prefix)
  }
  sub("^[^:]+: ?", "", as.character(unlist(metadata_values[idx, ], use.names = FALSE)))
}

array_ids <- colnames(expression_probe)
geo_accessions <- get_metadata_row("!Sample_geo_accession")
sample_titles <- get_metadata_row("!Sample_title")

if (!identical(array_ids, geo_accessions)) {
  stop("Expression columns do not match metadata GEO accessions.")
}

sample_ids <- sub("_Rep[0-9]+$", "", sample_titles)
unique_sample_ids <- unique(sample_ids)

annotation <- annotation_raw %>%
  transmute(
    probe_id = ID,
    symbol = trimws(Symbol),
    ilmn_gene = trimws(ILMN_Gene),
    entrez_gene_id = suppressWarnings(as.integer(Entrez_Gene_ID)),
    refseq_id = RefSeq_ID,
    probe_type = Probe_Type
  ) %>%
  mutate(symbol = if_else(is.na(symbol) | symbol == "", NA_character_, symbol))

annotation <- annotation[match(probe_ids, annotation$probe_id), , drop = FALSE]

if (!identical(annotation$probe_id, probe_ids)) {
  stop("Could not align platform annotation to expression probes.")
}

keep_probes <- !is.na(annotation$symbol)
expression_probe <- expression_probe[keep_probes, , drop = FALSE]
annotation <- annotation[keep_probes, , drop = FALSE]

expression_sample <- sapply(unique_sample_ids, function(sample_id) {
  cols <- which(sample_ids == sample_id)
  if (length(cols) == 1) {
    expression_probe[, cols]
  } else {
    rowMeans(expression_probe[, cols, drop = FALSE], na.rm = TRUE)
  }
})
expression_sample <- as.matrix(expression_sample)
rownames(expression_sample) <- rownames(expression_probe)

probe_stats <- tibble(
  probe_id = rownames(expression_sample),
  symbol = annotation$symbol,
  median_expression = apply(expression_sample, 1, median, na.rm = TRUE),
  mad_expression = apply(expression_sample, 1, mad, na.rm = TRUE)
) %>%
  mutate(probe_index = row_number()) %>%
  group_by(symbol) %>%
  arrange(desc(mad_expression), desc(median_expression), probe_id, .by_group = TRUE) %>%
  mutate(selected_for_gene = row_number() == 1) %>%
  ungroup()

selected_probes <- probe_stats %>%
  filter(selected_for_gene) %>%
  arrange(symbol)

expression <- expression_sample[selected_probes$probe_index, , drop = FALSE]
rownames(expression) <- selected_probes$symbol

if (any(duplicated(rownames(expression)))) {
  stop("Duplicated gene symbols after probe selection.")
}

array_clinical <- tibble(
  Sample_ID = sample_ids,
  Array_GEO_accession = geo_accessions,
  Array_title = sample_titles,
  Predicted_class_raw = get_characteristic("predicted class"),
  Class_confidence = suppressWarnings(as.numeric(get_characteristic("class confidence"))),
  Age = suppressWarnings(as.numeric(get_characteristic("age"))),
  Gender = get_characteristic("Sex"),
  Treatment_raw = get_characteristic("treated status"),
  Follow_up_status = get_characteristic("follow-up status"),
  OS_years_raw = suppressWarnings(as.numeric(get_characteristic("follow-up years")))
)

clinical_conflicts <- array_clinical %>%
  group_by(Sample_ID) %>%
  summarise(
    across(
      c(Predicted_class_raw, Class_confidence, Age, Gender, Treatment_raw, Follow_up_status, OS_years_raw),
      ~ n_distinct(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(if_any(-Sample_ID, ~ .x > 1))

if (nrow(clinical_conflicts) > 0) {
  stop("Conflicting clinical metadata across technical replicates.")
}

clinical <- array_clinical %>%
  group_by(Sample_ID) %>%
  summarise(
    GEO_accessions = paste(Array_GEO_accession, collapse = ";"),
    Array_titles = paste(Array_title, collapse = ";"),
    N_arrays = n(),
    Predicted_class_raw = first(Predicted_class_raw),
    Class_confidence = first(Class_confidence),
    Age = first(Age),
    Gender = first(Gender),
    Treatment_raw = first(Treatment_raw),
    Follow_up_status = first(Follow_up_status),
    OS_years_raw = first(OS_years_raw),
    .groups = "drop"
  ) %>%
  mutate(
    COO_class = case_when(
      Predicted_class_raw == "ABC" ~ "ABC",
      Predicted_class_raw == "GCB" ~ "GCB",
      Predicted_class_raw == "TypeIII" ~ "UNC",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      Treatment_raw == "1" ~ "R-CHOP",
      Treatment_raw == "0" ~ "not-curative-intent",
      TRUE ~ NA_character_
    ),
    OS_years = if_else(OS_years_raw >= 0, OS_years_raw, NA_real_),
    OS_event = case_when(
      is.na(OS_years) ~ NA_integer_,
      Follow_up_status == "Dead" ~ 1L,
      Follow_up_status == "Alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(
    Sample_ID,
    GEO_accessions,
    Array_titles,
    N_arrays,
    COO_class,
    Predicted_class_raw,
    Class_confidence,
    Age,
    Gender,
    Treatment,
    Follow_up_status,
    OS_years,
    OS_event,
    OS_years_raw
  )

if (!identical(colnames(expression), clinical$Sample_ID)) {
  clinical <- clinical[match(colnames(expression), clinical$Sample_ID), , drop = FALSE]
}

if (!identical(colnames(expression), clinical$Sample_ID)) {
  stop("Curated expression and clinical sample IDs are not aligned.")
}

selected_probe_annotation <- selected_probes %>%
  left_join(annotation, by = c("probe_id", "symbol"))

dataset_info <- list(
  dataset_id = "barrans_gse32918",
  display_name = "Barrans GSE32918",
  reference_label = "Barrans et al. GSE32918",
  geo_accession = "GSE32918",
  platform = "GPL8432 Illumina HumanRef-8 WG-DASL v3.0",
  pubmed_ids = c("22970711", "24875472"),
  source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32918",
  expression_type = "Illumina DASL FFPE expression, technical replicates averaged",
  gene_id_type = "gene_symbol",
  curation_notes = c(
    "Expression arrays were collapsed from 249 technical arrays to 172 patient-level samples by averaging replicate arrays.",
    "Illumina probes were mapped to gene symbols using GPL8432-11703.tsv.",
    "For genes represented by multiple probes, the probe with highest MAD across patient-level samples was selected.",
    "Only OS is available in the inspected GEO metadata.",
    "One raw negative follow-up time was set to missing for OS analyses.",
    "OS_event is coded as 1 = dead/event, 0 = alive/censored."
  ),
  n_arrays = ncol(expression_probe),
  n_samples = ncol(expression),
  n_genes = nrow(expression),
  n_input_probes = nrow(matrix_raw),
  n_symbol_annotated_probes = nrow(expression_sample),
  n_duplicate_symbol_probes = sum(duplicated(annotation$symbol)),
  n_negative_follow_up = sum(array_clinical$OS_years_raw < 0, na.rm = TRUE),
  cohorts = list(
    total = list(label = "Total", filter_column = NULL, filter_value = NULL, n = nrow(clinical)),
    rchop = list(label = "R-CHOP", filter_column = "Treatment", filter_value = "R-CHOP", n = sum(clinical$Treatment == "R-CHOP", na.rm = TRUE)),
    non_curative = list(label = "Not curative intent", filter_column = "Treatment", filter_value = "not-curative-intent", n = sum(clinical$Treatment == "not-curative-intent", na.rm = TRUE))
  ),
  outcomes = list(
    OS = list(time = "OS_years", event = "OS_event", label = "Overall survival", xlab = "Time (years)")
  )
)

saveRDS(expression, file.path(output_dir, "expression.rds"))
saveRDS(clinical, file.path(output_dir, "clinical.rds"))
saveRDS(dataset_info, file.path(output_dir, "dataset_info.rds"))
write_csv(selected_probe_annotation, file.path(output_dir, "selected_probe_annotation.csv"))
write_csv(probe_stats, file.path(output_dir, "probe_annotation_selection.csv"))

summary_table <- tibble(
  metric = c(
    "raw_probe_rows",
    "symbol_annotated_probe_rows",
    "curated_gene_rows",
    "raw_arrays",
    "curated_samples",
    "technical_replicated_samples",
    "rchop_samples",
    "non_curative_samples",
    "os_complete",
    "os_events",
    "negative_follow_up_set_missing",
    "coo_abc",
    "coo_gcb",
    "coo_unc"
  ),
  value = c(
    nrow(matrix_raw),
    nrow(expression_sample),
    nrow(expression),
    ncol(expression_probe),
    nrow(clinical),
    sum(clinical$N_arrays > 1),
    sum(clinical$Treatment == "R-CHOP", na.rm = TRUE),
    sum(clinical$Treatment == "not-curative-intent", na.rm = TRUE),
    sum(!is.na(clinical$OS_years) & !is.na(clinical$OS_event)),
    sum(clinical$OS_event == 1, na.rm = TRUE),
    sum(is.na(clinical$OS_years) & !is.na(clinical$OS_years_raw)),
    sum(clinical$COO_class == "ABC", na.rm = TRUE),
    sum(clinical$COO_class == "GCB", na.rm = TRUE),
    sum(clinical$COO_class == "UNC", na.rm = TRUE)
  )
)
write_csv(summary_table, file.path(output_dir, "curation_summary.csv"))

cat("Wrote curated Barrans dataset to ", output_dir, "\n", sep = "")
cat("Curated expression: ", nrow(expression), " genes x ", ncol(expression), " samples\n", sep = "")
cat("Raw arrays collapsed: ", ncol(expression_probe), " arrays to ", ncol(expression), " samples\n", sep = "")
cat("Complete OS samples: ", sum(!is.na(clinical$OS_years) & !is.na(clinical$OS_event)), "\n", sep = "")
