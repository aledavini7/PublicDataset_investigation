#!/usr/bin/env Rscript

# Batch expression-expression correlation analyses using the shared app-ready functions.

source("scripts/analysis_functions.R")

output_dir <- "results/gene_correlations"

summary_all <- list()
sample_values_all <- list()
output_index <- list()
skipped <- list()

for (dataset_name in names(analysis_datasets())) {
  for (gene in default_genes) {
    for (target_gene in default_target_genes) {
      if (gene == target_gene) {
        next
      }

      result <- tryCatch(
        run_correlation_analysis(
          dataset_name = dataset_name,
          gene = gene,
          target_gene = target_gene,
          save_outputs = TRUE,
          output_dir = output_dir
        ),
        error = function(e) {
          list(valid = FALSE, reason = conditionMessage(e))
        }
      )

      if (!isTRUE(result$valid)) {
        skipped[[length(skipped) + 1]] <- tibble(
          dataset = dataset_name,
          gene = gene,
          target_gene = target_gene,
          reason = result$reason
        )
        next
      }

      summary_all[[length(summary_all) + 1]] <- result$summary
      sample_values_all[[length(sample_values_all) + 1]] <- result$sample_values
      output_index[[length(output_index) + 1]] <- tibble(
        dataset = dataset_name,
        gene = gene,
        target_gene = target_gene,
        plot_pdf = result$output_paths[["plot_pdf"]],
        plot_png = result$output_paths[["plot_png"]],
        summary_csv = result$output_paths[["summary_csv"]],
        sample_values_csv = result$output_paths[["sample_values_csv"]]
      )
    }
  }
}

table_dir <- file.path(output_dir, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

if (length(summary_all) > 0) {
  write_csv(bind_rows(summary_all), file.path(table_dir, "gene_expression_correlation_summary.csv"))
  message("Wrote aggregate summary table to ", file.path(table_dir, "gene_expression_correlation_summary.csv"))
}

if (length(sample_values_all) > 0) {
  write_csv(bind_rows(sample_values_all), file.path(table_dir, "gene_expression_correlation_sample_values.csv"))
  message("Wrote aggregate sample values to ", file.path(table_dir, "gene_expression_correlation_sample_values.csv"))
}

if (length(output_index) > 0) {
  write_csv(bind_rows(output_index), file.path(table_dir, "gene_expression_correlation_output_index.csv"))
  message("Wrote output index to ", file.path(table_dir, "gene_expression_correlation_output_index.csv"))
}

skipped_table <- if (length(skipped) > 0) {
  bind_rows(skipped)
} else {
  tibble(dataset = character(), gene = character(), target_gene = character(), reason = character())
}
write_csv(skipped_table, file.path(table_dir, "gene_expression_correlation_skipped.csv"))
message("Wrote skipped analysis table to ", file.path(table_dir, "gene_expression_correlation_skipped.csv"))

message("Wrote structured plots to ", file.path(output_dir, "plots/<gene>/<dataset>/<target_gene>"))
message("Wrote structured tables to ", file.path(table_dir, "<gene>/<dataset>/<target_gene>"))
