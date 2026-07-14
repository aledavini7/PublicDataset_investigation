#!/usr/bin/env Rscript

# Batch boxplot association analyses using the shared app-ready functions.

source("scripts/analysis_functions.R")

output_dir <- "results/gene_boxplots"

summary_all <- list()
pairwise_all <- list()
sample_groups_all <- list()
output_index <- list()
skipped <- list()

for (dataset_name in names(analysis_datasets())) {
  for (gene in default_genes) {
    for (comparison_name in names(boxplot_comparisons)) {
      result <- tryCatch(
        run_boxplot_analysis(
          dataset_name = dataset_name,
          gene = gene,
          comparison_name = comparison_name,
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
          comparison = comparison_name,
          reason = result$reason
        )
        next
      }

      summary_all[[length(summary_all) + 1]] <- result$summary
      pairwise_all[[length(pairwise_all) + 1]] <- result$pairwise
      sample_groups_all[[length(sample_groups_all) + 1]] <- result$sample_groups
      output_index[[length(output_index) + 1]] <- tibble(
        dataset = dataset_name,
        gene = gene,
        comparison = comparison_name,
        plot_pdf = result$output_paths[["plot_pdf"]],
        plot_png = result$output_paths[["plot_png"]],
        summary_csv = result$output_paths[["summary_csv"]],
        pairwise_csv = result$output_paths[["pairwise_csv"]],
        sample_groups_csv = result$output_paths[["sample_groups_csv"]]
      )
    }
  }
}

table_dir <- file.path(output_dir, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

if (length(summary_all) > 0) {
  write_csv(bind_rows(summary_all), file.path(table_dir, "gene_boxplot_summary.csv"))
  message("Wrote aggregate summary table to ", file.path(table_dir, "gene_boxplot_summary.csv"))
}

if (length(pairwise_all) > 0) {
  write_csv(bind_rows(pairwise_all), file.path(table_dir, "gene_boxplot_pairwise.csv"))
  message("Wrote aggregate pairwise table to ", file.path(table_dir, "gene_boxplot_pairwise.csv"))
}

if (length(sample_groups_all) > 0) {
  write_csv(bind_rows(sample_groups_all), file.path(table_dir, "gene_boxplot_sample_groups.csv"))
  message("Wrote aggregate sample groups to ", file.path(table_dir, "gene_boxplot_sample_groups.csv"))
}

if (length(output_index) > 0) {
  write_csv(bind_rows(output_index), file.path(table_dir, "gene_boxplot_output_index.csv"))
  message("Wrote output index to ", file.path(table_dir, "gene_boxplot_output_index.csv"))
}

skipped_table <- if (length(skipped) > 0) {
  bind_rows(skipped)
} else {
  tibble(dataset = character(), gene = character(), comparison = character(), reason = character())
}
write_csv(skipped_table, file.path(table_dir, "gene_boxplot_skipped.csv"))
message("Wrote skipped analysis table to ", file.path(table_dir, "gene_boxplot_skipped.csv"))

message("Wrote structured plots to ", file.path(output_dir, "plots/<gene>/<dataset>/<comparison>"))
message("Wrote structured tables to ", file.path(table_dir, "<gene>/<dataset>/<comparison>"))
