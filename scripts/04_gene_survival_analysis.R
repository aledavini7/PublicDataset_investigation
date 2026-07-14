#!/usr/bin/env Rscript

# Batch Kaplan-Meier survival analyses using the shared app-ready functions.

source("scripts/analysis_functions.R")

output_dir <- "results/gene_survival"

all_summaries <- list()
all_sample_groups <- list()
all_plot_paths <- list()
skipped <- list()

for (dataset_name in names(analysis_datasets())) {
  for (gene in default_genes) {
    for (grouping in default_grouping_methods) {
      for (analysis_type in names(survival_stratifications)) {
        for (outcome_name in names(survival_outcomes)) {
          result <- tryCatch(
            run_survival_analysis(
              dataset_name = dataset_name,
              gene = gene,
              grouping = grouping,
              outcome_name = outcome_name,
              analysis_type = analysis_type,
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
              analysis_type = analysis_type,
              grouping = grouping,
              outcome = outcome_name,
              reason = result$reason
            )
            next
          }

          all_summaries[[length(all_summaries) + 1]] <- result$summary
          all_sample_groups[[length(all_sample_groups) + 1]] <- result$sample_groups
          all_plot_paths[[length(all_plot_paths) + 1]] <- tibble(
            dataset = dataset_name,
            gene = gene,
            analysis_type = analysis_type,
            grouping = grouping,
            outcome = outcome_name,
            plot_pdf = result$output_paths[["plot_pdf"]],
            plot_png = result$output_paths[["plot_png"]],
            risk_table_pdf = result$output_paths[["risk_table_pdf"]],
            risk_table_png = result$output_paths[["risk_table_png"]],
            summary_csv = result$output_paths[["summary_csv"]],
            sample_groups_csv = result$output_paths[["sample_groups_csv"]]
          )
        }
      }
    }
  }
}

table_dir <- file.path(output_dir, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

summary_table <- bind_rows(all_summaries)
sample_groups <- bind_rows(all_sample_groups)
plot_index <- bind_rows(all_plot_paths)
skipped_table <- if (length(skipped) > 0) {
  bind_rows(skipped)
} else {
  tibble(
    dataset = character(),
    gene = character(),
    analysis_type = character(),
    grouping = character(),
    outcome = character(),
    reason = character()
  )
}

write_csv(summary_table, file.path(table_dir, "gene_survival_summary.csv"))
write_csv(sample_groups, file.path(table_dir, "gene_survival_sample_groups.csv"))
write_csv(plot_index, file.path(table_dir, "gene_survival_output_index.csv"))
write_csv(skipped_table, file.path(table_dir, "gene_survival_skipped.csv"))

cat("Wrote aggregate summary table to ", file.path(table_dir, "gene_survival_summary.csv"), "\n", sep = "")
cat("Wrote aggregate sample groups to ", file.path(table_dir, "gene_survival_sample_groups.csv"), "\n", sep = "")
cat("Wrote output index to ", file.path(table_dir, "gene_survival_output_index.csv"), "\n", sep = "")
cat("Wrote structured plots to ", file.path(output_dir, "plots/<gene>/<dataset>/<analysis_type>"), "\n", sep = "")
cat("Wrote structured tables to ", file.path(table_dir, "<gene>/<dataset>/<analysis_type>"), "\n", sep = "")
