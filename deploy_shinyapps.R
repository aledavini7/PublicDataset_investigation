#!/usr/bin/env Rscript

# Deploy from the repository root so shinyapps.io receives only the Shiny app,
# shared analysis functions, and curated public RDS files required at runtime.

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install rsconnect first with install.packages('rsconnect').")
}

app_files <- c(
  "DESCRIPTION",
  "LICENSE",
  "README.md",
  "app/app.R",
  "scripts/analysis_functions.R",
  list.files(
    "sha/curated_gene_symbols",
    pattern = "\\.rds$",
    full.names = TRUE
  ),
  list.files(
    "curated_datasets/lenz",
    pattern = "\\.rds$",
    full.names = TRUE
  ),
  list.files(
    "curated_datasets/chapuy",
    pattern = "\\.rds$",
    full.names = TRUE
  ),
  list.files(
    "curated_datasets/barrans",
    pattern = "\\.rds$",
    full.names = TRUE
  )
)

missing_files <- app_files[!file.exists(app_files)]
if (length(missing_files) > 0) {
  stop("Missing deployment file(s): ", paste(missing_files, collapse = ", "))
}

message("Deploying ", length(app_files), " files.")

rsconnect::deployApp(
  appDir = ".",
  appFiles = app_files,
  appPrimaryDoc = "app/app.R",
  appName = "dlbcl-public-dataset-explorer"
)
