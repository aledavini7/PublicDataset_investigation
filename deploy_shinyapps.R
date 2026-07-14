#!/usr/bin/env Rscript

# Deploy from the repository root so shinyapps.io receives app/, scripts/,
# and sha/curated_gene_symbols/. The .rscignore file excludes generated outputs
# and raw working files.

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install rsconnect first with install.packages('rsconnect').")
}

rsconnect::deployApp(
  appDir = ".",
  appPrimaryDoc = "app/app.R",
  appName = "dlbcl-public-dataset-explorer"
)

