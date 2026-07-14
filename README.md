# DLBCL Public Dataset Explorer

Interactive Shiny app and reproducible R scripts for exploring gene-level associations in the public Sha/REMoDLB diffuse large B-cell lymphoma dataset.

The app lets users type a gene symbol and dynamically regenerate:

- Kaplan-Meier survival plots and risk tables
- expression boxplots by clinical or molecular metadata
- expression-expression correlation plots
- summary and per-sample tables

## Dataset

This repository uses the curated, gene-symbol-modernized Sha/REMoDLB files in:

```text
sha/curated_gene_symbols/
```

Current curated objects:

```text
clin_sha_tot.rds
clin_sha_rchop.rds
clin_sha_rbchop.rds
sha_total_expression.rds
sha_rchop_expression.rds
sha_rbchop_expression.rds
```

The total cohort contains 928 samples and 20,816 gene-symbol rows. Treatment subsets are R-CHOP and RB-CHOP.

Dataset provenance:

- GEO accession: GSE117556
- Platform: GPL14951, Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip
- Sha et al., Journal of Clinical Oncology, 2019
- DOI: 10.1200/JCO.18.01314
- PubMed ID: 30523719

## Repository Structure

```text
app/app.R                         Shiny app
scripts/analysis_functions.R      Shared app and batch-analysis functions
scripts/04_gene_survival_analysis.R
scripts/05_gene_boxplot_association_analysis.R
scripts/06_gene_expression_correlation_analysis.R
sha/curated_gene_symbols/         Curated public RDS files used by the app
reports/                          Dataset validation and preprocessing notes
results/                          Generated outputs, ignored by git
```

## Run Locally

From the repository root:

```r
shiny::runApp("app")
```

Or from a shell:

```bash
Rscript -e 'shiny::runApp("app", host = "127.0.0.1", port = 3842)'
```

## Batch Analyses

The Shiny app regenerates one selected analysis at a time. To rebuild the full output folders:

```bash
Rscript scripts/04_gene_survival_analysis.R
Rscript scripts/05_gene_boxplot_association_analysis.R
Rscript scripts/06_gene_expression_correlation_analysis.R
```

Generated files are written to `results/`, which is intentionally ignored by git.

## Deploy to shinyapps.io

1. Install packages:

```r
install.packages(c(
  "shiny", "DT", "dplyr", "ggplot2", "readr",
  "survival", "survminer", "tibble", "tidyr", "rsconnect"
))
```

2. Configure your shinyapps.io account in RStudio or with `rsconnect::setAccountInfo()`.

3. Deploy from the repository root:

```r
source("deploy_shinyapps.R")
```

The deployment uses `app/app.R` as the primary app file and includes `scripts/` plus `sha/curated_gene_symbols/`. The `.rscignore` file excludes generated outputs, reports, raw working RDS files, and supplementary files.

## Notes for Public Sharing

- The curated RDS files are public dataset derivatives and are required by the app.
- Generated result plots and tables are not stored in git because users can regenerate them.
- The app is exploratory and should be cited together with the original public dataset publication.

