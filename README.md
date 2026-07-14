# DLBCL Public Dataset Explorer

Interactive Shiny app and reproducible R curation scripts for exploring gene-level associations across public diffuse large B-cell lymphoma (DLBCL) expression datasets.

## Live App

The deployed Shiny app is available at:

https://aledavini.shinyapps.io/dlbcl-public-dataset-explorer/

The app lets users select a public dataset, cohort, gene symbol, and available metadata annotations to dynamically regenerate:

- Kaplan-Meier survival plots and risk tables
- expression boxplots by clinical or molecular metadata
- expression-expression correlation plots
- downloadable PNG/PDF plots with editable export size
- downloadable summary and per-sample CSV tables

## Datasets

The app currently supports four curated public DLBCL expression datasets.

| App dataset | Source | Samples used | Expression features | Outcomes | Main available annotations |
| --- | --- | ---: | ---: | --- | --- |
| Sha/REMoDLB | GSE117556, GPL14951 | 928 | 20,816 gene-symbol rows | OS, PFS | treatment arm, COO, MYC/BCL2/BCL6 rearrangements, hit status, double-expressor annotations, MYC/BCL2 RNA and IHC |
| Lenz GSE10846 | GSE10846, GPL570 | 414 | 21,300 gene symbols | OS | treatment, COO, stage, gender |
| Chapuy GSE98588 | GSE98588, GPL23432 | 137 | 12,636 gene symbols | OS, PFS | R-CHOP-like treatment, COO, CCC, Chapuy cluster, genomic complexity, MYC rearrangement, MYC/BCL2 RNA, TP53, 6q deletion |
| Barrans GSE32918 | GSE32918, GPL8432 | 172 | 18,402 gene symbols | OS | R-CHOP status, COO, gender |

### Dataset-Specific Notes

- **Sha/REMoDLB** is the original curated dataset for this project. Probe labels were verified against GPL14951 and modernized to gene symbols.
- **Lenz GSE10846** is OS-only in the inspected public/local metadata. Ambiguous multi-symbol rows containing `///` were excluded from direct app search.
- **Chapuy GSE98588** already had expression rows as clean gene symbols in the local object. Mutation/rearrangement supplements were inspected but are not yet part of the mutation/CNA app chapter.
- **Barrans GSE32918** was curated from Illumina probe-level expression. Technical replicate arrays were averaged to patient-level samples, probes were mapped with GPL8432, and duplicate probes per gene were resolved by highest MAD. Five negative raw follow-up times were set to missing for OS analyses.

## Curated Data Layout

Curated app-ready objects are stored locally in:

```text
sha/curated_gene_symbols/
curated_datasets/lenz/
curated_datasets/chapuy/
curated_datasets/barrans/
```

Large raw and curated data files are intentionally ignored by git, but they are included in the shinyapps.io deployment bundle by `deploy_shinyapps.R`.

## Reproducible Curation

Dataset curation scripts:

```text
scripts/07_curate_lenz_dataset.R
scripts/08_curate_chapuy_dataset.R
scripts/09_curate_barrans_dataset.R
```

Dataset inspection reports:

```text
reports/sha_dataset_validation.md
reports/sha_gene_symbol_probe_annotation_check.md
reports/lenz_dataset_inspection.md
reports/chapuy_dataset_inspection.md
reports/barrans_dataset_inspection.md
```

The shared app registry and analysis functions live in:

```text
scripts/analysis_functions.R
```

## Repository Structure

```text
app/app.R                         Shiny app
deploy_shinyapps.R                shinyapps.io deployment helper
scripts/analysis_functions.R      Shared app and batch-analysis functions
scripts/04_gene_survival_analysis.R
scripts/05_gene_boxplot_association_analysis.R
scripts/06_gene_expression_correlation_analysis.R
scripts/07_curate_lenz_dataset.R
scripts/08_curate_chapuy_dataset.R
scripts/09_curate_barrans_dataset.R
reports/                          Dataset validation and preprocessing notes
external_datasets_raw/            Local raw dataset staging area, ignored by git
curated_datasets/                 Local curated multi-dataset RDS files, ignored by git
sha/curated_gene_symbols/         Local curated Sha/REMoDLB RDS files
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

The Shiny app regenerates one selected analysis at a time. The original Sha-focused batch scripts can rebuild full output folders for selected genes:

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

The deployment uses `app/app.R` as the primary app file and explicitly includes the curated RDS files required at runtime.

## Notes for Public Sharing

- The app is exploratory and should be cited together with the original public dataset publications.
- Users should compare results across datasets with care because platforms, preprocessing, cohorts, and available metadata differ.
- Generated result plots and tables are not stored in git because users can regenerate them from the app.
- Future planned extensions may include additional expression datasets and a separate mutation/CNA visualization chapter.
