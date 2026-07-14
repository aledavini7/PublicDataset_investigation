# Public_dataset_investigation.R Review

## Scope

Reviewed `Public_dataset_investigation.R` as a reference only. The file was not modified.

## Main Findings

1. The script points to hardcoded HPC paths on lines 14-15. New scripts should use project-relative paths such as `sha/sha_total_expression.rds`.

2. The expression and clinical files contain the same sample IDs, but not in identical order. New scripts must reorder metadata with `meta <- meta[match(colnames(expr), meta$Sample_ID), ]` before modeling.

3. Lines 32-34, 125-127, and 224-226 assign `gene` three times. Only the last value, `CD82`, is used unless sections are manually rerun after editing.

4. Lines 42-53, 135-152, and 234-245 create CXCR4/CXCL12/CD82 groups using the same `pzhigh`/`pzlow` object from the currently selected gene. This can silently assign CD82 groups to all three gene labels.

5. Lines 56-65, 155-164, and 278-287 overwrite `fit_os` and `fit_pfs` repeatedly. Only the last model survives for plotting and p-value reporting.

6. Lines 119-120, 218-219, and 342-343 call `surv_pvalue(..., data = clin)`, but `clin` is not defined in the script.

7. Output paths on lines 108-112, 207-211, and 331-335 assume `figtabs/km_curves/update/` already exists. New scripts should create output directories before saving.

8. The x-axis label says "Time (years)", but the variables are named `Mesi_OS` and `Mesi_PFS`. Confirm whether these are months or already converted years before final plotting.

9. `DESeq2` is loaded on line 3, but the current expression matrix is normalized microarray-like data, not RNA-seq counts. For the current analyses, `limma`, `survival`, and `survminer` are more appropriate than `DESeq2`.

10. The MYC-combined analysis on lines 270-272 overwrites `meta_sha_ff` three times, so only the final CD82 filter is retained.

## Gene Symbol Check

Confirmed locally:

- 20,816 expression rows.
- No duplicated gene row names.
- No empty row names.
- No Illumina probe IDs remain as row names.
- Target genes `CXCR4`, `CXCL12`, and `CD82` are present.

Important nuance:

- 25 row names contain spaces, including `March 1`, `March 10`, `Septin 1`, and `Selenoprotein 15`.
- Offline `org.Hs.eg.db` recognized 16,813 / 20,816 row names as current SYMBOL keys and did not recognize 4,003. Many unrecognized names are likely legacy aliases rather than true errors.
- For the app, we should support alias/current-symbol lookup instead of assuming every user-entered gene name matches the curated row name exactly.

## Recommended New Script Structure

1. Load curated expression and clinical data.
2. Validate and reorder clinical metadata to expression column order.
3. Use functions for:
   - gene lookup
   - median high/low grouping
   - quartile high/low grouping
   - Kaplan-Meier fitting
   - plot generation
   - tabular summary export
4. Iterate over a vector of genes instead of manually rewriting blocks.
5. Save one results table and one plot set per gene/outcome/subset.

## First Reusable Analysis To Build

Build `scripts/02_gene_survival_analysis.R` for:

- input genes: `CXCR4`, `CXCL12`, `CD82`
- datasets: total, R-CHOP, RB-CHOP
- grouping: median and upper/lower quartiles
- outcomes: OS and PFS
- outputs: PDF plots and CSV summary tables

