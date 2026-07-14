# Lenz GSE10846 Dataset Inspection

Generated: 2026-07-14

## Scope

Raw folder inspected:

```text
external_datasets_raw/lenz/
```

Files:

```text
clin_lenz_chop.rds
clin_lenz_rchop.rds
lenz_chop_expression.rds
lenz_rchop_expression.rds
lenz_sha_tot.rds
lenz_total_expression.rds
```

Curated files were generated in:

```text
curated_datasets/lenz/
```

## Provenance

The local files correspond to GEO series **GSE10846**.

GEO record:

- Title: *Prediction of survival in diffuse large B cell lymphoma treated with chemotherapy plus Rituximab*
- Organism: Homo sapiens
- Experiment type: expression profiling by array
- Platform: **GPL570**, Affymetrix Human Genome U133 Plus 2.0 Array
- GEO design summary: 181 CHOP-treated and 233 R-CHOP-treated clinical samples
- GEO citations: PubMed IDs **19038878** and **21546504**
- GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846

The local expression and clinical objects match the GEO design sample split:

| cohort | expression dim | clinical dim |
| --- | ---: | ---: |
| total | 22,485 x 414 | 414 x 59 |
| CHOP | 22,485 x 181 | 181 x 59 |
| R-CHOP | 22,485 x 233 | 233 x 59 |

## Integrity Checks

Expression matrix:

- Rows: 22,485
- Columns/samples: 414
- Duplicate row names: 0
- Duplicate sample names: 0
- Any expression NA: false
- Expression range: 0 to 17.977
- Sample columns match clinical `geo_accession`: true
- Total sample order already identical to clinical `geo_accession`: true

Treatment split:

- CHOP expression/clinical objects match internally.
- R-CHOP expression/clinical objects match internally.
- CHOP and R-CHOP sample sets do not overlap.
- CHOP + R-CHOP sample sets cover the total expression object.

## Clinical Fields

Available app-relevant fields:

| field | type | missing | values/notes |
| --- | --- | ---: | --- |
| `geo_accession` | character | 0 | sample ID |
| `DLBCL_subtype` | character | 0 | ABC, GCB, Unclassified |
| `Status_OS` | numeric | 0 | 0/1; 1 corresponds to DEAD/event |
| `OS_y` | numeric | 0 | overall survival in years |
| `Treatment` | character | 0 | CHOP, R-CHOP |
| `Stage` | character | 8 | stages 1-4 plus missing text `NA` |
| `Gender` | character | 18 | male/female plus missing text `NA` |
| `Age` | character | 0 | numeric after conversion |
| `Clinical info:ch1` | character | 0 | semicolon-packed GEO clinical string |

No PFS field was identified in the current RDS objects. Local field-name searches found no PFS/progression/event-free survival variables, and GEO sample metadata exposes follow-up status and follow-up years only. Therefore, this dataset should be treated as **OS-only** unless an additional non-GEO clinical supplement is provided.

## Cohort Counts

Subtype by treatment:

| treatment | ABC | GCB | Unclassified | total |
| --- | ---: | ---: | ---: | ---: |
| CHOP | 74 | 76 | 31 | 181 |
| R-CHOP | 93 | 107 | 33 | 233 |
| total | 167 | 183 | 64 | 414 |

OS event status by treatment:

| treatment | alive/censored (`0`) | dead/event (`1`) | total |
| --- | ---: | ---: | ---: |
| CHOP | 76 | 105 | 181 |
| R-CHOP | 173 | 60 | 233 |
| total | 249 | 165 | 414 |

Stage by treatment:

| treatment | 1 | 2 | 3 | 4 | missing | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CHOP | 28 | 55 | 39 | 58 | 1 | 181 |
| R-CHOP | 38 | 67 | 58 | 63 | 7 | 233 |
| total | 66 | 122 | 97 | 121 | 8 | 414 |

OS follow-up:

| treatment | min years | median years | max years |
| --- | ---: | ---: | ---: |
| CHOP | 0.00 | 2.68 | 21.78 |
| R-CHOP | 0.00 | 2.12 | 10.29 |
| total | 0.00 | 2.38 | 21.78 |

## Gene Symbol Status

The expression object is already row-labeled by gene-like symbols rather than probe IDs.

Project genes checked and present:

```text
MYC, BCL2, BCL6, CD82, CXCR4, CXCL12, TP53, NFKBIA, STAT3
```

However:

- 1,185 row names contain multi-symbol annotations using `///`.
- Examples: `ABCB1 /// ABCB4`, `AKR1C1 /// AKR1C2`.
- If reduced to the first symbol before `///`, 562 duplicated symbols would be introduced.

Recommendation for app curation:

1. Keep only unambiguous single-symbol rows for direct gene search.
2. Save excluded multi-symbol rows to an annotation/exclusion table.
3. Do not silently choose the first symbol from ambiguous rows.
4. Preserve the original expression object path/metadata in the curation report.

## App Availability Assessment

Available for Lenz:

- Dataset/cohort selection: total, CHOP, R-CHOP
- Gene expression lookup: yes, after single-symbol filtering
- Survival: OS only
- Survival stratifications:
  - basic gene high/low
  - treatment, if using total cohort
  - DLBCL subtype / COO-like grouping
  - stage
  - gender
  - possible age group, if we define an age cutoff
- Boxplots:
  - gene expression by treatment
  - gene expression by DLBCL subtype
  - gene expression by stage
  - gene expression by gender
- Correlations:
  - gene vs gene expression correlations

Not available in the current Lenz files:

- PFS
- MYC/BCL2/BCL6 rearrangement
- double-hit/triple-hit status
- MYC/BCL2 RNA high/average metadata fields equivalent to Sha/REMoDLB
- IHC double-expressor metadata
- mutation/CNA data

## Next Recommended Step

The curation script is:

```text
scripts/07_curate_lenz_dataset.R
```

It writes:

```text
curated_datasets/lenz/expression.rds
curated_datasets/lenz/clinical.rds
curated_datasets/lenz/dataset_info.rds
curated_datasets/lenz/excluded_ambiguous_gene_rows.csv
curated_datasets/lenz/curation_summary.csv
```

Proposed standardized clinical fields:

| standardized field | source field | transformation |
| --- | --- | --- |
| `Sample_ID` | `geo_accession` | as character |
| `OS_months` | `OS_y` | `OS_y * 12` |
| `OS_event` | `Status_OS` | keep numeric; 1 = dead/event |
| `Treatment` | `Treatment` | CHOP / R-CHOP |
| `COO_class` | `DLBCL_subtype` | ABC / GCB / UNC |
| `Stage` | `Stage` | convert text `NA` to missing |
| `Gender` | `Gender` | convert text `NA` to missing |
| `Age` | `Age` | numeric |

## Curation Result

| output | result |
| --- | ---: |
| raw expression rows | 22,485 |
| curated expression rows | 21,300 |
| excluded ambiguous `///` rows | 1,185 |
| samples | 414 |
| CHOP samples | 181 |
| R-CHOP samples | 233 |
| OS events | 165 |
| OS censored | 249 |

Post-curation validation:

- expression/clinical sample alignment: true
- duplicated curated gene symbols: 0
- remaining `///` ambiguous rows: 0
- available outcome: OS
- unavailable outcome: PFS
