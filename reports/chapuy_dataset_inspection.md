# Chapuy GSE98588 Dataset Inspection

## Reference and Source

- Dataset: Chapuy GSE98588
- Reference: Chapuy et al., Cancer Cell 2018
- GEO accession: GSE98588
- GEO title: Genetically-defined Diffuse Large B-cell Lymphoma Subsets Arise by Distinct Pathogenetic Mechanisms and Predicts Outcome
- GEO platform: GPL23432, Affymetrix Human Genome U133 Plus 2.0 Array with Brainarray ENSG v18 CDF
- PubMed ID: 29713087
- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98588

GEO describes 137 expression samples, combining 52 newly generated DLBCL profiles with 85 previously generated Affymetrix U133 Plus 2.0 profiles from GSE34171.

## Raw Files Inspected

- `external_datasets_raw/chapuy/gep_expression_geo_GSE98588.rds`
- `external_datasets_raw/chapuy/clin_137_gep.rds`
- `external_datasets_raw/chapuy/clin_total.xlsx`
- `external_datasets_raw/chapuy/mutations_rearrangements.xlsx`

## Expression Integrity

- Raw expression dimensions: 12,636 genes x 137 samples
- Raw expression object: data frame, converted to numeric matrix in curation
- Expression row labels: gene symbols
- Duplicated gene labels: 0
- Ambiguous `///` gene labels: 0
- Duplicated sample labels: 0
- Missing expression values: 0
- Expression columns exactly match clinical `pair_id`

## Clinical Integrity

- Raw clinical dimensions: 137 samples x 73 fields
- Curated clinical dimensions: 137 samples x 47 fields
- Complete PFS records: 100
- Complete OS records: 105
- Event coding preserved as `1 = event`, `0 = censored`
- Survival times are kept in months

Key curated annotations:

- R-CHOP-like treatment: yes 101, no 29, missing 7
- COO by GEP: ABC 63, GCB 54, UNC 20
- Consensus clustering class: BCR 39, HR 40, OxPhos 43, missing 15
- Chapuy cluster: C0 6, C1 29, C2 32, C3 28, C4 21, C5 21
- Genomic complexity: clean 44, complex 87, missing 6
- MYC rearrangement: rearranged 11, not rearranged 126
- TP53 status: mutated 15, not-mutated 122
- 6q deletion: deleted 23, normal 114

## Curated Outputs

- `curated_datasets/chapuy/expression.rds`
- `curated_datasets/chapuy/clinical.rds`
- `curated_datasets/chapuy/dataset_info.rds`
- `curated_datasets/chapuy/curation_summary.csv`

The curated dataset is ready for app-level expression analyses. Mutation and rearrangement supplement files were inspected but are left for the future mutation/CNA chapter rather than being fully integrated in this expression-focused pass.
