# Barrans GSE32918 Dataset Inspection

## Reference and Source

- Dataset: Barrans GSE32918
- Reference label: Barrans et al. GSE32918
- GEO accession: GSE32918
- GEO title: Whole genome expression profiling based on paraffin embedded tissue can be used to classify diffuse large b-cell lymphoma and predict clinical outcome
- Platform: GPL8432, Illumina HumanRef-8 WG-DASL v3.0
- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32918
- PubMed IDs listed by GEO: 22970711, 24875472

GEO describes RNA from FFPE biopsies for 172 DLBCL patients profiled on Illumina DASL. The raw matrix contains 249 arrays because several patients have technical replicate arrays.

## Raw Files Inspected

- `external_datasets_raw/barrans/barrans_matrix.tsv`
- `external_datasets_raw/barrans/barrans_metadata.tsv`
- `external_datasets_raw/barrans/GPL8432-11703.tsv`

## Expression and Probe Annotation

- Raw probe rows: 24,526
- Raw arrays: 249
- Patient-level samples after averaging technical replicates: 172
- Platform annotation matched all matrix probe IDs
- Symbol-annotated probes: 24,526
- Curated gene symbols after duplicate-probe selection: 18,402
- Duplicated curated gene symbols: 0
- Missing curated expression values: 0

Technical replicate handling:

- Base sample IDs were derived by removing `_RepN` suffixes from GEO sample titles.
- Replicate arrays were averaged per patient/sample before probe-to-gene reduction.
- 39 patient-level samples had more than one technical array.

Probe-to-gene handling:

- Probes were mapped to gene symbols using `GPL8432-11703.tsv`.
- For genes represented by multiple probes, the probe with highest MAD across patient-level samples was selected.
- Selected probe annotations are stored in `curated_datasets/barrans/selected_probe_annotation.csv`.

## Clinical Integrity

- Curated clinical samples: 172
- R-CHOP treated: 140
- Not treated with curative intent: 32
- COO classes: ABC 53, GCB 82, UNC/TypeIII 37
- Follow-up status: Alive 79, Dead 93
- Complete usable OS records: 167
- OS events after filtering invalid times: 88

Survival handling:

- OS is available as follow-up years.
- No PFS/progression field was present in the inspected metadata.
- Five patients had negative raw follow-up years and were set to missing for OS analyses:
  - 54: -0.00548
  - 117: -0.03830
  - 133: -0.00548
  - 159: -0.01370
  - 160: -0.03560
- `OS_event` is coded as `1 = dead/event`, `0 = alive/censored`.

## Curated Outputs

- `curated_datasets/barrans/expression.rds`
- `curated_datasets/barrans/clinical.rds`
- `curated_datasets/barrans/dataset_info.rds`
- `curated_datasets/barrans/curation_summary.csv`
- `curated_datasets/barrans/selected_probe_annotation.csv`
- `curated_datasets/barrans/probe_annotation_selection.csv`

This dataset is ready for app-level expression analyses with OS only.
