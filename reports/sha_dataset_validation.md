# Sha/REMoDLB Curated Dataset Validation

Generated: 2026-07-13 14:46:17 CEST

## Provenance Anchor

- Local files match the Sha et al. 2019 JCO / REMoDLB public dataset structure.
- GEO accession: GSE117556.
- PubMed ID: 30523719.
- DOI: 10.1200/JCO.18.01314.
- Platform: GPL14951, Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip.

## Object Summaries

file | class | rows | columns | duplicated_rows | duplicated_columns | any_na
--- | --- | --- | --- | --- | --- | ---
sha/sha_total_expression.rds | data.frame | 20816 | 928 | 0 | 0 | FALSE
sha/sha_rchop_expression.rds | data.frame | 20816 | 469 | 0 | 0 | FALSE
sha/sha_rbchop_expression.rds | data.frame | 20816 | 459 | 0 | 0 | FALSE
sha/clin_sha_tot.rds | data.frame | 928 | 55 | 0 | 0 | TRUE
sha/clin_sha_rchop.rds | data.frame | 469 | 55 | 0 | 0 | TRUE
sha/clin_sha_rbchop.rds | data.frame | 459 | 55 | 0 | 0 | TRUE

## Internal Consistency Checks

check | result
--- | ---
total expression columns match total clinical Sample_ID | TRUE
R-CHOP expression columns match R-CHOP clinical Sample_ID | TRUE
RB-CHOP expression columns match RB-CHOP clinical Sample_ID | TRUE
total clinical Sample_ID has no duplicates | TRUE
total expression column names have no duplicates | TRUE
total expression row names have no duplicates | TRUE
total expression and clinical are already in identical order | FALSE

Note: expression and clinical samples match, but the total files are not in identical order. Analysis scripts should reorder clinical metadata to expression columns before modeling.

## Expression Summary

- Total expression range: 7.0976 to 16.1067
- Total expression contains NA: FALSE

## Treatment Counts

Treatment | n
--- | ---
R-CHOP | 469
RB-CHOP | 459

## Molecular Class Counts

Class | n
--- | ---
ABC | 249
GCB | 468
MHG | 83
UNC | 128

## Survival Field Missingness

field | missing_n
--- | ---
Mesi_OS | 0
Evento_OS | 0
Mesi_PFS | 0
Evento_PFS | 0

## Clinical Missingness By Field

field | missing_n | missing_pct
--- | --- | ---
TP53 | 528 | 56.9
EOT_CT_RESP_DT | 71 | 7.65
RESP_ASSESS | 71 | 7.65
STAGE | 1 | 0.11
AGE | 0 | 0
BCL2 | 0 | 0
BCL2_IHC | 0 | 0
BCL2_RNA | 0 | 0
BCL2_class | 0 | 0
BCL2_rearrangement | 0 | 0
BCL6_rearrangement | 0 | 0
COO_class | 0 | 0
Class | 0 | 0
DEATH_IND | 0 | 0
DEXP_mRNA | 0 | 0
DOB_DT | 0 | 0
ECOG | 0 | 0
EXTRA_BON_MAR_INV_IND | 0 | 0
Evento_OS | 0 | 0
Evento_PFS | 0 | 0
FU_DT | 0 | 0
FU_TIME | 0 | 0
GENDER | 0 | 0
GEO_sample_name | 0 | 0
Hit_rearrangement | 0 | 0
INV_EXTRANODAL_BAS_IND | 0 | 0
IPI | 0 | 0
IPI_SCORE | 0 | 0
LDH | 0 | 0
MAX_TUMOUR | 0 | 0
MBN_Signature | 0 | 0
MBN_Signature_class | 0 | 0
MYC | 0 | 0
MYC_IHC | 0 | 0
MYC_RNA | 0 | 0
MYC_rearrangement | 0 | 0
Mesi_OS | 0 | 0
Mesi_PFS | 0 | 0
NFKBIA | 0 | 0
NFKBIA_class | 0 | 0
PROG_DT | 0 | 0
PROG_IND | 0 | 0
PROG_TIME | 0 | 0
REG_DT | 0 | 0
STAT3 | 0 | 0
STAT3_class | 0 | 0
Sample_ID | 0 | 0
TRT_ARM | 0 | 0
Treatment | 0 | 0
expressor_IHC | 0 | 0
expressor_RNA | 0 | 0
retrospective_COO_class | 0 | 0
retrospective_class | 0 | 0
sample | 0 | 0
trial_class | 0 | 0

## Supplementary Files

- DS_JCO.18.01314-1.pdf
- DS_JCO.18.01314-2.xlsx
- DS_JCO.18.01314-3.xlsx
- DS_JCO.18.01314-4.xlsx
- DS_JCO.18.01314-5.xlsx
- DS_JCO.18.01314-6.xlsx
- DS_JCO.18.01314-7.xlsx
- ~$DS_JCO.18.01314-4.xlsx

## Notes

- Temporary Excel lock file detected: ~$DS_JCO.18.01314-4.xlsx
- This file should be ignored by analysis code.
