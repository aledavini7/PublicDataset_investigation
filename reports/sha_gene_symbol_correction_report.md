# Sha Gene Symbol Correction Report

Generated: 2026-07-13 14:53:05 CEST

## Output

- Corrected RDS folder: sha/curated_gene_symbols
- Original RDS files were not modified.
- Clinical RDS files were copied unchanged into the corrected folder for convenience.

## Correction Summary

                      file  rows columns duplicated_rows space_containing_rows
  sha_total_expression.rds 20816     928               0                     0
  sha_rchop_expression.rds 20816     469               0                     0
 sha_rbchop_expression.rds 20816     459               0                     0

## Mapping

       old_symbol new_symbol
          March 1    MARCHF1
          March 2    MARCHF2
          March 3    MARCHF3
          March 4    MARCHF4
          March 5    MARCHF5
          March 6    MARCHF6
          March 7    MARCHF7
          March 8    MARCHF8
          March 9    MARCHF9
         March 10   MARCHF10
         March 11   MARCHF11
 Selenoprotein 15    SELENOF
         Septin 1    SEPTIN1
         Septin 2    SEPTIN2
         Septin 3    SEPTIN3
         Septin 4    SEPTIN4
         Septin 5    SEPTIN5
         Septin 6    SEPTIN6
         Septin 7    SEPTIN7
         Septin 9    SEPTIN9
        Septin 10   SEPTIN10
        Septin 11   SEPTIN11
        Septin 12   SEPTIN12
        Septin 13  SEPTIN7P2
        Septin 14   SEPTIN14
                                                                 reason
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
      Legacy MARCH family name converted to current MARCHF-style symbol
                  Legacy selenoprotein name converted to current symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
            Legacy septin name converted to current SEPTIN-style symbol
 Legacy SEPT13 alias maps to SEPTIN7P2 in local org.Hs.eg.db annotation
            Legacy septin name converted to current SEPTIN-style symbol

## Notes

- These corrections target the 25 obvious space-containing legacy names found in the expression row names.
- `Septin 13` is handled conservatively using the local `org.Hs.eg.db` alias mapping from `SEPT13` to `SEPTIN7P2`.
- Broader alias modernization across all 20,816 rows should be treated as a separate step because many old symbols are biologically valid historical aliases and may require careful one-to-many handling.
