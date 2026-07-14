# Sha Probe-To-Gene Label Verification

Generated: 2026-07-13 15:01:37 CEST

## Source

- Platform: GPL14951, Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip.
- Probe annotation fields used: `Symbol` and `ILMN_Gene`.

## Summary

                                                              check     n
            original row names present in GPL14951 Symbol/ILMN_Gene 20816
              modern row names present in GPL14951 Symbol/ILMN_Gene 20791
 modern row names absent from GPL14951 but intentionally modernized    25
                            original row names absent from GPL14951     0
         row names changed between original and modernized datasets    25
 denominator
       20816
       20816
       20816
       20816
       20816

## Interpretation

- All original curated expression row names are present in GPL14951 `Symbol` or `ILMN_Gene`.
- Therefore, the original curated gene labels are consistent with the probe annotation used by the platform.
- The modernized dataset differs only in the 25 intentional legacy-symbol updates.
- Those 25 modern symbols are not direct GPL14951 labels; they are current-symbol replacements for downstream usability.

## Intentional Modernizations

  original_symbol modern_symbol
          March 1       MARCHF1
         March 10      MARCHF10
         March 11      MARCHF11
          March 2       MARCHF2
          March 3       MARCHF3
          March 4       MARCHF4
          March 5       MARCHF5
          March 6       MARCHF6
          March 7       MARCHF7
          March 8       MARCHF8
          March 9       MARCHF9
 Selenoprotein 15       SELENOF
         Septin 1       SEPTIN1
        Septin 10      SEPTIN10
        Septin 11      SEPTIN11
        Septin 12      SEPTIN12
        Septin 13     SEPTIN7P2
        Septin 14      SEPTIN14
         Septin 2       SEPTIN2
         Septin 3       SEPTIN3
         Septin 4       SEPTIN4
         Septin 5       SEPTIN5
         Septin 6       SEPTIN6
         Septin 7       SEPTIN7
         Septin 9       SEPTIN9
