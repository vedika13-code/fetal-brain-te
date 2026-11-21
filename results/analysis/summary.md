# TE expression in prenatal brain: summary

n = 1 per stage; descriptive only.

| sample | gene read pairs | TE (unique only) | TE (with multimappers) | gain from multimappers |
|---|---|---|---|---|
| early_CS20 | 2,766,673 | 286,828 | 298,703 | +4% |
| mid_14pcw | 2,713,275 | 278,944 | 290,372 | +4% |

## TE reads by class (% of TE reads, multimappers included)

| class      |   early_CS20 |   mid_14pcw |
|:-----------|-------------:|------------:|
| DNA        |         12.1 |        12.1 |
| LINE       |         31.7 |        31.1 |
| LTR        |         10.8 |        10.9 |
| Retroposon |          0.3 |         0.3 |
| SINE       |         45.1 |        45.6 |

## Young subfamilies: how much unique-only counting keeps

| subfamily   | class      |   mean_div |   uniq_total |   multi_total |   pct_recovered_by_unique |
|:------------|:-----------|-----------:|-------------:|--------------:|--------------------------:|
| L1HS        | LINE       |       19.4 |        209.5 |         529.6 |                      39.6 |
| L1PA2       | LINE       |       23.8 |        666.8 |        1103.6 |                      60.4 |
| AluYa5      | SINE       |       14.6 |        311.2 |         462.4 |                      67.3 |
| AluYb8      | SINE       |       23.5 |        242.3 |         314.5 |                      77   |
| SVA_F       | Retroposon |      118.1 |        173.7 |         256.7 |                      67.7 |

## Expressed TE loci (>=5 reads) by genomic context

| context            |   early_CS20 |   mid_14pcw |
|:-------------------|-------------:|------------:|
| exonic_sense       |         6307 |        5893 |
| exonic_antisense   |          149 |         170 |
| intronic_sense     |         2118 |        2549 |
| intronic_antisense |          343 |         325 |
| intergenic         |          659 |         626 |

## Top 15 subfamilies (CPM, relative to gene-assigned read pairs)

| subfamily   | class   |   cpm_early_CS20 |   cpm_mid_14pcw |
|:------------|:--------|-----------------:|----------------:|
| MIRb        | SINE    |           4645.1 |          4426.8 |
| AluJb       | SINE    |           3841.3 |          3819.4 |
| AluSx       | SINE    |           3724.1 |          3617.1 |
| AluSx1      | SINE    |           3303.2 |          3395.3 |
| MIR         | SINE    |           3222.9 |          3315.4 |
| AluSz       | SINE    |           3068.8 |          3183.6 |
| L2a         | LINE    |           2969.1 |          3092.4 |
| L2c         | LINE    |           2828.1 |          2631.4 |
| MIR3        | SINE    |           2755.9 |          2635.2 |
| AluY        | SINE    |           2574.4 |          2624.7 |
| AluJr       | SINE    |           2566.9 |          2538.1 |
| MIRc        | SINE    |           2535.8 |          2770.5 |
| AluJo       | SINE    |           2033   |          2091.7 |
| AluSp       | SINE    |           1973.7 |          1929.9 |
| AluSq2      | SINE    |           1912.5 |          1867.2 |

## SFARI autism genes

- Expressed TE loci (>=5 reads in a sample) inside SFARI genes: 1125, in 504 genes
- By context: exonic_sense=564, intronic_sense=522, exonic_antisense=20, intronic_antisense=19
- Top genes: PTPRD, RBFOX1, EIF3F, CELF2, PTPRS, ELAVL3, PCM1, DPYD, CUX2, FBXO11
