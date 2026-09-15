# eQTL input

## Functions

| Function | Description |
|---|---|
| `cater_prepare_eqtl_table()` | Materializes merged full cis-eQTL summaries to per-gene files |
| `cater_mr_from_table()` | Materializes required genes and calls `cater_mr()` |
| `print.cater_eqtl_cache()` | Prints eQTL-cache metadata |

## Cis summary schema

| Field | Accepted names |
|---|---|
| SNP | `SNP`, `rsid`, `variant_id`, `variant`, `ID` |
| chromosome | `CHR`, `chrom`, `chromosome` |
| position | `BP`, `POS`, `position`, `base_pair_location` |
| effect allele | `A1`, `EA`, `effect_allele`, `ALT` |
| other allele | `A2`, `NEA`, `other_allele`, `non_effect_allele`, `REF` |
| effect | `b`, `beta`, `BETA`, `effect`, `estimate` |
| standard error | `se`, `SE`, `stderr`, `standard_error` |
| P value | `p`, `P`, `pval`, `p_value`, `pvalue` |
| effect-allele frequency | `freq`, `EAF`, `eaf`, `effect_allele_frequency`, `AF` |
| sample size | `N`, `n`, `samplesize`, `sample_size` |

Per-gene files:

```text
eqtl_dir/<GENE>.txt.gz
```

## Significant-only trans schema

Gene column:

```text
GENE, gene, gene_name, gene_symbol, symbol, SYMBOL
```

Association fields use the cis-summary aliases above.

Contract:

```text
trans_data_mode = "significant_only"
trans_eqtl != NULL
instrument_p <= trans_reporting_p
```

## Outcome schema

```text
SNP + effect allele + other allele + effect + standard error
```

Optional:

```text
P value
EAF
N
chromosome
position
```
