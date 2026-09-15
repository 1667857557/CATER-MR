# eQTL input

## Cis eQTL

`eqtl_dir/<GENE>.txt.gz` contains the complete tested cis summary for one gene.

Accepted columns:

| Field | Accepted names |
|---|---|
| SNP | `SNP`, `rsid`, `variant_id`, `variant`, `ID` |
| chromosome | `CHR`, `chrom`, `chromosome` |
| position | `BP`, `POS`, `position`, `base_pair_location` |
| effect allele | `A1`, `EA`, `effect_allele`, `ALT` |
| other allele | `A2`, `NEA`, `other_allele`, `non_effect_allele`, `REF` |
| effect | `b`, `beta`, `BETA`, `effect`, `estimate` |
| standard error | `se`, `SE`, `stderr`, `standard_error` |
| P value | `p`, `P`, `pval`, `p_value`, `pvalue`; optional |
| effect-allele frequency | `freq`, `EAF`, `eaf`, `effect_allele_frequency`, `AF`; optional |
| sample size | `N`, `n`, `samplesize`, `sample_size`; optional when `qtl_n` is supplied |

## Significant-only trans eQTL

`trans_eqtl` is a data frame or tab-delimited file with one gene column plus the summary-statistic fields above.

Accepted gene names:

```text
GENE, gene, gene_name, gene_symbol, symbol, SYMBOL
```

Contract:

```text
trans_data_mode = "significant_only"
trans_eqtl != NULL
instrument_p <= trans_reporting_p
```

## Outcome

`outcome` uses the same SNP, allele, effect, standard-error and optional P/EAF/N aliases. Chromosome and position are optional.

## Full-summary table adapter

```r
source("CATER_MR.R")
source("CATER_EQTL_INPUT.R")

cache <- cater_prepare_eqtl_table(
  eqtl_table = "full_cis_eqtl.tsv.gz",
  gene_col = "GENE",
  outdir = "eqtl/cis"
)
```

`cater_prepare_eqtl_table()` materializes `<GENE>.txt.gz` files and does not perform significance filtering.

## Direct run from a merged full-summary table

```r
fit <- cater_mr_from_table(
  grn = grn,
  eqtl_table = "full_cis_eqtl.tsv.gz",
  outcome = outcome,
  gene_col = "GENE",
  trans_eqtl = trans_hits,
  gene_annotation = annotation,
  ld_bfile = "ld/EUR_GRCh38",
  instrument_p = 5e-8,
  trans_reporting_p = 5e-8
)
```

## Functions

| Function | Description |
|---|---|
| `cater_prepare_eqtl_table()` | Materializes a merged full-summary table to per-gene files |
| `cater_mr_from_table()` | Materializes required genes and calls `cater_mr()` |
| `print.cater_eqtl_cache()` | Prints cache metadata |
