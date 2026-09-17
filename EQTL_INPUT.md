# eQTL input

## Required data structure

Lean CATER-MR is designed for the common cell-type QTL setting in which:

- **full cis-eQTL summary statistics** are available per gene;
- **trans-eQTL data are significance-censored**, typically genome-wide-significant SNP-gene pairs only;
- the GRN is a cell-type-matched directed TF→target network.

The full cis files are used twice: first to construct target cis instruments, and second to verify that each candidate trans SNP actually perturbs its proposed parent TF in cis.

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

For every parent TF that may contribute trans instruments, its per-gene file must contain the tested cis-region SNPs, including non-significant associations. Missing SNPs are treated as unavailable evidence, never as zero effects.

## Significant-only trans schema

Gene column:

```text
GENE, gene, gene_name, gene_symbol, symbol, SYMBOL
```

Association fields use the cis-summary aliases above.

Recommended strict contract:

```text
cis_full_summary = TRUE
trans_data_mode = "significant_only"
trans_eqtl != NULL
instrument_p <= trans_reporting_p
```

For significance-censored trans data, absence of a SNP-gene pair means only that it was not present in the reported catalog. It does **not** imply `beta = 0`.

## Lean trans qualification

For a target `X`, a reported trans SNP `G` first has to map to the cis region of a direct GRN parent `TF -> X`. CATER-MR then queries the full cis file for that TF and requires the same SNP to satisfy:

```text
p_tf < tf_anchor_p
```

with default:

```text
tf_anchor_p = instrument_p = 5e-8
```

The default core set additionally requires exactly one qualifying parent TF and:

```text
n_reported_trans_targets <= max_reported_trans_targets
max_reported_trans_targets = 1
```

This is an observed-specificity rule only. It does not demonstrate absence of sub-threshold trans effects or other horizontal-pleiotropic pathways.

Use:

```text
trans_set = "extended"
```

to retain all unique-parent, TF-cis-anchored trans instruments regardless of the number of other reported trans targets. This is intended for sensitivity analysis.

## Full-summary trans data

`trans_data_mode="full_summary"` remains supported for compatibility. The reported-target-count filter is primarily intended for significance-censored catalogs; with genuinely full trans summary data, users should interpret cross-gene effects directly rather than treating non-significance as absence of effect.

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
