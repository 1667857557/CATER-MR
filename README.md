# CATER-MR

**Cis And Trans eQTLs guided by Regulatory networks for drug-target Mendelian randomization.**

Current workflow:

```text
full cis eQTL + reported trans eQTL + cell-type GRN
    -> common IV eligibility threshold
    -> target-cis / direct-parent-TF locus assignment
    -> allele/LD QC + deterministic LD pruning
    -> signed-LD cis / trans / combined MR
    -> TF-anchor + sibling + hotspot + locus diagnostics
    -> MVMR sensitivity when the complete SNP x exposure matrix is available
```

## Files

| File | Function |
|---|---|
| `CATER_MR.R` | Main entrypoint |
| `R/CATER_MR_core.R` | CATER-MR implementation |
| `R/LD_helpers.R` | PLINK signed-LD extraction |
| `CATER_EQTL_INPUT.R` | Full-summary eQTL table adapter |
| `SCMORE_GRN.R` | scMORE-to-CATER GRN adapter |
| `MATHEMATICS.md` | Mathematical specification |
| `EQTL_INPUT.md` | eQTL input schema |
| `SCMORE_GRN.md` | scMORE GRN adapter interface |

## Requirements

- R
- PLINK with `--r square` support
- ancestry-matched PLINK `bed/bim/fam` LD reference
- optional: scMORE, Pando, Seurat, Signac for `SCMORE_GRN.R`

## Inputs

| Argument | Input |
|---|---|
| `grn` | Data frame containing `TF`, `Target`; coordinates may also be supplied in GRN columns |
| `eqtl_dir` | Directory of complete cis summaries: `<GENE>.txt.gz` |
| `trans_eqtl` | Data frame or tab-delimited file of reported trans SNP-gene associations |
| `outcome` | Outcome GWAS summary-statistics data frame |
| `gene_annotation` | `symbol`, `chr`, `tss`; optional when GRN contains coordinates |
| `ld_bfile` | PLINK reference prefix |
| `cross_effect_lookup` | Optional complete SNP-exposure lookup for MVMR |
| `input_manifest` | Optional data-contract metadata; required when `strict_input_contract=TRUE` |

Summary-statistic aliases accepted by the parser are listed in [`EQTL_INPUT.md`](EQTL_INPUT.md).

For significance-censored trans input:

```text
trans_data_mode = "significant_only"
trans_eqtl must be supplied
instrument_p <= trans_reporting_p
```

## Public functions

| Function | Description |
|---|---|
| `cater_mr()` | Runs CATER-MR for selected targets |
| `cater_prepare_eqtl_table()` | Materializes a merged full-summary eQTL table to gene files |
| `cater_mr_from_table()` | Runs `cater_mr()` from a merged full-summary eQTL table |
| `cater_build_scmore_grn()` | Builds cell-type GRNs with audited scMORE `createRegulon()` |
| `cater_get_scmore_grn()` | Extracts one CATER-ready cell-type GRN |

## Minimal run

```r
source("CATER_MR.R")

grn <- read.delim("grn.tsv")
trans_hits <- read.delim("trans_hits.tsv.gz")
outcome <- read.delim("outcome.tsv.gz")
annotation <- read.delim("gene_annotation.tsv")

fit <- cater_mr(
  grn = grn,
  eqtl_dir = "eqtl/cis",
  trans_eqtl = trans_hits,
  outcome = outcome,
  gene_annotation = annotation,
  ld_bfile = "ld/EUR_GRCh38",
  instrument_p = 5e-8,
  trans_reporting_p = 5e-8,
  ld_clump_r2 = 0.01,
  primary_policy = "cis_anchor",
  input_manifest = list(
    grn_build = "GRCh38",
    eqtl_build = "GRCh38",
    ld_build = "GRCh38",
    eqtl_n_unit = "donors",
    eqtl_ancestry = "EUR",
    ld_ancestry = "EUR",
    cis_full_summary = TRUE,
    trans_data_mode = "significant_only"
  ),
  strict_input_contract = TRUE
)
```

## Results

```text
CATER_MR_results/
├── cater_mr_results.tsv
├── cater_mr_target_summary.tsv
├── instruments/
├── diagnostics/
├── mechanism/
├── pleiotropy/
└── mvmr/
```

Models:

| Model | Role |
|---|---|
| `cis` | Primary target-specific anchor |
| `trans` | GRN-constrained secondary evidence |
| `combined` | cis + GRN-trans augmentation |
| `network` | Conditional MVMR sensitivity |

## Mathematics

See [`MATHEMATICS.md`](MATHEMATICS.md).

## References

1. Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC. *Genet Epidemiol.* 2017;41:714-725. https://doi.org/10.1002/gepi.22077
2. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. *Int J Epidemiol.* 2019;48:713-727. https://doi.org/10.1093/ije/dyy262
3. Fleck JS, Jansen SMJ, Wollny D, et al. *Nature.* 2023;621:365-372. https://doi.org/10.1038/s41586-022-05279-8
4. Ma Y, Yao Y, Zhou Y, et al. *Nature Aging.* 2026;6:270-289. https://doi.org/10.1038/s43587-025-01027-5
