# CATER-MR

**Cis And Trans eQTLs guided by Regulatory networks for drug-target Mendelian randomization.**

## Lean v0.9 design

CATER-MR uses a conservative cis MR estimate as the primary causal estimate and uses one-hop, cell-type-specific TF→target GRN edges to qualify mechanistically supported trans-eQTL instruments.

For a target gene `X`, a trans SNP is eligible for the default core set only when:

1. `G -> X` is a reported genome-wide-significant trans-eQTL;
2. `G` lies in the cis region of a direct GRN parent `TF -> X`;
3. the same SNP has a genome-wide-significant cis-eQTL effect on that parent TF in the full cis summary data;
4. exactly one direct parent TF passes the cis-anchor criterion;
5. by default, the SNP has no other reported genome-wide-significant trans target in the supplied cell-type trans catalog.

The last rule is an **observed specificity filter**, not proof of no pleiotropy. A non-reported trans association is censored by the source reporting threshold and must never be interpreted as a zero effect.

The core models are:

- `cis`: primary MR estimate;
- `trans`: qualified trans-IV triangulation estimate;
- `augmented`: signed-LD GIVW using cis plus qualified trans IVs.

For backward compatibility, the result table also retains the historical `combined` row as an alias of `augmented`.

Sibling co-perturbation screening, MVMR, and the experimental correlated conditional-F diagnostic remain available as advanced sensitivity analyses, but are disabled by default and cannot replace the cis primary estimate.

## Why the GRN is used as a gate rather than a weight

The RNA+ATAC GRN supplies a cell-type-specific one-hop regulatory prior (`TF -> target`). CATER-MR does not multiply MR effects by GRN coefficients and does not interpret GRN edge magnitude as an IV-validity probability. The GRN and TF cis-eQTL evidence reduce mechanism ambiguity, but they do **not** prove the MR exclusion restriction.

## Main functions

| Function | Description |
|---|---|
| `cater_mr()` | Lean cis-primary MR with GRN-mediated trans qualification, signed-LD trans/augmented MR, and optional advanced sensitivity analyses |
| `cater_prepare_eqtl_table()` | Materializes merged full cis-eQTL summaries to `<GENE>.txt.gz` files |
| `cater_mr_from_table()` | Runs `cater_mr()` from a merged full cis-eQTL table |
| `cater_build_scmore_grn()` | Builds cell-type GRNs with audited `scMORE::createRegulon()` |
| `cater_get_scmore_grn()` | Extracts a CATER-ready cell-type GRN or associated scMORE evidence |

## Default Lean controls

```r
enable_sibling_screen = FALSE
enable_mvmr = FALSE
instrument_p = 5e-8
tf_anchor_p = instrument_p
max_reported_trans_targets = 1L
trans_set = "core"
ld_clump_r2 = 0.01
```

`trans_set="extended"` retains unique-parent, TF-cis-anchored trans IVs even when the same SNP has other reported significant trans targets. This is intended as a sensitivity analysis.

## Outputs added in v0.9

Target summaries report:

- `n_trans_raw`
- `n_trans_anchor_pass`
- `n_trans_unique_parent`
- `n_trans_core`
- `n_trans_extended`
- `n_parent_tf_core`
- `cis_effective_F`
- `trans_effective_F`
- `augmented_effective_F`
- `primary_effective_F`
- `augmented_se_reduction`
- `augmented_precision_fraction`

The legacy `effective_F` field now maps to the actual primary model (`cis`) rather than automatically reporting the augmented model strength.

Per-target qualification tables are written under `mechanism/*_trans_qualification.tsv`.

## Interfaces

- [`EQTL_INPUT.md`](EQTL_INPUT.md)
- [`SCMORE_GRN.md`](SCMORE_GRN.md)
- [`MATHEMATICS.md`](MATHEMATICS.md)

## References

1. Aguet F, Brown AA, Castel SE, et al. Genetic effects on gene expression across human tissues. *Nature*. 2017;550:204-213. https://doi.org/10.1038/nature24277
2. Zheng J, Haberland V, Baird D, et al. Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics*. 2020;52:1122-1131. https://doi.org/10.1038/s41588-020-0682-6
3. Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC. Mendelian randomization with fine-mapped genetic data: Choosing from large numbers of correlated instrumental variables. *Genetic Epidemiology*. 2017;41:714-725. https://doi.org/10.1002/gepi.22077
4. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *International Journal of Epidemiology*. 2019;48:713-727. https://doi.org/10.1093/ije/dyy262
