# CATER-MR v0.8 methodology

CATER-MR v0.8 is designed for the asymmetric cell-type QTL data contract in which **full cis association summaries are available, whereas trans associations may be available only after genome-wide-significance ascertainment**.

## Core estimand and analysis roles

For a target gene `X`, CATER-MR reports four distinct analysis roles:

- `cis`: the primary target-specific MR anchor;
- `trans`: secondary evidence from GRN-constrained trans instruments;
- `combined`: cis + GRN-trans augmentation / triangulation;
- `network`: conditional MVMR sensitivity when measured competing sibling pathways and complete cross-exposure SNP effects are available.

The GRN is a biological eligibility gate. GRN edge weights do **not** rescale eQTL effects.

## Input contract

### Full cis QTL

`eqtl_dir/<GENE>.txt.gz` supplies the association rows available for each gene. Under the v0.8 strict contract, the cis region must be complete for every analyzed gene (`cis_full_summary=TRUE`).

### Significant-only trans QTL

`trans_eqtl` may be either a data frame or a tab-delimited file containing reported trans SNP-gene pairs, with a gene column plus standard SNP summary-statistic fields.

When the manifest declares `trans_data_mode="significant_only"`, an explicit `trans_eqtl` input is required. CATER-MR fails closed rather than silently falling back to legacy full-summary trans extraction from `eqtl_dir`.

The default analysis threshold is:

```text
instrument_p = 5e-8
```

This is a common **instrument-eligibility threshold**, not a claim that cis and trans discovery scans have identical family-wise multiple-testing properties.

`trans_reporting_p` describes the source-study availability threshold. For significant-only input, the requested eligibility threshold must satisfy:

```text
instrument_p <= trans_reporting_p
```

If `trans_reporting_p < instrument_p`, eligible trans associations between the two thresholds are unobserved by construction, so CATER-MR stops rather than analyzing an incomplete instrument pool. If `trans_reporting_p > instrument_p`, CATER-MR simply applies the more stringent `instrument_p` selection to the available catalog.

`CATER_EQTL_INPUT.R` remains a legacy full-summary storage adapter. It must not be used to make a significance-censored trans catalog appear to be a full gene-level summary table.

## Candidate construction

For target `X`, cis candidates satisfy:

```text
P(G -> X) < instrument_p
and G lies in the physical target-cis component.
```

Trans candidates satisfy:

```text
P(G -> X) < instrument_p
and G lies in a locus of a direct GRN parent TF T of X
and the TF locus is not connected to the target-cis physical component.
```

The existing connected-interval rule is retained: any connected component containing the target cis window is classified as cis in full, preventing one local LD neighborhood from being represented simultaneously as cis and trans evidence.

Association provenance is retained (`FULL_CIS`, `REPORTED_TRANS`, or legacy `FULL_SUMMARY_TRANS`).

## Genetic selection and LD

Manc-COJO is no longer the default IV selector. Eligible cis and trans candidates are processed using the same deterministic significance + LD-selection rule.

Candidates are ordered by target-association P value. Within the configured genomic window, a candidate is removed when its squared LD with an already retained candidate exceeds `ld_clump_r2` (default `0.01`). The retained signed LD matrix is still propagated into downstream generalized IVW, so pruning is a redundancy/numerical-stability operation rather than an assumption that residual LD is zero.

Default genetic QC continues to include allele compatibility, LD-reference availability, non-finite LD handling, palindromic-variant removal and hg38 MHC exclusion.

## MR estimation

CATER-MR continues to use marginal SNP-exposure and SNP-outcome effects. For correlated instruments:

```text
Omega_Y = D_Y R D_Y
beta = (gamma' Omega_Y^-1 Gamma) / (gamma' Omega_Y^-1 gamma)
```

where `R` is signed LD. COJO conditional effect estimates are not required for this estimator.

## TF-anchor evidence

For each selected GRN-trans IV `G` mapped to parent TF `T`, CATER-MR queries the cis-QTL evidence `G -> T` when available. This remains mechanistic evidence rather than a hard default exclusion gate because lack of a strong TF-expression cis-eQTL does not establish absence of regulation through TF activity, splicing or other mechanisms.

## Sibling and hotspot diagnostics

For a selected trans IV, CATER-MR screens direct siblings of the same parent TF while excluding target descendants when preserving the total target effect.

A missing trans SNP-gene pair in a significance-censored catalog is **not** interpreted as a zero effect. A negative sibling conclusion is permitted only when every requested cross-effect is numerically available and the omnibus test is valid. Otherwise the screen is marked incomplete.

For significant-only trans input, CATER-MR additionally reports how many distinct reported trans targets are associated with each selected trans IV as a hotspot diagnostic. Hotspot status is a pleiotropy flag, not an automatic invalid-IV rule.

## MVMR

MVMR is retained. Significance-based IV selection is not itself incompatible with MVMR.

For target `X` and active competing sibling exposures, CATER-MR forms the union of eligible instruments and then requires a complete matrix:

```text
selected SNPs x all MVMR exposures -> beta and SE
```

A SNP does not need to be significant for every exposure, but its association estimate with every included exposure must be available. Therefore significance-censored missing trans cross-effects cannot be set to zero and cannot be handled by silently dropping incomplete rows.

If any required cross-effect is unavailable, CATER-MR returns:

```text
MVMR_CROSS_EXPOSURE_EFFECTS_UNAVAILABLE
```

`cross_effect_lookup` may be supplied as a function or table to retrieve non-significant cross-exposure beta/SE values from the original QTL mapping or another complete lookup source.

The correlated-IV GLS MVMR point estimator is retained. The current custom correlated conditional-F calculation is labelled:

```text
EXPERIMENTAL_CORRELATED_IV_CONDITIONAL_F
```

Network MVMR cannot become the primary result by default. Promotion requires all of the following explicit gates in addition to the numerical eligibility checks:

```text
primary_policy = "screened_cater"
accept_experimental_conditional_f = TRUE
allow_network_primary = TRUE
sibling_screen_independent = TRUE
```

The last flag is an explicit user assertion that the sibling-screen evidence used to trigger adjustment is suitable for primary-model selection rather than merely exploratory reuse of the same ascertainment process. These gates remain `FALSE` by default. The custom conditional-strength implementation should be externally benchmarked before conventional conditional-F cutoffs are given strong interpretation.

## Primary-result policy

The default remains conservative:

```text
valid cis MR -> primary cis estimate
trans MR -> secondary network-constrained evidence
combined MR -> CATER augmentation/sensitivity
MVMR -> conditional sensitivity when identifiable
```

Hit-only absence is never sufficient to promote a combined result as if trans pleiotropy had been excluded.

## Important unresolved assumptions

The following are deliberately not treated as settled facts:

- `ld_clump_r2=0.01` is a conservative default, not a proven optimum; sensitivity analyses should assess alternative thresholds.
- A predicted cell-type GRN edge is not automatically a causal regulatory edge; perturbational/independent validation can change confidence in trans qualification.
- TF-anchor evidence is not yet validated as a hard inclusion/exclusion rule.
- A trans hotspot is not automatically horizontal pleiotropy; it can include vertical downstream network propagation.
- The custom correlated conditional-F implementation is experimental until externally benchmarked/calibrated.
- Significant-only trans beta estimates are subject to ascertainment/winner's-curse considerations.

These uncertainties should be tested by simulation and negative-control analyses rather than resolved by implementation convention.
