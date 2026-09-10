# CATER-MR

**Cis And Trans eQTLs guided by Regulatory networks for drug-target Mendelian Randomization**

CATER-MR is a target-centric framework that augments conventional cis molecular/drug-target MR with biologically constrained trans-eQTL signals from **one-hop upstream TF loci** in a cell-type-specific GRN.

## V0.5 design

The implementation follows the current Occam-style design:

```text
cell-type-specific GRN
    -> target cis + direct parent-TF loci
    -> one target-level Manc-COJO selection
    -> independent cis + GRN-trans signals
    -> allele-aligned IVW/GIVW
    -> TF exact-IV mechanism annotation
    -> complete same-IV one-hop sibling co-perturbation screen
    -> local MVMR sensitivity analysis only when measured sibling co-perturbation is detected
    -> BH-FDR across target genes
```

The core principle is:

> **SNPs are the instruments; the GRN is a biological gate; Manc-COJO is an independent-signal selector.**

GRN edge weights never rescale eQTL effects.

## Primary inputs

### Cell-type-specific GRN

`grn` is a `data.frame` with at least:

```text
TF  Target
```

or a list containing such a table in `$grn`.

Only direct `TF -> Target` edges are used. No two-hop or recursive expansion is performed.

### Cell-type-specific full-summary eQTL directory

Each measured gene is stored as:

```text
<eqtl_dir>/<SYMBOL>.txt.gz
```

Required/recognized information:

```text
SNP / rsid / variant_id
CHR
BP / POS
A1 / effect_allele
A2 / other_allele
EAF / freq
BETA / b
SE
P        # optional; derived from beta/se if absent
N        # or provide qtl_n
```

`A1` is the effect allele and `EAF` is the frequency of A1.

### Outcome

`outcome` is an in-memory summary-statistics `data.frame` with at least:

```text
SNP A1 A2 BETA SE
```

### Gene coordinates

Provide either:

```text
symbol chr tss
```

through `gene_annotation`, or include:

```text
TF_chr TF_tss Target_chr Target_tss
```

in the GRN.

### LD reference

`ld_bfile` is a PLINK `bed/bim/fam` prefix. The preferred reference is the genotype sample used for the eQTL mapping; otherwise use a large ancestry-matched panel.

## Manc-COJO

CATER-MR uses **Manc-COJO** (`light156/multi-ancestry-COJO`) in single-cohort mode for the current one-summary-file-per-gene input contract.

For each target `X`, CATER-MR performs **one** conditional-selection problem:

\[
\mathcal R_X=L_X\cup\bigcup_{T\in P_1(X)}L_T
\]

and runs:

```text
full X eQTL summary + --extract(candidate SNPs in R_X) + --cojo-slct
```

The target cis and all direct parent-TF loci are therefore conditioned jointly. This avoids treating LD-correlated cis and nearby TF-locus signals as independent merely because they were analyzed in separate windows.

When multiple SNPs are selected, CATER-MR runs a second joint-only Manc-COJO call on those selected SNPs to obtain the compact `.ldr.cojo` signed LD matrix.

## Allele orientation

Manc-COJO LD is defined using the reference-panel allele coding. CATER-MR therefore:

1. orients the target eQTL beta to the Manc-COJO A1;
2. flips beta/EAF when the QTL alleles are reversed;
3. harmonizes the outcome to the same A1/A2 system;
4. drops palindromic SNPs by default.

This is required for mathematically valid use of signed LD.

## MR estimators

CATER-MR reports:

```text
cis
trans
combined
network   # only when local MVMR is triggered
```

COJO is used for signal selection. The MR layer uses the original marginal SNP-exposure and SNP-outcome associations.

For correlated selected instruments:

\[
\Omega_Y=D_Y R D_Y
\]

and generalized IVW is:

\[
\hat\theta=
\frac{\gamma^T\Omega_Y^{-1}\Gamma}
     {\gamma^T\Omega_Y^{-1}\gamma}.
\]

No hidden ridge regularization is used: singular LD/information matrices are returned as explicit failure states.

## Instrument strength

Single-SNP strength is:

\[
F_j=(\beta_{jX}/SE_{jX})^2.
\]

For a correlated IV set CATER-MR also reports:

\[
I_X=\gamma^T(D_X R D_X)^{-1}\gamma
\]

and:

\[
F_{\mathrm{eff}}=I_X/K.
\]

When \(R=I\), this reduces to mean single-SNP F.

### Trans information fraction

Rather than counting trans SNPs, CATER-MR measures their incremental information:

\[
f_{\mathrm{trans}}
=
\frac{I_{\mathrm{all}}-I_{\mathrm{cis}}}{I_{\mathrm{all}}}.
\]

## TF-anchor mechanism annotation

For each selected GRN-trans IV `G` assigned to parent TF `T`, CATER-MR queries the exact SNP in `T.txt.gz`:

\[
\beta_{GT},SE_{GT}.
\]

This is a mechanism/evidence annotation, not a hard instrument gate. Target-level BH-FDR is applied across tested TF-anchor pairs.

## Same-IV one-hop sibling co-perturbation screen

For each actual selected trans IV \(G\) from a parent-TF locus, CATER-MR asks whether the **same IV** affects other direct targets of that TF.

For each valid sibling \(Z\):

\[
T\to X,\qquad T\to Z,
\]

while excluding direct descendants \(X\to Z\), CATER-MR queries:

\[
\beta_{GZ},SE_{GZ}
\]

directly from `Z.txt.gz`.

For multiple relevant trans IVs:

\[
Q_{TZ}
=
\beta_{GZ}^T
(D_Z R D_Z)^{-1}
\beta_{GZ}.
\]

For one SNP this is exactly the squared Wald z-statistic. BH-FDR is applied across the one-hop sibling candidate set for each target.

The screen is **coverage-aware**. For every sibling CATER-MR records the number of requested actual trans IVs, the number successfully queried and allele-aligned, and the number missing. A sibling is considered completely screened only when every requested IV is available and the omnibus test is numerically valid. A significant partial test can still identify measured sibling co-perturbation and trigger MVMR, but a non-significant test is only "no detected co-perturbation"; it is not proof of the MR exclusion restriction. Partial or failed coverage can never be interpreted as evidence that the remaining pathway is clean.

This direct same-IV test is intentionally used instead of requiring the sibling gene to select the identical COJO sentinel: an IV can affect a sibling through LD even when the sibling's own conditional analysis selects a different lead SNP.

## Triggered local MVMR

MVMR is triggered only if the exposure-side sibling screen detects at least one measured co-perturbed sibling. It remains a standard LD-aware GLS-MVMR sensitivity analysis and is not interpreted as a complete correction for unmeasured pleiotropy, weak-instrument bias, or sample overlap.

For active siblings:

\[
\mathcal I_X^{MVMR}
=
\mathcal S_X^{cis}
\cup
\mathcal S_X^{trans}
\cup
\bigcup_Z\mathcal S_Z^{cis}.
\]

Sibling cis signals are obtained with Manc-COJO only **after** the sibling is triggered.

The local model is:

\[
\Gamma_Y=B\theta+\epsilon,
\]

where the exposure matrix contains `X` and active siblings only. Parent TF is not included by default because its SNP-association vector can be nearly collinear with that of `X`.

### Conditional-F safeguard

Because gene-eQTL summaries within a cell type can come from the same donors, the errors in SNP-exposure estimates are correlated across exposures. With residual LD between instruments, the covariance also extends across SNPs. CATER-MR therefore uses the separable working covariance

\[
\operatorname{Cov}(\hat\gamma_k,\hat\gamma_l)
\approx
\rho_{kl}D_k R D_l,
\]

where \(D_k=\operatorname{diag}(SE_{1k},\ldots,SE_{mk})\). `exposure_corr` is the named matrix containing \(\rho_{kl}\), and is validated for names, symmetry, unit diagonal and positive semidefiniteness.

For exposure \(i\), conditional strength is evaluated from

\[
r_i(\delta)=\hat\gamma_i-\hat\Gamma_{-i}\delta,
\]

with residual covariance

\[
\Omega_i(\delta)
=
\sum_{k,l}q_kq_l\,\rho_{kl}D_kRD_l,
\qquad
q_i=1,\;q_{-i}=-\delta.
\]

The nuisance coefficients are obtained by feasible covariance-weighted GLS under \(\Omega_i(\delta)\), not by an unweighted regression. The reported statistic is

\[
\boxed{
F_{cond,i}
=
\frac{r_i^T\Omega_i^{-1}r_i}{m-p+1}
}
\]

for \(m\) SNPs and \(p\) exposures. The \(m-p+1\) denominator is the residual degrees of freedom after fitting the \(p-1\) nuisance exposure-association vectors.

For the MVMR sensitivity estimate, CATER-MR reports the following identification diagnostics:

- `exposure_corr` is supplied;
- the covariance-weighted conditional-F iteration converges;
- target conditional F is at least `min_cond_F` (default 10);
- the LD-weighted exposure design is full rank;
- its condition number is below `mvmr_max_condition`.

`mvmr_max_r2` is now an **optional** user-specified sensitivity gate and defaults to `NULL`. Residual LD is already represented explicitly by signed \(R\) in GLS, so CATER-MR no longer imposes an arbitrary default \(r^2<0.01\) rule on an LD-aware estimator.

The network coefficient is reported as sensitivity output. `exposure_corr` is used for conditional-strength assessment; it does not turn the standard GLS point estimator into a weak-instrument or sample-overlap bias-corrected estimator.

## Diagnostics

For each target CATER-MR reports:

- number of COJO-selected cis/trans signals **before** outcome harmonization;
- number of LD-allele-aligned signals and final MR-usable cis/trans IVs;
- `min_F`, `mean_F`, `effective_F`;
- trans information fraction;
- TF-anchor evidence;
- one-hop sibling pleiotropy results, requested/tested/missing actual-IV counts, and screen completeness;
- cis/trans heterogeneity;
- per-TF-locus incremental information;
- leave-one-TF-locus-out effect change;
- local-MVMR status, rank, weighted-design condition number, conditional F/convergence, IV-union coverage and residual maximum \(r^2\) where available.

Outcome-based quantities are never used to select trans instruments or sibling exposures.

### Primary-result semantics

Trans evidence is promoted conservatively:

- no trans IVs: use the cis estimate when available;
- default `primary_policy="cis_anchor"`: use cis MR as the primary anchor whenever a valid cis estimate is available;
- `primary_policy="screened_cater"` is an explicit research opt-in: a complete screen with no detected sibling co-perturbation can promote the combined estimate;
- detected sibling co-perturbation: keep cis as primary when available and report local MVMR as sensitivity;
- sibling screen disabled or incomplete: fall back to cis when available;
- detected/unexcluded trans pleiotropy with no valid cis fallback: leave `primary_model`, `primary_beta`, `primary_se`, `primary_p`, and therefore `primary_q` unset.

Thus a deliberately unresolved combined/trans sensitivity estimate remains in `cater_mr_results.tsv` but is never silently copied into the primary fields.

## Multiple testing

COJO significance controls instrument selection only.

Across all tested target genes, CATER-MR separately applies BH-FDR to:

- each long-form MR model (`cis`, `trans`, `combined`, `network`);
- the final primary target-level effect.

## Minimal usage

```r
source("CATER_MR.R")

res <- cater_mr(
  grn = microglia_grn,
  eqtl_dir = "/data/microglia/full_eqtl",
  outcome = outcome_gwas,
  gene_annotation = gene_annotation,
  ld_bfile = "/data/ld/microglia_donors",
  manc_cojo_bin = "manc_cojo",
  cell_type = "Microglia",
  trait = "Disease",
  qtl_n = 982,
  outdir = "microglia_CATER"
)
```

If local MVMR is allowed to become the primary model, provide donor-level gene-expression correlations for the possible exposures:

```r
res <- cater_mr(
  ...,
  exposure_corr = gene_expression_correlation,
  min_cond_F = 10
)
```

## Output

```text
<CATER_MR_results>/
├── cater_mr_results.tsv
├── cater_mr_target_summary.tsv
├── instruments/
├── cojo/
├── diagnostics/
├── mechanism/
├── pleiotropy/
└── cojo_mvmr/
```

`cater_mr_results.tsv` is long-form model output. `cater_mr_target_summary.tsv` contains the final target-level primary model, BH-FDR and design diagnostics.

## Status semantics

Statistical non-identifiability is not treated as a software error. Examples include:

```text
NO_EQTL_FILE
NO_TARGET_ANNOTATION
NO_CANDIDATE_SNP
NO_COJO_SIGNAL
NO_LD_ALLELE_MATCH
NO_HARMONIZED_IV
LD_SINGULAR
TRANS_PLEIOTROPY_UNRESOLVED
TRANS_PLEIOTROPY_UNRESOLVED_NO_CIS
TRANS_UNSCREENED_CIS_FALLBACK
TRANS_UNSCREENED_NO_CIS
SIBLING_SCREEN_INCOMPLETE_CIS_FALLBACK
SIBLING_SCREEN_INCOMPLETE_NO_CIS
OK_CIS_ONLY
OK_CATER
OK_NETWORK_ADJUSTED
```

## Design validation

`tests/smoke.R` checks mathematical invariants and regression failure modes without requiring Manc-COJO itself, including:

- one-hop region assignment;
- Manc-COJO LD block parsing;
- allele orientation;
- GIVW reduction to ordinary IVW when \(R=I\);
- \(F_{\mathrm{eff}}\) reduction to mean F when \(R=I\);
- one-SNP sibling omnibus reduction to \(z^2\);
- cis/trans heterogeneity under equal slopes;
- trans-only information fraction equals 1 rather than `NA`;
- covariance-weighted conditional F differs from unweighted OLS when standard errors differ;
- conditional F uses signed LD and residual degrees of freedom `m-p+1`;
- partial actual-IV sibling coverage is marked incomplete;
- unresolved or incompletely screened trans sensitivity estimates never populate `primary_*`;
- exact recovery of known coefficients in a deterministic local-MVMR example.

## Scope deliberately not added

V0.5 still does not make the following default components:

- two-hop/recursive network expansion;
- parent-TF MVMR exposure;
- MR-BMA;
- whole-GRN joint MVMR;
- fine-mapping/colocalization as a hard gate;
- MR-link-2 as the primary estimator;
- outcome-driven instrument selection.

Those remain ablation/sensitivity extensions.

## References

- Wang X, Wang Y, Visscher PM, Wray NR, Yengo L. Multi-ancestry conditional and joint analysis (Manc-COJO). bioRxiv. 2026.
- Yang J, Ferreira T, Morris AP, et al. Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nature Genetics*. 2012;44:369-375.
- Sanderson E, et al. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. *International Journal of Epidemiology*. 2019;48:713-727.

Manc-COJO tutorial: https://light156.github.io/multi-ancestry-COJO-docs/tutorial/
