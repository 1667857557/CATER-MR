# CATER-MR

**Cis And Trans eQTLs guided by Regulatory networks for drug-target Mendelian Randomization**

CATER-MR is a target-centric framework for combining cell-type-specific cis-eQTL signals with biologically constrained trans-eQTL signals from one-hop upstream TF loci.

## V0.1 design

The first version deliberately follows an Occam-style workflow:

```text
cell-type GRN
    -> target cis + one-hop TF genomic regions
    -> target full-summary eQTL
    -> GCTA-COJO selection within the candidate SNP set
    -> cis / trans / combined MR
```

There is one main script: [`CATER_MR.R`](CATER_MR.R).

V0.1 **does not** yet implement the optional QTL-triggered sibling MVMR layer, fine-mapping/colocalization, MR-BMA, or whole-network models. Those are intentionally deferred until real-data ablations show that they are necessary.

## Primary inputs

### 1. Cell-type-specific GRN object

`grn` may be either:

- a `data.frame` containing `TF` and `Target`; or
- a list containing such a table in `$grn`.

All unique TF and Target symbols are scanned by default. Regulation is strictly **one-hop**: for target `X`, only rows `TF -> X` define candidate trans loci.

If the GRN also contains `TF_chr`, `TF_tss`, `Target_chr`, and `Target_tss`, CATER-MR can derive gene coordinates from the GRN. Otherwise provide `gene_annotation` with:

```text
symbol  chr  tss
```

### 2. Cell-type-specific full-summary eQTL directory

Each measured gene is stored as:

```text
<eqtl_dir>/<SYMBOL>.txt.gz
```

For example:

```text
eqtl/
├── APOE.txt.gz
├── SPI1.txt.gz
├── TREM2.txt.gz
└── ...
```

CATER-MR accepts common aliases for these required columns:

```text
SNP / rsid / variant_id
CHR
BP / POS
A1 / effect_allele
A2 / other_allele
EAF / freq
BETA / b
SE
P          # optional; derived from beta/se if absent
N          # required for COJO, or pass qtl_n=<constant donor N>
```

`A1` must be the effect allele and EAF/freq must be the frequency of A1.

### 3. Outcome object

`outcome` is an in-memory `data.frame`. Minimal columns are:

```text
SNP, A1, A2, BETA, SE
```

`P`, `EAF`, chromosome and position are optional for the V0.1 MR step.

## Technical inputs required by COJO

CATER-MR V0.1 uses GCTA-COJO as the default independent-signal selector, so it also requires:

- `gcta64` (or another GCTA executable supplied via `gcta_bin`);
- an ancestry-matched PLINK LD reference supplied as the `bed/bim/fam` prefix `ld_bfile`.

The best LD reference is the QTL donor genotype sample itself. If that is unavailable, use a sufficiently large ancestry-matched reference.

CATER-MR writes the **full target-gene eQTL summary** to the GCTA COJO input and uses `--extract` to restrict model selection to target cis plus one-hop TF regions. This follows GCTA's requirement that the COJO summary input retain genome-wide summary information while an extraction list limits the region/SNP set under analysis.

## Minimal usage

```r
source("CATER_MR.R")

res <- cater_mr(
  grn = microglia_grn,
  eqtl_dir = "/data/microglia/full_eqtl",
  outcome = outcome_gwas,
  gene_annotation = gene_annotation,
  ld_bfile = "/data/ld/microglia_donors",
  gcta_bin = "gcta64",
  cis_window = 1e6,
  cojo_p = 5e-8,
  outdir = "microglia_trait_CATER"
)
```

If every QTL file has the same donor sample size but does not contain an `N` column:

```r
res <- cater_mr(
  grn = microglia_grn,
  eqtl_dir = "/data/microglia/full_eqtl",
  outcome = outcome_gwas,
  gene_annotation = gene_annotation,
  ld_bfile = "/data/ld/microglia_donors",
  qtl_n = 982
)
```

## What the script does for each target X

1. Identify target cis region around `TSS_X`.
2. Find only direct `TF -> X` GRN parents and define their local gene regions.
3. Read `X.txt.gz` once.
4. Build the candidate SNP set as the union of X-cis and one-hop parent-TF loci; overlapping TF-locus SNPs inside X-cis remain classified as cis.
5. Run one joint COJO selection on that candidate set using the full X summary as the COJO summary input.
6. Classify selected independent sentinels as cis or GRN-trans by genomic position.
7. Harmonize selected SNPs against the outcome object; palindromic SNPs are dropped by default.
8. Estimate cis-only, trans-only and combined effects. The `.jma.ldr` signed LD matrix is retained, so residual LD among selected COJO signals is handled with generalized IVW rather than pretending all selected SNPs are independent.

The primary V0.1 combined model is:

\[
\hat\theta_X =
\frac{\gamma_X^T\Omega_Y^{-1}\Gamma_Y}
     {\gamma_X^T\Omega_Y^{-1}\gamma_X},
\qquad
\Omega_Y=D_Y R D_Y.
\]

The SNP-to-X effects used in MR are the original marginal eQTL effects. COJO is used for **signal selection**, not to manufacture or rescale eQTL effects.

## Output

```text
<CATER_MR_results>/
├── cater_mr_results.tsv
├── instruments/
│   └── <SYMBOL>.tsv
└── cojo/
    ├── <SYMBOL>.ma
    ├── <SYMBOL>.candidate.snplist
    ├── <SYMBOL>.jma
    ├── <SYMBOL>.jma.ldr
    └── logs...
```

`cater_mr_results.tsv` contains three rows per successfully analyzed target (`cis`, `trans`, `combined`) with:

```text
target, model, n_iv, beta, se, p, Q, Q_p, mean_F, min_F, status
```

Missing QTL files, missing annotation, failed COJO, absent signals and failed harmonization are returned as explicit statuses rather than silently converted to zero effects.

## Important V0.1 assumptions

- SNP IDs and alleles must be compatible between QTL, outcome and LD reference.
- GRN topology is used only to define candidate one-hop trans loci. GRN edge weights never multiply eQTL beta values.
- `cojo_p` is a single threshold in V0.1. Separate cis/trans thresholds can be added after real-data requirements justify the extra branch.
- V0.1 does not claim that a COJO sentinel is the causal variant; it represents a conditionally independent association signal.
- V0.1 reports raw cis+GRN-trans MR. Structured TF-regulon pleiotropy is the next planned module and should be added only after the simplest pipeline is benchmarked on real data.

## Smoke test

The repository includes a base-R smoke test for input normalization, one-hop region assignment, allele harmonization and GIVW:

```bash
Rscript tests/smoke.R
```

The smoke test does not require GCTA or an LD panel; full end-to-end analysis does.

## Reference for COJO

Yang J, Ferreira T, Morris AP, et al. Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nature Genetics*. 2012;44:369-375. doi:10.1038/ng.2213.
