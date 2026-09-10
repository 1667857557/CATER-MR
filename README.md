# CATER-MR

**Cis And Trans eQTLs guided by Regulatory networks for drug-target Mendelian Randomization**

CATER-MR is a target-centric framework for combining cell-type-specific cis-eQTL signals with biologically constrained trans-eQTL signals from one-hop upstream TF loci.

## V0.2 design

The first working version deliberately follows an Occam-style workflow:

```text
cell-type GRN
    -> target cis + one-hop TF genomic regions
    -> target full-summary eQTL
    -> Manc-COJO conditional/joint signal selection
    -> cis / trans / combined MR
```

There is one main script: [`CATER_MR.R`](CATER_MR.R).

V0.2 **does not** yet implement the optional QTL-triggered sibling MVMR layer, fine-mapping/colocalization, MR-BMA, or whole-network models. Those are intentionally deferred until real-data ablations show that they are necessary.

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

CATER-MR accepts common aliases for these columns:

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
N          # required by COJO, or pass qtl_n=<constant donor N>
```

`A1` must be the effect allele and `EAF/freq` must be the frequency of A1.

### 3. Outcome object

`outcome` is an in-memory `data.frame`. Minimal columns are:

```text
SNP, A1, A2, BETA, SE
```

`P`, `EAF`, chromosome and position are optional for the current MR step.

## Manc-COJO backend

CATER-MR uses **Manc-COJO** (`light156/multi-ancestry-COJO`) instead of GCTA-COJO.

Manc-COJO is a C++ implementation of conditional and joint analysis. It supports both single-ancestry and multi-ancestry analyses. The current CATER-MR input contract provides one cell-type-specific eQTL summary file per gene, so V0.2 intentionally uses Manc-COJO in **single-cohort mode**. Multi-cohort eQTL input can be added later without changing the CATER-MR conceptual framework.

Install Manc-COJO following its official documentation, for example through Bioconda, and make sure the executable is available as:

```bash
manc_cojo --help
```

The COJO summary input uses the required first eight columns:

```text
SNP A1 A2 freq b se p N
```

The LD reference is supplied as an ancestry-matched PLINK `bed/bim/fam` prefix through `ld_bfile`. Ideally this is the genotype sample used for the eQTL mapping.

### Why two Manc-COJO calls per target?

CATER-MR first performs independent-signal selection:

```bash
manc_cojo \
  --bfile <LD> \
  --cojo-file <full target summary> \
  --extract <cis + one-hop TF candidate SNPs> \
  --cojo-slct \
  --out <target.select>
```

This produces `<target.select>.jma.cojo`.

Manc-COJO only writes `.ldr.cojo` when `--output-all` is enabled. Enabling `--output-all` during stepwise selection can also generate a large `.cma.cojo`, which CATER-MR does not need. Therefore, when more than one SNP is selected, CATER-MR follows the Manc-COJO tutorial and runs a second **joint-only** call on the selected SNPs:

```bash
manc_cojo \
  --bfile <LD> \
  --cojo-file <full target summary> \
  --extract <target.select.jma.cojo> 2 header \
  --cojo-joint \
  --output-all \
  --out <target.joint>
```

The selected-SNP LD matrix is then read from `<target.joint>.ldr.cojo`. For SNPs on different chromosomes, CATER-MR combines the Manc-COJO chromosome blocks as a block-diagonal LD matrix, with cross-chromosome correlations set to zero.

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
  cojo_threads = 4L,
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

1. Identify the target cis region around `TSS_X`.
2. Find only direct `TF -> X` GRN parents and define their local gene regions.
3. Read `X.txt.gz` once.
4. Build the candidate SNP set as the union of X-cis and one-hop parent-TF loci; overlapping TF-locus SNPs inside X-cis remain classified as cis.
5. Run one Manc-COJO stepwise selection on that candidate set using the full X summary as the COJO summary input.
6. Classify selected independent sentinels as cis or GRN-trans by genomic position.
7. If multiple SNPs are selected, run Manc-COJO joint mode only on those SNPs to obtain the compact signed LD matrix.
8. Harmonize selected SNPs against the outcome object; palindromic SNPs are dropped by default.
9. Estimate cis-only, trans-only and combined effects using IVW/GIVW.

The primary combined model is:

\[
\hat\theta_X =
\frac{\gamma_X^T\Omega_Y^{-1}\Gamma_Y}
     {\gamma_X^T\Omega_Y^{-1}\gamma_X},
\qquad
\Omega_Y=D_Y R D_Y.
\]

The SNP-to-X effects used in MR are the original marginal eQTL effects. Manc-COJO is used for **independent-signal selection**, not to manufacture or rescale eQTL effects.

## Output

```text
<CATER_MR_results>/
├── cater_mr_results.tsv
├── instruments/
│   └── <SYMBOL>.tsv
└── cojo/
    ├── <SYMBOL>.sumstat
    ├── <SYMBOL>.candidate.snplist
    ├── <SYMBOL>.select.jma.cojo
    ├── <SYMBOL>.select.log
    ├── <SYMBOL>.joint.jma.cojo
    ├── <SYMBOL>.joint.ldr.cojo
    └── stdout/stderr logs
```

`cater_mr_results.tsv` contains three rows per successfully analyzed target (`cis`, `trans`, `combined`) with:

```text
target, model, n_iv, beta, se, p, Q, Q_p, mean_F, min_F, status
```

Missing QTL files, missing annotation, failed COJO, absent signals and failed harmonization are returned as explicit statuses rather than silently converted to zero effects.

## Important V0.2 assumptions

- SNP IDs and alleles must be compatible between QTL, outcome and LD reference.
- GRN topology is used only to define candidate one-hop trans loci. GRN edge weights never multiply eQTL beta values.
- `cojo_p` is a single threshold in this first implementation. Separate cis/trans thresholds can be added only if real-data requirements justify the extra branch.
- Manc-COJO defaults (`GCTA` selection/effect modes) are retained unless future ablations justify another mode.
- A COJO sentinel represents a conditionally independent association signal; CATER-MR does not claim that it is the causal variant.
- V0.2 reports raw cis+GRN-trans MR. Structured TF-regulon pleiotropy is the next planned module and should be added only after the simplest pipeline is benchmarked on real data.

## Smoke test

The repository includes base-R smoke tests for input normalization, one-hop region assignment, allele harmonization, GIVW, and parsing Manc-COJO chromosome-block `.ldr.cojo` output:

```bash
Rscript tests/smoke.R
```

The smoke test does not require a Manc-COJO executable or an LD panel; full end-to-end analysis does.

## COJO references

- Wang X, Wang Y, Visscher PM, Wray NR, Yengo L. **Multi-ancestry conditional and joint analysis (Manc-COJO) applied to GWAS summary statistics.** bioRxiv. 2026.01.30.702783.
- Yang J, Ferreira T, Morris AP, et al. **Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits.** *Nature Genetics*. 2012;44:369-375. doi:10.1038/ng.2213.

Manc-COJO documentation: https://light156.github.io/multi-ancestry-COJO-docs/tutorial/
