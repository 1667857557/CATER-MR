# Combined full-summary eQTL table input

CATER-MR supports the original per-gene storage contract and, through `CATER_EQTL_INPUT.R`, a single merged cell-type-specific full-summary eQTL table.

## Two equivalent storage layouts

### Existing per-gene layout

```text
eqtl_dir/
  TP53.txt.gz
  EGFR.txt.gz
  STAT3.txt.gz
  ...
```

Each file contains genome-wide full-summary eQTL statistics for that gene.

### New merged-table layout

```text
GENE   SNP   CHR   BP   A1   A2   EAF   BETA   SE   P   N
TP53   ...
TP53   ...
EGFR   ...
STAT3  ...
...
```

`GENE` identifies which molecular exposure the QTL row belongs to. Common alternatives such as `gene`, `gene_symbol`, `symbol`, and `SYMBOL` are detected case-insensitively, or another name can be supplied with `gene_col=`.

The mathematical requirement is simply:

\[
Q_X = Q^{merged}[GENE=X].
\]

The gene column is a storage/index field only. It does not enter Manc-COJO, GIVW, sibling screening, or local MVMR.

## Why the first implementation materializes a per-gene cache

The current core estimator already has one validated QTL access contract:

```text
read gene X -> X.txt.gz -> standardize full summary -> CATER-MR
```

To keep the estimator mathematically identical, the merged-table adapter performs:

```text
one merged table
    -> group by GENE
    -> temporary <GENE>.txt.gz cache
    -> unchanged cater_mr()
```

Therefore target Manc-COJO, parent-TF exact-IV queries, same-IV sibling screening, and sibling cis-QTL queries in triggered local MVMR all use the same core code as directory mode.

This follows the Occam-first implementation strategy: storage compatibility is added without introducing a second estimator or a second QTL parser. If real datasets make the materialization step an I/O bottleneck, the adapter can later be replaced by bgzip/tabix, Arrow/Parquet, or DuckDB indexing without changing the MR mathematics.

## Usage with an in-memory merged table

```r
source("CATER_MR.R")
source("CATER_EQTL_INPUT.R")

res <- cater_mr_from_table(
  grn = celltype_grn,
  eqtl_table = all_eqtl,
  outcome = outcome_gwas,
  ld_bfile = "/data/ld/celltype_donors",
  cell_type = "Microglia",
  trait = "Disease"
)
```

`all_eqtl` can contain:

```text
GENE SNP CHR BP A1 A2 EAF BETA SE P N
```

`P` remains optional if it can be derived from `BETA/SE`; `N` can still be supplied globally through `qtl_n=` exactly as in ordinary CATER-MR.

## Usage with one merged `.txt.gz` file

```r
res <- cater_mr_from_table(
  grn = celltype_grn,
  eqtl_table = "/data/Microglia_all_genes_eqtl.txt.gz",
  outcome = outcome_gwas,
  ld_bfile = "/data/ld/celltype_donors",
  cell_type = "Microglia"
)
```

The merged file is expected to be tab-delimited.

## Custom gene-column name

```r
res <- cater_mr_from_table(
  grn = celltype_grn,
  eqtl_table = all_eqtl,
  gene_col = "EXPOSURE",
  outcome = outcome_gwas,
  ld_bfile = ld_bfile
)
```

## Explicit cache preparation

For repeated outcomes using the same cell-type eQTL data, it is more efficient to materialize once and reuse the resulting directory:

```r
cache <- cater_prepare_eqtl_table(
  eqtl_table = "/data/Microglia_all_genes_eqtl.txt.gz",
  outdir = "/data/CATER_cache/Microglia"
)

res1 <- cater_mr(
  grn = celltype_grn,
  eqtl_dir = cache$eqtl_dir,
  outcome = outcome1,
  ld_bfile = ld_bfile
)

res2 <- cater_mr(
  grn = celltype_grn,
  eqtl_dir = cache$eqtl_dir,
  outcome = outcome2,
  ld_bfile = ld_bfile
)
```

This avoids rebuilding the per-gene cache for each outcome.

## Input validation

The adapter:

- requires a nonmissing gene identifier for every row;
- rejects gene identifiers containing path separators because the core cache is named `<GENE>.txt.gz`;
- validates each gene subset with the exact `.cater_standardize_sumstats()` parser used by the estimator;
- preserves the original SNP/effect columns when writing the cache;
- does not perform significance filtering before Manc-COJO;
- does not change beta, SE, P, EAF or N values.

## Scope

This first implementation intentionally does **not** introduce a database backend. A single very large compressed table is read into memory before being materialized. This is acceptable for the initial compatibility layer but may be inappropriate for atlas-scale merged tables. The planned optimization point is the QTL-access layer only; the downstream CATER-MR estimator should remain unchanged.
