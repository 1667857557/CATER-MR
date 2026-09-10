# scMORE cell-type GRN construction for CATER-MR

CATER-MR can construct outcome-independent, cell-type-specific GRNs directly from a processed single-cell multiome Seurat object while delegating the actual GRN inference to the upstream `scMORE::createRegulon()` implementation.

## Why this adapter is needed

The current scMORE upstream code distinguishes two stages:

1. `createRegulon(single_cell, ...)` constructs a **global TF-peak-gene GRN for the cells supplied to the function** using Pando.
2. The top-level `scMore()` workflow subsequently computes cell-type specificity and integrates GWAS/MAGMA information to identify trait-relevant cell-type-specific eRegulons.

CATER-MR must not use GWAS/outcome information to construct or select its GRN, because the GRN is part of the exposure-side instrument-selection mechanism. Therefore CATER-MR does **not** call `scMore()` or `regulon2disease()` when constructing GRNs.

Instead, the CATER adapter defines the cell population first and then calls the unmodified upstream GRN engine:

```text
processed multiome Seurat
    -> split cells by cell type
    -> scMORE::createRegulon(cell-type subset)
    -> one GRN per cell type
    -> CATER-MR
```

The only CATER-specific operation before GRN inference is the cell subset. No scMORE/Pando regression, motif, module, filtering or edge-score formula is reimplemented in CATER-MR.

## Audited upstream revision

The adapter was audited against:

```text
repository: mayunlong89/scMORE
commit:     f614736b9f49631471b0f18dfa415ee35e4d7b66
version:    2.0.0
```

For strict reproducibility install exactly that GitHub revision:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github(
  "mayunlong89/scMORE@f614736b9f49631471b0f18dfa415ee35e4d7b66"
)
```

`cater_build_scmore_grn(..., strict_upstream = TRUE)` verifies the installed package `RemoteSha`. Set `strict_upstream = FALSE` only after independently validating a different upstream revision.

## Exact scMORE computation retained

For each cell-type subset, CATER-MR calls:

```r
scMORE::createRegulon(
  single_cell,
  n_targets = 5,
  peak2gene_method = "Signac",
  infer_method = "glm",
  tss_upstream = 100000,
  tss_downstream = 0,
  exclude_exon_regions = TRUE
)
```

When `conserved_regions` is not explicitly supplied, CATER-MR omits that argument so scMORE evaluates its own upstream default `phastConsElements20Mammals.UCSC.hg38`.

Inside the audited upstream function the sequence is:

```text
Seurat::FindVariableFeatures(assay = "RNA")
    -> Pando::initiate_grn(
         peak_assay = "peaks",
         rna_assay = "RNA",
         exclude_exons = exclude_exon_regions,
         regions = conserved_regions)
    -> Pando::find_motifs(
         pfm = scMORE motifs,
         genome = BSgenome.Hsapiens.UCSC.hg38)
    -> Pando::infer_grn(
         peak_to_gene_method = peak2gene_method,
         method = infer_method,
         upstream = tss_upstream,
         downstream = tss_downstream,
         alpha = 0.5,
         family = "gaussian",
         adjust_method = "fdr",
         scale = FALSE,
         verbose = TRUE)
    -> Pando::find_modules(
         p_thresh = 0.1,
         nvar_thresh = 2,
         min_genes_per_module = 1,
         rsq_thresh = 0.05)
    -> Pando::NetworkModules()
    -> scMORE internal extract_grn()
    -> retain TFs with at least n_targets GRN rows
```

CATER-MR does not change these fixed internal values.

## Required Seurat object state

The processed multiome Seurat object must satisfy scMORE's actual code assumptions:

- human GRCh38/hg38 data;
- RNA assay named exactly `RNA`;
- chromatin assay named exactly `peaks`;
- the `peaks` ChromatinAssay already has gene annotation assigned with `Annotation(object[["peaks"]])`;
- cell-type labels are available either in a metadata column or in `Idents(object)`.

The adapter intentionally does not rename assays, normalize RNA, rebuild peaks, add annotations, or alter the object silently.

## Usage

```r
source("SCMORE_GRN.R")

scmore_grns <- cater_build_scmore_grn(
  single_cell = A,
  celltype_col = "cell_type",
  outdir = "scMORE_celltype_GRN"
)

scmore_grns$summary
```

Use one fitted cell-type GRN in the CATER estimator:

```r
source("CATER_MR.R")

microglia_grn <- cater_get_scmore_grn(
  scmore_grns,
  cell_type = "Microglia"
)

res <- cater_mr(
  grn = microglia_grn,
  eqtl_dir = "/data/microglia/full_eqtl",
  outcome = outcome_gwas,
  gene_annotation = gene_annotation,
  ld_bfile = "/data/ld/microglia_donors",
  cell_type = "Microglia",
  trait = "Disease"
)
```

To fit only selected populations:

```r
scmore_grns <- cater_build_scmore_grn(
  A,
  celltype_col = "cell_type",
  cell_types = c("Microglia", "Astrocyte")
)
```

If `celltype_col = NULL`, the function uses `Idents(A)`.

## Returned object

`cater_build_scmore_grn()` returns a `cater_scmore_grn` list containing:

```text
grns            # per-cell-type GRN tables; ready for cater_mr()
scmore_outputs  # raw upstream createRegulon() output for each cell type
summary         # n_cells / n_edges / n_TFs / n_targets / status
provenance      # scMORE, Pando, Seurat and Signac versions + upstream SHA
celltype_source # metadata column or Idents
createRegulon_args
mode
```

The raw upstream object is preserved separately. The `grns` tables only add a `cell_type` metadata column; `TF`, `Target`, `Regions`, and `Pval`/`Corr`/`Gain` are not recalculated.

## Failure semantics

Default behavior is strict:

```r
on_error = "stop"
```

If one cell type fails inside scMORE/Pando, CATER-MR stops and exposes the original error. It never substitutes another GRN algorithm.

For large atlases, explicit partial collection is available:

```r
on_error = "record"
```

Failed cell types are then marked `SCMORE_FAILED` in the summary and have no GRN object.

Cell types are processed sequentially. No outer parallel layer is added because Pando/scMORE may perform expensive internal work and nested parallelism can amplify memory use.

## Important upstream implementation boundary

The audited `createRegulon()` documentation mentions several inference methods, but the current upstream internal `extract_grn()` contains explicit extraction branches for:

- `glm`;
- `cv.glmnet` / `glmnet`;
- `xgb`.

Other `infer_method` values currently reach the upstream `Invalid inference method` branch during extraction. CATER-MR deliberately does not patch this behavior; it forwards the requested method to scMORE unchanged.

## Why CTS/GWAS-based scMORE filtering is not used here

The top-level scMORE trait workflow integrates cell-type specificity with GWAS/MAGMA-derived gene relevance. That is appropriate for scMORE's trait-regulon objective, but using outcome-derived information to decide which CATER-MR trans loci or GRN edges enter MR would make instrument/network selection outcome-informed.

Therefore the CATER-MR GRN-construction layer stops at the outcome-independent `createRegulon()` output. Any later use of scMORE CTS/TRS should be treated as external annotation or sensitivity analysis, not as the primary CATER instrument-selection gate.
