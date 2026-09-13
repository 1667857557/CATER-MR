# scMORE → CATER-MR cell-type GRN construction

CATER-MR can construct its required **cell-type-specific, direct one-hop GRN input** from a processed single-cell multiome Seurat object while delegating the actual regulatory-network inference to the audited upstream `scMORE::createRegulon()` implementation.

## CATER-MR GRN contract

The MR estimator does not need a TF–peak–gene row table as its graph. For each cell type it needs a direct topology

```text
TF  Target
```

with one unique direct `TF -> Target` edge, because for target `X` the candidate region is

\[
\mathcal R_X=L_X\cup\bigcup_{T\in P_1(X)}L_T,
\]

where `P1(X)` is the set of direct upstream TFs.

To make the GRN directly consumable by `cater_mr()` without another annotation step, the generated table contains:

```text
TF
Target
TF_chr
TF_tss
Target_chr
Target_tss
cell_type
```

Only this direct one-hop topology is used by CATER-MR. No two-hop or recursive expansion is performed, and scMORE edge scores never rescale eQTL effects.

## Two-layer architecture

Upstream scMORE and CATER-MR have different output contracts, so the adapter deliberately keeps two layers:

```text
processed multiome Seurat
    -> subset one cell type
    -> scMORE::createRegulon() unchanged
         -> raw TF/Target/Regions/Pval-or-Corr-or-Gain evidence
    -> deterministic CATER adapter
         -> unique direct TF -> Target edges
         -> remove TF == Target self-loops from the CATER trans graph
         -> append hg38 TF/Target chr + TSS
    -> cater_mr()
```

The raw scMORE output is never overwritten. It is retained in `scmore_outputs` and `edge_evidence`. The CATER-facing graph is stored separately in `grns`.

This separation is important: collapsing repeated TF–Target rows is a graph-contract transformation, not a new statistical aggregation. CATER-MR does **not** invent a combined P value, correlation, or edge weight across scMORE peak/evidence rows.

## Why per-cell-type refitting is used

The audited upstream `createRegulon()` constructs a global GRN for the cells supplied to it. CATER-MR therefore defines the cell population first and calls that same upstream engine separately for each cell type:

```text
all processed cells
  -> cells of type c
  -> scMORE::createRegulon(cells of type c)
  -> GRN_c
```

The internal scMORE/Pando computation is unchanged. This is a CATER-MR adapter around the upstream global-per-input engine; it is not a claim that scMORE's top-level `scMore()` function natively refits one Pando network per cell type.

## Why `scMore()` / `regulon2disease()` are not used

The top-level scMORE workflow adds cell-type specificity and GWAS/MAGMA-derived trait relevance. CATER-MR uses the GRN to gate trans-instrument search, so outcome/GWAS information must not determine which GRN edges enter the primary MR analysis.

Therefore the CATER-MR GRN layer stops at outcome-independent `createRegulon()` inference. CTS/TRS or disease relevance can be used later as external annotation or sensitivity analysis, but not as the primary trans-IV gate.

## Audited upstream revision

```text
repository: mayunlong89/scMORE
commit:     f614736b9f49631471b0f18dfa415ee35e4d7b66
version:    2.0.0
```

Recommended reproducible install:

```r
remotes::install_github(
  "mayunlong89/scMORE@f614736b9f49631471b0f18dfa415ee35e4d7b66"
)
```

`strict_upstream = TRUE` verifies the installed GitHub `RemoteSha`.

## Upstream calculation preserved

For each cell-type subset, CATER-MR calls `scMORE::createRegulon()` with its audited public calculation defaults except for the candidate regulatory-region prior:

```r
n_targets = 5
peak2gene_method = "Signac"
infer_method = "glm"
tss_upstream = 100000
tss_downstream = 0
exclude_exon_regions = TRUE
```

By default, CATER-MR loads `phastConsElements20Mammals.UCSC.hg38` and `SCREEN.ccRE.UCSC.hg38` from Pando, combines them with `GenomicRanges::union()`, and passes that `GRanges` object as `conserved_regions`. Consequently, the downstream `Pando::initiate_grn(regions=...)` call uses the conserved-plus-SCREEN union rather than scMORE's upstream phastCons-only default.

To reproduce the audited scMORE upstream region behavior instead, explicitly set:

```r
conserved_regions = NULL
```

CATER-MR then omits the `conserved_regions` argument when calling `scMORE::createRegulon()`, allowing scMORE itself to evaluate its `phastConsElements20Mammals.UCSC.hg38` default.

The upstream sequence remains:

```text
FindVariableFeatures(RNA)
 -> Pando::initiate_grn()
 -> Pando::find_motifs()
 -> Pando::infer_grn()
 -> Pando::find_modules()
 -> Pando::NetworkModules()
 -> scMORE internal extract_grn()
 -> n_targets TF filter
```

No Pando regression, motif score, module statistic, P value, correlation, or gain is recalculated in CATER-MR.

## Gene coordinates

CATER-MR requires one genomic TSS for every TF and target node because target cis and parent-TF loci are defined around those TSS values.

By default the adapter derives:

```text
symbol  chr  tss
```

from `Signac::Annotation(single_cell[["peaks"]])` using strand-aware TSS:

\[
TSS=
\begin{cases}
start,& strand=+\\
end,& strand=-
\end{cases}
\]

If a symbol maps to multiple distinct chromosome/TSS values, the adapter stops rather than silently choosing one transcript. In that case provide an explicit gene-level annotation:

```r
gene_annotation = data.frame(
  symbol = ...,
  chr = ...,
  tss = ...
)
```

The strict default also errors if a CATER edge lacks TF or target coordinates. `missing_coordinate="drop"` is available only as an explicit opt-in and records the number of discarded edges.

## Genome-build safeguard

scMORE motif inference uses `BSgenome.Hsapiens.UCSC.hg38`. Therefore explicitly annotated non-hg38 builds are rejected. If annotation genome metadata is absent, the build is recorded as unknown rather than guessed; the user remains responsible for ensuring the processed object is GRCh38/hg38.

## Self loops

`TF == Target` rows are preserved in raw scMORE evidence but are removed from the default CATER-ready graph:

```r
drop_self_loops = TRUE
```

A self locus is the target's own cis locus, not an independent one-hop upstream trans locus. Keeping it as a CATER parent TF would add no new trans region and could complicate sibling-path interpretation.

## Usage

```r
source("SCMORE_GRN.R")
source("CATER_MR.R")

fit <- cater_build_scmore_grn(
  single_cell = A,
  celltype_col = "cell_type",
  outdir = "scMORE_celltype_GRN"
)

microglia_grn <- cater_get_scmore_grn(fit, "Microglia")
```

`microglia_grn` is already a valid CATER-MR input, including coordinates:

```r
res <- cater_mr(
  grn = microglia_grn,
  eqtl_dir = "/data/microglia/full_eqtl",
  outcome = outcome_gwas,
  ld_bfile = "/data/ld/microglia_donors",
  cell_type = "Microglia",
  trait = "Disease"
)
```

No separate `gene_annotation` argument is required in this path because `cater_mr()` extracts coordinates from `TF_chr/TF_tss/Target_chr/Target_tss`.

## Returned object

```text
grns
  per-cell-type CATER-ready direct one-hop graphs

edge_evidence
  raw scMORE GRN evidence rows + cell_type only

scmore_outputs
  complete raw createRegulon() list for each cell type

gene_annotation
  standardized symbol/chr/tss table used for the graph contract

summary
  n_cells
  n_raw_rows
  n_unique_scMORE_edges
  n_cater_edges
  n_tfs
  n_targets
  n_self_loops_dropped
  n_missing_coordinate_edges
  status

provenance
celltype_source
gene_annotation_source
createRegulon_args
mode
```

Retrieve other layers with:

```r
cater_get_scmore_grn(fit, "Microglia")                  # CATER-ready graph
cater_get_scmore_grn(fit, "Microglia", evidence=TRUE)   # raw edge evidence table
cater_get_scmore_grn(fit, "Microglia", raw=TRUE)        # full createRegulon output
```

## Output files

For each cell type, filenames include a deterministic unique suffix so labels that sanitize to the same text cannot overwrite each other.

```text
scMORE_raw_<celltype>__NNN.rds
scMORE_edge_evidence_<celltype>__NNN.tsv
CATER_grn_<celltype>__NNN.tsv
CATER_gene_annotation.tsv
scMORE_grn_summary.tsv
scMORE_provenance.rds
```

Empty scMORE GRNs are recorded as `EMPTY_GRN` rather than causing a scalar-column assignment error. Recorded upstream error messages are sanitized before TSV output.

## What CATER-MR ultimately consumes

For a target `X`, only direct edges in the generated graph are used:

\[
P_1(X)=\{T:(T\to X)\in E_c\}.
\]

Then the estimator constructs:

\[
\mathcal R_X=L_X\cup\bigcup_{T\in P_1(X)}L_T,
\]

queries the target's own full-summary eQTL in these regions, and performs one target-level Manc-COJO selection. The scMORE peak regions are biological evidence for how the direct edge was inferred; they are **not** substituted for the TF genomic locus used by CATER-MR.

That distinction is the reason the adapter outputs both raw scMORE evidence and a separate direct, coordinate-complete CATER graph.
