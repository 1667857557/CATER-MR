# scMORE GRN adapter

## Upstream

Audited implementation:

```text
repository: mayunlong89/scMORE
version: 2.0.0
commit: f614736b9f49631471b0f18dfa415ee35e4d7b66
```

## Input

`cater_build_scmore_grn()` requires a processed Seurat object with:

```text
RNA assay
peaks ChromatinAssay
gene annotation on peaks
cell-type labels from Idents() or celltype_col
hg38 / GRCh38 coordinates
```

Optional `gene_annotation` columns:

```text
symbol
chr
tss
```

## Functions

| Function | Description |
|---|---|
| `cater_build_scmore_grn()` | Refits `scMORE::createRegulon()` per requested cell type and returns CATER-ready direct TF-target graphs |
| `cater_get_scmore_grn()` | Extracts one CATER-ready GRN, raw scMORE output, or edge-evidence table |
| `print.cater_scmore_grn()` | Prints build summary |

## Output

CATER-ready GRN columns:

```text
TF
Target
TF_chr
TF_tss
Target_chr
Target_tss
cell_type
```

`cater_scmore_grn` components:

```text
grns
edge_evidence
scmore_outputs
gene_annotation
summary
provenance
celltype_source
gene_annotation_source
createRegulon_args
mode
```

## Minimal run

```r
source("SCMORE_GRN.R")

x <- cater_build_scmore_grn(
  single_cell = obj,
  celltype_col = "cell_type",
  cell_types = c("B", "T"),
  strict_upstream = TRUE,
  outdir = "scmore_grn"
)

grn_B <- cater_get_scmore_grn(x, "B")
```

## References

1. Ma Y, Yao Y, Zhou Y, et al. *Nature Aging.* 2026;6:270-289. https://doi.org/10.1038/s43587-025-01027-5
2. Fleck JS, Jansen SMJ, Wollny D, et al. *Nature.* 2023;621:365-372. https://doi.org/10.1038/s41586-022-05279-8
