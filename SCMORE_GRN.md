# scMORE GRN adapter

## Upstream

```text
repository: mayunlong89/scMORE
version: 2.0.0
commit: f614736b9f49631471b0f18dfa415ee35e4d7b66
```

## Functions

| Function | Description |
|---|---|
| `cater_build_scmore_grn()` | Refits `scMORE::createRegulon()` per cell type and returns CATER-ready direct TF-target graphs |
| `cater_get_scmore_grn()` | Extracts one CATER-ready GRN, raw scMORE output, or edge-evidence table |
| `print.cater_scmore_grn()` | Prints GRN-build metadata |

## Input schema

```text
Seurat object
RNA assay
peaks ChromatinAssay
peak gene annotation
cell-type labels
hg38 / GRCh38 coordinates
```

Optional gene annotation:

```text
symbol
chr
tss
```

## Output schema

```text
TF
Target
TF_chr
TF_tss
Target_chr
Target_tss
cell_type
```

Object components:

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

## References

1. Ma Y, Yao Y, Zhou Y, et al. *Nature Aging.* 2026;6:270-289. https://doi.org/10.1038/s43587-025-01027-5
2. Fleck JS, Jansen SMJ, Wollny D, et al. *Nature.* 2023;621:365-372. https://doi.org/10.1038/s41586-022-05279-8
