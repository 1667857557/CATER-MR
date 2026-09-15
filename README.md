# CATER-MR

**Cis And Trans eQTLs guided by Regulatory networks for drug-target Mendelian randomization.**

## Functions

| Function | Description |
|---|---|
| `cater_mr()` | GRN-constrained cis/trans/combined MR with signed-LD estimation and optional MVMR sensitivity |
| `cater_prepare_eqtl_table()` | Materializes merged full cis-eQTL summaries to `<GENE>.txt.gz` files |
| `cater_mr_from_table()` | Runs `cater_mr()` from a merged full cis-eQTL table |
| `cater_build_scmore_grn()` | Builds cell-type GRNs with audited `scMORE::createRegulon()` |
| `cater_get_scmore_grn()` | Extracts a CATER-ready cell-type GRN or associated scMORE evidence |

## Interfaces

- [`EQTL_INPUT.md`](EQTL_INPUT.md)
- [`SCMORE_GRN.md`](SCMORE_GRN.md)
- [`MATHEMATICS.md`](MATHEMATICS.md)

## References

1. Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC. *Genet Epidemiol.* 2017;41:714-725. https://doi.org/10.1002/gepi.22077
2. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. *Int J Epidemiol.* 2019;48:713-727. https://doi.org/10.1093/ije/dyy262
3. Fleck JS, Jansen SMJ, Wollny D, et al. *Nature.* 2023;621:365-372. https://doi.org/10.1038/s41586-022-05279-8
4. Ma Y, Yao Y, Zhou Y, et al. *Nature Aging.* 2026;6:270-289. https://doi.org/10.1038/s43587-025-01027-5
