# Pre-COJO LD consistency diagnosis

CATER-MR performs an LD/summary-statistic consistency gate before every Manc-COJO selection step by default.

Trans candidates are defined by local parent-TF windows. TF windows on the same chromosome that overlap are merged into one physical locus component before candidate SNP assignment (for example, overlapping A and B windows share `locus_id=TF:A+B`). The original windows are retained for `parent_tf`, so a SNP can still be attributed to A, B, or A;B according to its actual position while LD/SuSiE treats the full overlapping component as one locus.

## Method

For each physical candidate cis/TF locus, CATER-MR:

1. forms local diagnostic loci. The target cis locus is kept as its own locus. Trans candidates are defined by parent-TF windows; when two or more TF windows overlap, SNPs annotated to the shared TF set connect those TFs and the full connected component is merged into one physical LD-diagnosis locus. Thus `TF:A`, `TF:A;B`, and `TF:B` are analysed together rather than as three artificial loci;
2. extracts candidate variants for that local locus from the configured PLINK LD reference;
3. computes a signed allele-count correlation matrix with PLINK 1.9 (`--r square`);
4. aligns eQTL z-score signs to direct PLINK A1/A2 same/swap coding; strand-complement relationships are recorded for audit but excluded because the downstream CATER-MR COJO/MR allele contract also requires direct A1/A2 agreement;
5. estimates the SuSiE-RSS consistency parameter with `susieR::estimate_s_rss()` and obtains conditional z-score diagnostics with `susieR::kriging_rss()`;
6. removes variants satisfying `logLR > 2 & abs(z) > 2` before Manc-COJO runs.

The original `parent_tf` and `locus_id` fields are retained for biological attribution. The diagnostic output additionally records `ld_locus_id`, which is the merged physical locus used for PLINK/SuSiE-RSS. This merge affects only LD diagnosis; it does not redefine the downstream TF attribution or the MR model.

Variants absent from the LD reference, variants with non-finite LD rows, and variants whose allele pair cannot be reconciled with the direct A1/A2 contract are also removed and explicitly reported. A one-variant locus cannot be conditionally diagnosed and is retained with status `NOT_DIAGNOSABLE_SINGLETON`.

The PLINK matrix is checked for dimensions, finite entries, symmetry, unit diagonal, and correlation range. CATER-MR deliberately does not impose an additional strict positive-semidefinite eigenvalue gate before SuSiE-RSS, because `susieR` performs the RSS eigenvalue handling itself and text-formatted PLINK matrices can contain small rounding artifacts.

PLINK subset generation is fail-closed: a non-zero PLINK exit status is treated as an error except for the explicit `No variants remaining` case, which is reported as LD-reference missing.

This is a QC gate only. It does not use SuSiE PIP or credible sets for instrument selection; Manc-COJO remains the instrument-selection method.

## Dependencies and options

With the default `enable_ld_diagnosis=TRUE`, PLINK 1.9 and the R package `susieR` are required. `plink_bin` can be an executable name on `PATH` or an explicit path. The default detection thresholds are `ld_diag_loglr=2` and `ld_diag_abs_z=2`.

Set `enable_ld_diagnosis=FALSE` only to reproduce the legacy COJO path without this QC gate.

## Output

For target `GENE`, the full SNP-level diagnostic table is written to:

`<outdir>/cojo/GENE.ld_diagnosis.tsv`

The target summary also reports the number removed by each major reason and the maximum estimated SuSiE-RSS consistency parameter across diagnosed loci.

Reference implementation used for the diagnostic logic: xinhe-lab/mapgen, commit `393e66fe62442d2c499e4cab10a7f48309da216a`, `R/LD_diagnosis.R`.
