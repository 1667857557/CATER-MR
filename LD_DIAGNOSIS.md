# Pre-COJO LD consistency diagnosis

CATER-MR performs an LD/summary-statistic consistency gate before every Manc-COJO selection step by default.

## Method

For each candidate cis/trans locus, CATER-MR:

1. extracts candidate variants from the configured PLINK LD reference;
2. computes a signed allele-count correlation matrix with PLINK 1.9 (`--r square`);
3. aligns eQTL z-score signs to PLINK A1/A2 coding, including strand-complement matches;
4. estimates the SuSiE-RSS consistency parameter with `susieR::estimate_s_rss()` and obtains conditional z-score diagnostics with `susieR::kriging_rss()`;
5. removes variants satisfying `logLR > 2 & abs(z) > 2` before Manc-COJO runs.

Variants absent from the LD reference, variants with non-finite LD rows, and variants whose allele pair cannot be reconciled with the LD reference are also removed and explicitly reported. A one-variant locus cannot be conditionally diagnosed and is retained with status `NOT_DIAGNOSABLE_SINGLETON`.

This is a QC gate only. It does not use SuSiE PIP or credible sets for instrument selection; Manc-COJO remains the instrument-selection method.

## Dependencies and options

With the default `enable_ld_diagnosis=TRUE`, PLINK 1.9 and the R package `susieR` are required. `plink_bin` can be an executable name on `PATH` or an explicit path. The default detection thresholds are `ld_diag_loglr=2` and `ld_diag_abs_z=2`.

Set `enable_ld_diagnosis=FALSE` only to reproduce the legacy COJO path without this QC gate.

## Output

For target `GENE`, the full SNP-level diagnostic table is written to:

`<outdir>/cojo/GENE.ld_diagnosis.tsv`

The target summary also reports the number removed by each major reason and the maximum estimated SuSiE-RSS consistency parameter across diagnosed loci.

Reference implementation used for the diagnostic logic: xinhe-lab/mapgen, commit `393e66fe62442d2c499e4cab10a7f48309da216a`, `R/LD_diagnosis.R`.
