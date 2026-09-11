# Pre-COJO LD consistency diagnosis

CATER-MR performs an LD/summary-statistic consistency gate before every Manc-COJO selection step by default.

Candidate loci are defined from the target cis window and all parent-TF windows together. Overlapping intervals on the same chromosome are merged into connected physical components before candidate SNP assignment. Any component containing the target cis window is treated in full as the cis locus (`locus_id=cis`), including TF-window extensions connected to it; those SNPs are not reused as trans instruments. Components containing only TF windows remain local trans loci such as `TF:A+B`.

## Method

For each physical local cis/TF locus, CATER-MR:

1. merges the target cis window and parent-TF windows by genomic interval into connected physical loci; any component containing the target cis window is the full cis locus, while TF-only components remain trans loci;
2. assigns candidate SNPs to these physical loci; `parent_tf` is retained only for SNPs in TF-only trans components, whereas SNPs absorbed into the cis component have no trans parent assignment;
3. extracts candidate variants for that locus from the configured PLINK LD reference;
4. computes a signed allele-count correlation matrix with PLINK 1.9 (`--r square`);
5. aligns eQTL z-score signs to direct PLINK A1/A2 same/swap coding for non-palindromic variants; strand-complement relationships are recorded for audit but excluded because the downstream CATER-MR COJO/MR allele contract also requires direct A1/A2 agreement;
6. estimates the SuSiE-RSS consistency parameter with `susieR::estimate_s_rss()` and obtains conditional z-score diagnostics with `susieR::kriging_rss()`;
7. removes variants satisfying `logLR > 2 & abs(z) > 2` before Manc-COJO runs.

`locus_id` is the physical cis/TF locus used for LD diagnosis and is also propagated with the retained instruments. TF-only components retain `parent_tf` for biological attribution. If a TF component overlaps the target cis component, the entire connected component is classified as cis and its SNPs are deliberately not assigned a trans parent, avoiding cis/trans double interpretation within one local LD neighbourhood.

Variants absent from the LD reference, variants with non-finite LD rows, and variants whose allele pair cannot be reconciled with the direct A1/A2 contract are also removed and explicitly reported. All variants whose original PLINK LD row contains any non-finite entry are removed simultaneously, so this QC is not variant-order dependent. Palindromic A/T and C/G variants are not assigned a signed z score from allele labels alone; without frequency-based strand disambiguation they are excluded from SuSiE-RSS diagnosis, retained for downstream COJO, and reported as `NOT_DIAGNOSABLE_PALINDROMIC`. A one-variant diagnosable locus is retained with status `NOT_DIAGNOSABLE_SINGLETON`.

The PLINK matrix is checked for dimensions, finite entries, symmetry, unit diagonal, and correlation range. CATER-MR deliberately does not impose an additional strict positive-semidefinite eigenvalue gate before SuSiE-RSS, because `susieR` performs the RSS eigenvalue handling itself and text-formatted PLINK matrices can contain small rounding artifacts.

PLINK subset generation is fail-closed: a non-zero PLINK exit status is treated as an error except for the explicit `No variants remaining` case, which is reported as LD-reference missing.

This is a QC gate only. It does not use SuSiE PIP or credible sets for instrument selection; Manc-COJO remains the instrument-selection method.

## Dependencies and options

With the default `enable_ld_diagnosis=TRUE`, PLINK 1.9 and the R package `susieR` are required. `plink_bin` can be an executable name on `PATH` or an explicit path. The default detection thresholds are `ld_diag_loglr=2` and `ld_diag_abs_z=2`. `enable_ld_diagnosis` must be a scalar logical (`TRUE` or `FALSE`); numeric or string truthy values are rejected rather than silently disabling the QC gate.

Set `enable_ld_diagnosis=FALSE` only to reproduce the legacy COJO path without this QC gate.

## Output

For target `GENE`, the full SNP-level diagnostic table is written to:

`<outdir>/cojo/GENE.ld_diagnosis.tsv`

The target summary also reports the number removed by each major reason and the maximum estimated SuSiE-RSS consistency parameter across diagnosed loci.

The Codex review hardening is covered by smoke tests for unresolved palindromic variants, simultaneous removal of every non-finite LD row, strict logical validation of `enable_ld_diagnosis`, and physical merging of overlapping TF windows.

Reference implementation used for the diagnostic logic: xinhe-lab/mapgen, commit `393e66fe62442d2c499e4cab10a7f48309da216a`, `R/LD_diagnosis.R`.
