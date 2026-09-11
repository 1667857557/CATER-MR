from pathlib import Path
import re

p = Path('README.md')
s = p.read_text()


def sub_once(pattern, repl, text, label):
    out, n = re.subn(pattern, repl, text, count=1, flags=re.S)
    if n != 1:
        raise SystemExit(f'{label}: expected exactly one replacement, got {n}')
    return out

s = sub_once(
    r"## V0\.5 design\n.*?\n## Primary inputs",
    """## Current design (v0.7.3)

CATER-MR keeps the analysis target-centric and deliberately separates biological gating, genetic QC, independent-signal selection, and MR estimation.

```text
cell-type-specific GRN + full-summary eQTL + ancestry-matched LD
    -> default genetic QC
    -> target-cis + direct parent-TF physical loci
    -> SuSiE-RSS LD/summary consistency QC
    -> Manc-COJO independent-signal selection
    -> LD-aware cis / trans / combined MR
    -> trans mechanism and pleiotropy sensitivity analyses
```

The core principles are:

> **SNPs are the instruments; the GRN is a biological gate; SuSiE-RSS is a pre-COJO consistency QC; Manc-COJO is the independent-signal selector.**

GRN edge weights do not rescale eQTL effects. Outcome associations are not used to select exposure instruments.

## Primary inputs""",
    s,
    'design section'
)

s = sub_once(
    r"### LD reference\n\n`ld_bfile` is a PLINK `bed/bim/fam` prefix\. The preferred reference is the genotype sample used for the eQTL mapping; otherwise use a large ancestry-matched panel\.\n\n## Manc-COJO",
    """### LD reference

`ld_bfile` is a PLINK `bed/bim/fam` prefix. The preferred reference is the genotype sample used for the eQTL mapping; otherwise use a large ancestry-matched panel.

### Genome build and default genetic QC

The current coordinate contract is **GRCh38/hg38**. CATER-MR does not perform liftOver internally.

By default:

- palindromic A/T and C/G SNPs are removed before locus construction (`drop_palindromic=TRUE`);
- the extended MHC is excluded using hg38 `chr6:25,000,000-36,000,000` (`exclude_mhc=TRUE`);
- targets whose TSS lies inside that interval are skipped;
- parent TFs whose TSS lies inside that interval do not contribute trans loci.

The same variant QC is propagated to sibling-cis selection used by the local MVMR sensitivity path.

## Physical loci, LD consistency QC, and Manc-COJO""",
    s,
    'input/QC section'
)

s = sub_once(
    r"## Physical loci, LD consistency QC, and Manc-COJO\n.*?\n## Allele orientation",
    """## Physical loci, LD consistency QC, and Manc-COJO

For each target `X`, CATER-MR starts from the target cis window and the local windows of its direct parent TFs. Overlapping windows on the same chromosome are merged into connected physical components before SNP assignment.

If a connected component contains the target cis window, the entire component is classified as `cis`. Components containing only parent-TF windows remain local `trans` loci. This prevents one local LD neighborhood from being interpreted simultaneously as direct cis and TF-mediated trans evidence.

Before COJO, each physical locus is checked against the configured PLINK LD reference with the SuSiE-RSS consistency diagnostic (`estimate_s_rss()` + `kriging_rss()`). Variants with LD/summary-statistic inconsistency, missing LD-reference support, non-finite LD rows, or irreconcilable allele coding are removed before conditional selection.

CATER-MR then performs one target-level Manc-COJO selection across the retained candidate SNPs from the cis and GRN-constrained trans loci. When multiple signals are retained, the joint Manc-COJO step supplies the signed LD matrix used by the downstream LD-aware MR estimators.

### Why COJO remains the IV selector

The current framework needs **conditionally independent exposure-association signals** for MR. Manc-COJO directly serves that role. SuSiE is therefore used here as an LD/summary-consistency QC, not as a fine-mapping selector: CATER-MR does not use PIP or credible sets to define IVs.

Fine-mapping can be added as a sensitivity or signal-interpretation layer, but it is not part of the default instrument-selection path.

## Allele orientation""",
    s,
    'locus/COJO section'
)

s = sub_once(
    r"## Allele orientation\n.*?\n## MR estimators",
    """## Allele orientation

Palindromic SNPs are removed before the default locus/COJO path. For selected non-palindromic variants, CATER-MR orients the target eQTL effect to the Manc-COJO reference-panel A1/A2 coding, flips beta/EAF when required, and then harmonizes the outcome to the same allele orientation.

This preserves a consistent signed-LD convention for the MR layer.

## MR estimators""",
    s,
    'allele section'
)

s = sub_once(
    r"## Diagnostics\n\nFor each target CATER-MR reports:\n\n.*?\n\nOutcome-based quantities are never used to select trans instruments or sibling exposures\.",
    """## Diagnostics

For each target CATER-MR reports the main QC, selection, strength and sensitivity diagnostics, including:

- numbers of hg38 MHC SNPs, palindromic SNPs and parent TFs excluded by the default genetic QC;
- pre-COJO LD-diagnosis removals and SuSiE-RSS consistency summaries;
- numbers of COJO-selected cis/trans signals and final MR-usable IVs;
- `min_F`, `mean_F`, `effective_F`, and trans information fraction;
- TF-anchor, sibling co-perturbation, cis/trans heterogeneity and leave-one-TF-locus diagnostics;
- local-MVMR identification diagnostics when that sensitivity analysis is triggered.

Outcome-based quantities are never used to select trans instruments or sibling exposures.""",
    s,
    'diagnostics section'
)

s = s.replace(
    "NO_CANDIDATE_SNP\nNO_COJO_SIGNAL",
    "TARGET_IN_MHC_EXCLUDED\nNO_CANDIDATE_SNP\nNO_CANDIDATE_AFTER_LD_DIAGNOSIS\nLD_DIAGNOSIS_FAILED\nNO_COJO_SIGNAL",
    1
)

s = sub_once(
    r"## Design validation\n\n`tests/smoke\.R` checks mathematical invariants and regression failure modes without requiring Manc-COJO itself, including:\n\n.*?\n\n## Scope deliberately not added",
    """## Design validation

The CI runs the core estimator smoke tests plus `tests/v06_hardening_smoke.R`. Together they cover the mathematical invariants and the current evidence-safety/QC contracts, including physical-locus merging, allele orientation, hg38 MHC exclusion, default palindromic removal, SuSiE-RSS diagnostic handling, non-finite LD handling, GIVW reductions, sibling-screen completeness, and local-MVMR safeguards.

## Scope deliberately not added""",
    s,
    'validation section'
)

s = s.replace(
    "V0.5 still does not make the following default components:",
    "The current implementation deliberately does not make the following default components:",
    1
)
s = s.replace(
    "- fine-mapping/colocalization as a hard gate;",
    "- fine-mapping/PIP/credible-set selection or colocalization as a hard IV gate;",
    1
)

if '## V0.5 design' in s or 'V0.5 still does not' in s:
    raise SystemExit('stale V0.5 wording remains')
if 'SuSiE-RSS is a pre-COJO consistency QC' not in s:
    raise SystemExit('new architecture summary missing')
if 'chr6:25,000,000-36,000,000' not in s:
    raise SystemExit('hg38 MHC contract missing')
if 'Fine-mapping can be added as a sensitivity' not in s:
    raise SystemExit('fine-mapping role missing')

p.write_text(s)
