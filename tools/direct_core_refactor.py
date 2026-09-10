from pathlib import Path
import re

CORE = Path("R/CATER_MR_core_v05.R")
OVERRIDE = Path("R/CATER_MR_v06_overrides.R")
ENTRY = Path("CATER_MR.R")
README = Path("README.md")
SMOKE = Path("tests/smoke.R")
HARDEN = Path("tests/v06_hardening_smoke.R")


def scan_function_blocks(text):
    pat = re.compile(r"(?m)^((?:`[^`]+`)|(?:[A-Za-z.][A-Za-z0-9._]*))\s*<-\s*function\s*\(")
    out = {}
    for m in pat.finditer(text):
        name = m.group(1)
        start = m.start()
        brace = text.find("{", m.end())
        if brace < 0:
            continue
        depth = 0
        quote = None
        escape = False
        comment = False
        end = None
        i = brace
        while i < len(text):
            ch = text[i]
            if comment:
                if ch == "\n":
                    comment = False
                i += 1
                continue
            if quote is not None:
                if escape:
                    escape = False
                elif ch == "\\":
                    escape = True
                elif ch == quote:
                    quote = None
                i += 1
                continue
            if ch in ('"', "'"):
                quote = ch
                i += 1
                continue
            if ch == "#":
                comment = True
                i += 1
                continue
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    while end < len(text) and text[end] in " \t":
                        end += 1
                    if end < len(text) and text[end] == "\n":
                        end += 1
                    break
            i += 1
        if end is not None:
            out[name] = (start, end, text[start:end])
    return out


def replace_function(text, name, block):
    blocks = scan_function_blocks(text)
    if name not in blocks:
        marker = scan_function_blocks(text).get("cater_mr")
        if marker is None:
            return text.rstrip() + "\n\n" + block.rstrip() + "\n"
        pos = marker[0]
        return text[:pos] + block.rstrip() + "\n\n" + text[pos:]
    s, e, _ = blocks[name]
    return text[:s] + block.rstrip() + "\n\n" + text[e:]


core = CORE.read_text()
over = OVERRIDE.read_text()
oblocks = scan_function_blocks(over)

# Consolidate the v0.6 hardening functions into the actual core.  Conditional-F
# is kept from the core because the override implementation is intentionally a
# wrapper around the old function; below we add the method label directly.
skip = {".cater_conditional_f"}
for name, (_, _, block) in oblocks.items():
    if name in skip:
        continue
    core = replace_function(core, name, block)

# `%||%` is a one-line function and is not captured by the brace-based parser.
if "`%||%` <- function" not in core:
    marker = scan_function_blocks(core)["cater_mr"][0]
    core = core[:marker] + '`%||%` <- function(x,y) if(is.null(x)) y else x\n\n' + core[marker:]

# Direct version marker/header.
core = re.sub(
    r"^# CATER-MR v0\.5\n# Cis And Trans eQTLs guided by Regulatory networks for drug-target MR\n# Target-centric Manc-COJO \+ GIVW \+ complete same-IV pleiotropy screen \+ triggered local MVMR\.\n",
    "# CATER-MR v0.7\n# Cis And Trans eQTLs guided by Regulatory networks for drug-target MR\n# Direct implementation: GRN-gated cis/trans selection + standard LD-aware MR.\n\n.CATER_VERSION <- \"0.7.0\"\n",
    core,
    count=1,
)

# The input manifest no longer carries qtl_outcome_overlap.  This field was
# provenance only and was never an estimator-level correction.
core = core.replace(',"outcome_ld_ancestry","qtl_outcome_overlap")', ',"outcome_ld_ancestry")')
core = core.replace(', "outcome_ld_ancestry", "qtl_outcome_overlap")', ', "outcome_ld_ancestry")')

# Add the conditional-F method label directly rather than through an override.
cf = scan_function_blocks(core).get(".cater_conditional_f")
if cf is not None:
    s, e, block = cf
    if 'attr(ans,"method")' not in block and 'attr(ans, "method")' not in block:
        block = re.sub(r"\n\s*ans\n}\s*$", '\n  attr(ans,"method") <- "CORRELATED_IV_CONDITIONAL_F"\n  ans\n}\n', block)
        core = core[:s] + block + core[e:]

# Reframe the same-IV sibling test as co-perturbation detection.  `active` is
# retained as an internal compatibility alias for MVMR triggering.
old = '''  tab$active <- is.finite(tab$q)&tab$q<sibling_fdr
  list(table=tab,active=tab$sibling[tab$active],n_candidate=nrow(tab),n_incomplete=sum(!tab$complete),n_qtl_missing=qtl_missing,n_iv_missing=sum(tab$n_iv_missing),complete=all(tab$complete),testable=any(tab$testable),status=if(all(tab$complete))"COMPLETE" else "INCOMPLETE")
'''
new = '''  tab$co_perturbation_detected<-is.finite(tab$q)&tab$q<sibling_fdr
  tab$interpretation<-ifelse(tab$co_perturbation_detected,"DETECTED_SIBLING_COPERTURBATION",ifelse(tab$testable,"NO_DETECTED_SIBLING_COPERTURBATION","NOT_TESTABLE"))
  tab$active<-tab$co_perturbation_detected
  active<-tab$sibling[tab$active]
  complete<-all(tab$complete)
  list(table=tab,active=active,co_perturbed=active,n_candidate=nrow(tab),n_incomplete=sum(!tab$complete),n_qtl_missing=qtl_missing,n_iv_missing=sum(tab$n_iv_missing),complete=complete,testable=any(tab$testable),status=if(!complete)"INCOMPLETE" else if(length(active))"COMPLETE_DETECTED_COPERTURBATION" else "COMPLETE_NO_DETECTED_COPERTURBATION")
'''
if old not in core:
    raise SystemExit("Could not find sibling-screen tail to replace")
core = core.replace(old, new, 1)

# Default public/internal primary interpretation is the conventional cis anchor.
core = core.replace('primary_policy=c("screened_cater","cis_anchor")', 'primary_policy=c("cis_anchor","screened_cater")')
core = core.replace('"MEASURED_SIBLING_SCREEN_NEGATIVE"', '"NO_DETECTED_SIBLING_COPERTURBATION"')
core = core.replace('"MEASURED_CO_PERTURBATION"', '"DETECTED_SIBLING_COPERTURBATION"')

# Expose the already-existing policy at the public entrypoint and keep the
# default simple: cis is primary; screened CATER is an explicit research opt-in.
sig_old = 'exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,\n                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE) {'
sig_new = 'exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,\n                     primary_policy=c("cis_anchor","screened_cater"),\n                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE) {\n  primary_policy<-match.arg(primary_policy)'
if sig_old not in core:
    raise SystemExit("Could not find cater_mr signature")
core = core.replace(sig_old, sig_new, 1)

call_old = '      sibling_screen_complete=if(has_trans) isTRUE(srow$sibling_screen_complete) else TRUE,\n      n_active_siblings=length(sib$active))'
call_new = '      sibling_screen_complete=if(has_trans) isTRUE(srow$sibling_screen_complete) else TRUE,\n      n_active_siblings=length(sib$active),primary_policy=primary_policy,\n      sibling_screen_testable=if(has_trans) isTRUE(sib$testable %||% FALSE) else TRUE)'
if call_old not in core:
    raise SystemExit("Could not find primary decision call")
core = core.replace(call_old, call_new, 1)

# Add a transparent TF-anchor evidence label; it is annotation only and never a
# hard IV gate.
anchor_old = '''      if(length(aa)) anchor$q_tf[aa]<-p.adjust(anchor$p_tf[aa],method="BH")
      utils::write.table(anchor,file.path(outdir,"mechanism",paste0(target,"_tf_anchor.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
'''
anchor_new = '''      if(length(aa)) anchor$q_tf[aa]<-p.adjust(anchor$p_tf[aa],method="BH")
      anchor$evidence<-ifelse(anchor$status!="OK","TF_ANCHOR_NOT_TESTABLE",ifelse(is.finite(anchor$q_tf)&anchor$q_tf<0.05,"TF_ANCHOR_SUPPORTED","TF_ANCHOR_NOT_SUPPORTED"))
      utils::write.table(anchor,file.path(outdir,"mechanism",paste0(target,"_tf_anchor.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
'''
if anchor_old not in core:
    raise SystemExit("Could not find TF anchor output block")
core = core.replace(anchor_old, anchor_new, 1)

# Network MVMR remains standard GLS and is reported as sensitivity output only.
network_old = 'information=NA,effective_F=net$conditional_F[target],mean_F=NA,min_F=NA,status=if(net$primary_eligible&&isTRUE(sib$complete))"OK" else "SENSITIVITY_ONLY",'
network_new = 'information=NA,joint_wald_per_df=NA,effective_F=net$conditional_F[target],precision_information=NA,mean_F=NA,min_F=NA,status="SENSITIVITY_ONLY",'
if network_old not in core:
    raise SystemExit("Could not find network long-row block")
core = core.replace(network_old, network_new, 1)

# Final guard: the direct core must not contain the deleted field or duplicate
# top-level function definitions.
if "qtl_outcome_overlap" in core:
    raise SystemExit("qtl_outcome_overlap remains in direct core")
names = list(scan_function_blocks(core))
if len(names) != len(set(names)):
    raise SystemExit("Duplicate top-level function definitions remain")

CORE.write_text(core)
OVERRIDE.unlink()
ENTRY.write_text('# CATER-MR stable entrypoint\n# Direct implementation: functions are defined once in the core; no runtime overrides.\n\nsource(file.path("R", "CATER_MR_core_v05.R"), local = FALSE)\n')

# Documentation: keep wording conservative and backward compatible with output
# paths, but remove claims that a negative screen proves exclusion restriction.
readme = README.read_text()
readme = readme.replace('complete same-IV one-hop sibling pleiotropy screen', 'complete same-IV one-hop sibling co-perturbation screen')
readme = readme.replace('local MVMR only when measured bypass pleiotropy is detected', 'local MVMR sensitivity analysis only when measured sibling co-perturbation is detected')
readme = readme.replace('## Same-IV one-hop pleiotropy screen', '## Same-IV one-hop sibling co-perturbation screen')
readme = readme.replace('A significant partial test can still identify a bypass and trigger MVMR, but partial or failed coverage can never be interpreted as evidence that the remaining pathway is clean.', 'A significant partial test can still identify measured sibling co-perturbation and trigger MVMR, but a non-significant test is only "no detected co-perturbation"; it is not proof of the MR exclusion restriction. Partial or failed coverage can never be interpreted as evidence that the remaining pathway is clean.')
readme = readme.replace('MVMR is not run by default.\n\nIt is triggered only if the exposure-side sibling screen detects at least one measured bypass gene.', 'MVMR is triggered only if the exposure-side sibling screen detects at least one measured co-perturbed sibling. It remains a standard LD-aware GLS-MVMR sensitivity analysis and is not interpreted as a complete correction for unmeasured pleiotropy, weak-instrument bias, or sample overlap.')
readme = readme.replace('A network estimate becomes **primary-eligible** only when:', 'For the MVMR sensitivity estimate, CATER-MR reports the following identification diagnostics:')
readme = readme.replace('Without exposure covariance, the network coefficient may still be reported as sensitivity output, but CATER-MR does not promote it to the primary result.', 'The network coefficient is reported as sensitivity output. `exposure_corr` is used for conditional-strength assessment; it does not turn the standard GLS point estimator into a weak-instrument or sample-overlap bias-corrected estimator.')
readme = readme.replace('- trans IVs with a **complete** sibling screen and no detected bypass: the combined CATER estimate can be primary;\n- detected bypass: use an identifiable/strong local network MVMR; otherwise fall back to cis;', '- default `primary_policy="cis_anchor"`: use cis MR as the primary anchor whenever a valid cis estimate is available;\n- `primary_policy="screened_cater"` is an explicit research opt-in: a complete screen with no detected sibling co-perturbation can promote the combined estimate;\n- detected sibling co-perturbation: keep cis as primary when available and report local MVMR as sensitivity;')
README.write_text(readme)

# Add regression checks for the new conservative semantics.
smoke = SMOKE.read_text()
needle = 'stopifnot(sib$table$n_iv_requested==2L,sib$table$n_iv_tested==1L,\n          sib$table$status=="PARTIAL_SNP_COVERAGE",isTRUE(sib$table$active))\n'
if needle in smoke and 'co_perturbation_detected' not in smoke:
    smoke = smoke.replace(needle, needle + 'stopifnot("co_perturbation_detected" %in% names(sib$table),\n          sib$table$interpretation=="DETECTED_SIBLING_COPERTURBATION")\n', 1)
SMOKE.write_text(smoke)

hard = HARDEN.read_text()
append = '''\n# 11. Direct-core Occam semantics: no qtl-outcome-overlap field and cis is the default primary anchor.\nmanifest_occam <- .cater_validate_input_manifest(good, TRUE)\nstopifnot(!"qtl_outcome_overlap" %in% names(manifest_occam))\nfits_occam <- list(cis=mkfit("OK",b=.11,se=.04,p=.01),trans=mkfit("OK"),combined=mkfit("OK",b=.20,se=.03,p=.001))\nd_occam <- .cater_primary_decision(fits_occam,has_trans=TRUE,sibling_screen_performed=TRUE,\n                                   sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE)\nstopifnot(d_occam$model=="cis",d_occam$status=="OK_CIS_ANCHOR_PRIMARY")\nstopifnot("primary_policy" %in% names(formals(cater_mr)))\n'''
if 'manifest_occam <-' not in hard:
    hard = hard.replace('cat("CATER-MR v0.6 evidence-safe hardening tests passed\\n")', append + '\ncat("CATER-MR direct-core evidence-safety tests passed\\n")')
HARDEN.write_text(hard)

print("Direct-core refactor prepared")
