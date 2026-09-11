from pathlib import Path

CORE = Path("R/CATER_MR_core_v05.R")
TEST = Path("tests/v06_hardening_smoke.R")
WF = Path(".github/workflows/fix-codex-review-once.yml")
SELF = Path("tools/fix_codex_review_p2.py")

core = CORE.read_text()

old_sig = '''                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,\n                     primary_policy=c("cis_anchor","screened_cater"),\n                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE) {\n  primary_policy<-match.arg(primary_policy)'''
new_sig = '''                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,\n                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,\n                     primary_policy=c("cis_anchor","screened_cater")) {\n  primary_policy<-match.arg(primary_policy)'''
if old_sig not in core:
    raise SystemExit("Could not find cater_mr public signature")
core = core.replace(old_sig, new_sig, 1)

old_decision = '''  if(primary_policy=="cis_anchor") {\n    use_fit("cis",fits$cis,"OK_CIS_ANCHOR_PRIMARY","CATER_SENSITIVITY_UNCALIBRATED")\n  } else if(identical(fits$combined$status,"OK")) {\n    use_fit("combined",fits$combined,"OK_CATER","NO_DETECTED_SIBLING_COPERTURBATION")\n  }\n  out\n}'''
new_decision = '''  if(primary_policy=="cis_anchor") {\n    use_fit("cis",fits$cis,"OK_CIS_ANCHOR_PRIMARY","CATER_SENSITIVITY_UNCALIBRATED")\n  } else if(identical(fits$combined$status,"OK")) {\n    use_fit("combined",fits$combined,"OK_CATER","NO_DETECTED_SIBLING_COPERTURBATION")\n  } else {\n    # screened_cater is an opt-in promotion policy, not permission to discard a\n    # valid conservative cis estimate when the augmented cis+trans fit fails.\n    use_fit("cis",fits$cis,"CATER_COMBINED_FAILED_CIS_FALLBACK","COMBINED_FIT_FAILED")\n  }\n  out\n}'''
if old_decision not in core:
    raise SystemExit("Could not find primary decision tail")
core = core.replace(old_decision, new_decision, 1)
CORE.write_text(core)

test = TEST.read_text()
needle = '''stopifnot("primary_policy" %in% names(formals(cater_mr)))\n'''
addition = '''stopifnot("primary_policy" %in% names(formals(cater_mr)))\n\n# 12. Backward compatibility: new options must not shift established trailing positional arguments.\nformal_names <- names(formals(cater_mr))\nstopifnot(identical(tail(formal_names, 4L), c("outdir","drop_palindromic","verbose","primary_policy")))\n\n# 13. screened_cater promotes a valid combined fit, but a failed augmented fit must fall back to cis.\nfits_combined_fail <- list(cis=mkfit("OK",b=.11,se=.04,p=.01),\n                           trans=mkfit("OK"),combined=mkfit("LD_SINGULAR"))\nd_combined_fail <- .cater_primary_decision(\n  fits_combined_fail,has_trans=TRUE,sibling_screen_performed=TRUE,\n  sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE,\n  primary_policy="screened_cater")\nstopifnot(d_combined_fail$model=="cis",\n          d_combined_fail$status=="CATER_COMBINED_FAILED_CIS_FALLBACK",\n          d_combined_fail$evidence_status=="COMBINED_FIT_FAILED")\n'''
if needle not in test:
    raise SystemExit("Could not find hardening test insertion point")
test = test.replace(needle, addition, 1)
TEST.write_text(test)

# Remove the one-shot migration machinery from the final branch state.
SELF.unlink()
WF.unlink()
print("Applied Codex review fixes and removed one-shot patch machinery")
