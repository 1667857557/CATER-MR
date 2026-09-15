from pathlib import Path

core = Path("R/CATER_MR_core_v05.R")
s = core.read_text()

# 1) MVMR must create the per-target directory before PLINK writes <prefix>.extract.
old = 'umap<-data.frame(snp=u$snp,source="mvmr",parent_tf="",locus_id="mvmr",provenance="MVMR_UNION",stringsAsFactors=FALSE);ldsel<-tryCatch(.cater_ld_select(u,umap,ld_bfile,plink_bin,ld_threads,file.path(outdir,"mvmr",target,"ld"),ld_clump_r2,ld_clump_kb,verbose),error=function(e)e);'
new = 'umap<-data.frame(snp=u$snp,source="mvmr",parent_tf="",locus_id="mvmr",provenance="MVMR_UNION",stringsAsFactors=FALSE);mvmr_dir<-file.path(outdir,"mvmr",target);dir.create(mvmr_dir,recursive=TRUE,showWarnings=FALSE);ldsel<-tryCatch(.cater_ld_select(u,umap,ld_bfile,plink_bin,ld_threads,file.path(mvmr_dir,"ld"),ld_clump_r2,ld_clump_kb,verbose),error=function(e)e);'
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit("MVMR directory patch target not found")

# 2) Centralize trans ascertainment/data-contract semantics.
needle = '.cater_validate_input_manifest <- function(input_manifest,strict=TRUE){legacy_keys<-c("grn_build","eqtl_build","ld_build","outcome_build","eqtl_n_unit","eqtl_full_summary","eqtl_ancestry","outcome_ancestry","ld_ancestry","trans_qtl_qc","outcome_ld_ancestry");if(is.null(input_manifest)){if(isTRUE(strict)).cater_stop("strict_input_contract=TRUE requires input_manifest");z<-as.list(setNames(rep(NA_character_,length(legacy_keys)),legacy_keys));z$eqtl_full_summary<-NA;return(z)};if(!is.list(input_manifest)).cater_stop("input_manifest must be a named list");z<-input_manifest;norm_build<-function(x)toupper(gsub("[^A-Z0-9]","",as.character(x)));bs<-vapply(z[intersect(c("grn_build","eqtl_build","ld_build"),names(z))],norm_build,character(1));bs<-bs[nzchar(bs)&!is.na(bs)];if(length(bs)>1L&&length(unique(bs))>1L).cater_stop("GRN/eQTL/LD genome builds must match");if(isTRUE(strict)&&!identical(tolower(as.character(z$eqtl_n_unit%||%"")),"donors")).cater_stop("eQTL N must denote genetically independent donors, not cells");new_contract<-all(c("cis_full_summary","trans_data_mode")%in%names(z));if(isTRUE(strict)){if(new_contract){if(!isTRUE(z$cis_full_summary)).cater_stop("cis_full_summary must be TRUE");if(!tolower(as.character(z$trans_data_mode))%in%c("significant_only","full_summary")).cater_stop("trans_data_mode must be significant_only or full_summary")}else if(!isTRUE(z$eqtl_full_summary)).cater_stop("Legacy strict manifest requires eqtl_full_summary=TRUE; use cis_full_summary=TRUE and trans_data_mode=\'significant_only\' for the v0.8 contract")};if(isTRUE(strict)&&all(c("eqtl_ancestry","ld_ancestry")%in%names(z))&&!identical(tolower(as.character(z$eqtl_ancestry)),tolower(as.character(z$ld_ancestry)))).cater_stop("LD ancestry must match eQTL ancestry");z}\n'
helper = needle + '\n.cater_resolve_trans_contract <- function(manifest,trans_eqtl,instrument_p,trans_reporting_p){\n  raw_mode<-manifest$trans_data_mode\n  mode<-if(is.null(raw_mode)||length(raw_mode)!=1L||is.na(raw_mode)||!nzchar(as.character(raw_mode))) "" else tolower(as.character(raw_mode))\n  if(identical(mode,"significant_only")&&is.null(trans_eqtl)) .cater_stop("trans_data_mode=\'significant_only\' requires explicit trans_eqtl; refusing legacy full-summary fallback")\n  effective_mode<-if(nzchar(mode)) mode else if(!is.null(trans_eqtl)) "significant_only" else "legacy_full_summary"\n  if(!is.null(trans_eqtl)&&!identical(effective_mode,"full_summary")){\n    if(length(trans_reporting_p)!=1L||!is.finite(trans_reporting_p)||trans_reporting_p<=0||trans_reporting_p>1) .cater_stop("trans_reporting_p must be in (0,1] for significant-only trans input")\n    if(trans_reporting_p<instrument_p) .cater_stop("trans_reporting_p (%g) is more stringent than instrument_p (%g); eligible trans IVs between the two thresholds are unobserved",trans_reporting_p,instrument_p)\n    if(trans_reporting_p>instrument_p) warning("trans_reporting_p is less stringent than instrument_p; CATER-MR will apply instrument_p to reported hits")\n  }\n  effective_mode\n}\n'
if '.cater_resolve_trans_contract <- function(' not in s:
    if needle not in s:
        raise SystemExit("manifest helper insertion target not found")
    s = s.replace(needle, helper, 1)

# 3) Expose explicit network-primary gates without changing conservative defaults.
old = '                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,\n                     input_manifest=NULL,strict_input_contract=FALSE) {'
new = '                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,\n                     allow_network_primary=FALSE,sibling_screen_independent=FALSE,\n                     input_manifest=NULL,strict_input_contract=FALSE) {'
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit("public API gate patch target not found")

old = '  if(length(mvmr_max_condition)!=1L||!is.finite(mvmr_max_condition)||mvmr_max_condition<=1).cater_stop("mvmr_max_condition must be >1")\n  manifest<-.cater_validate_input_manifest(input_manifest,strict_input_contract)\n'
new = '  if(length(mvmr_max_condition)!=1L||!is.finite(mvmr_max_condition)||mvmr_max_condition<=1).cater_stop("mvmr_max_condition must be >1")\n  if(!is.logical(allow_network_primary)||length(allow_network_primary)!=1L||is.na(allow_network_primary)).cater_stop("allow_network_primary must be TRUE or FALSE")\n  if(!is.logical(sibling_screen_independent)||length(sibling_screen_independent)!=1L||is.na(sibling_screen_independent)).cater_stop("sibling_screen_independent must be TRUE or FALSE")\n  manifest<-.cater_validate_input_manifest(input_manifest,strict_input_contract)\n  trans_data_mode<-.cater_resolve_trans_contract(manifest,trans_eqtl,instrument_p,trans_reporting_p)\n'
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit("contract resolver call patch target not found")

# 4) Remove the old warning-only threshold branch now covered by the resolver.
old = '  outcome<-.cater_standardize_sumstats(outcome,label="outcome",require_position=FALSE);trans_hits<-.cater_standardize_trans_hits(trans_eqtl,trans_gene_col,qtl_n)\n  if(!is.null(trans_hits)&&nrow(trans_hits)&&is.finite(trans_reporting_p)){if(trans_reporting_p>instrument_p)warning("trans_reporting_p is less stringent than instrument_p; CATER-MR will apply instrument_p to reported hits");if(trans_reporting_p<instrument_p)warning("trans_reporting_p is more stringent than instrument_p; the trans instrument pool is availability-censored at the source reporting threshold")}\n'
new = '  outcome<-.cater_standardize_sumstats(outcome,label="outcome",require_position=FALSE);trans_hits<-.cater_standardize_trans_hits(trans_eqtl,trans_gene_col,qtl_n)\n'
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit("old trans threshold warning block not found")

# 5) Pass explicit gates into the primary decision.
old = 'dec<-.cater_primary_decision(fits,net,target,has_trans,srow$sibling_screen_performed,if(has_trans)isTRUE(srow$sibling_screen_complete)else TRUE,length(sib$active),primary_policy,sibling_screen_testable=if(has_trans)isTRUE(sib$testable%||%FALSE)else TRUE);'
new = 'dec<-.cater_primary_decision(fits,net,target,has_trans,srow$sibling_screen_performed,if(has_trans)isTRUE(srow$sibling_screen_complete)else TRUE,length(sib$active),primary_policy,allow_network_primary=allow_network_primary,sibling_screen_independent=sibling_screen_independent,sibling_screen_testable=if(has_trans)isTRUE(sib$testable%||%FALSE)else TRUE);'
if old in s:
    s = s.replace(old, new, 1)
elif new not in s:
    raise SystemExit("primary decision gate patch target not found")

core.write_text(s)

# Regression tests for the Codex review findings.
test = Path("tests/v08_significant_trans_smoke.R")
t = test.read_text()
marker = '# Codex review regression tests (PR #13)'
if marker not in t:
    t += '''\n\n# Codex review regression tests (PR #13)\n# Significant-only manifests must never silently fall back to legacy full-summary trans.\nm_sig <- list(trans_data_mode="significant_only")\nstopifnot(inherits(try(.cater_resolve_trans_contract(m_sig,NULL,5e-8,5e-8),silent=TRUE),"try-error"))\nstopifnot(identical(.cater_resolve_trans_contract(m_sig,data.frame(),5e-8,5e-8),"significant_only"))\n\n# A source reporting threshold stricter than the requested instrument threshold is not identifiable.\nstopifnot(inherits(try(.cater_resolve_trans_contract(m_sig,data.frame(),5e-8,1e-10),silent=TRUE),"try-error"))\n\n# Network primary requires explicit policy + conditional-F acceptance + both public gates.\nmk_review_fit <- function(status="OK",b=.2,se=.05,p=.01) data.frame(n_iv=2L,beta=b,se=se,p=p,Q=NA,Q_p=NA,information=10,joint_wald_per_df=5,effective_F=5,precision_information=10,mean_F=20,min_F=15,status=status)\nfits_review <- list(cis=mk_review_fit(),trans=mk_review_fit(),combined=mk_review_fit())\nnet_review <- list(status="OK",primary_eligible=TRUE,beta=c(X=.3),se=c(X=.06),p=c(X=.001))\nd_blocked <- .cater_primary_decision(fits_review,net_review,"X",TRUE,TRUE,TRUE,1L,"screened_cater",allow_network_primary=FALSE,sibling_screen_independent=TRUE,sibling_screen_testable=TRUE)\nd_allowed <- .cater_primary_decision(fits_review,net_review,"X",TRUE,TRUE,TRUE,1L,"screened_cater",allow_network_primary=TRUE,sibling_screen_independent=TRUE,sibling_screen_testable=TRUE)\nstopifnot(d_blocked$model=="cis",d_allowed$model=="network")\n\n# The MVMR builder creates a per-target directory before calling LD selection.\nb_mvmr <- paste(deparse(body(.cater_build_mvmr)),collapse="\\n")\nstopifnot(grepl("dir.create\\\\(mvmr_dir",b_mvmr),grepl("file.path\\\\(mvmr_dir, \\\"ld\\\"\\\\)",b_mvmr))\n\n# Candidate frames are normalized to a common summary-stat schema before rbind.\nstopifnot(grepl("d <- d\\\\[, core, drop = FALSE\\\\]",gsub(";", "; ", b_mvmr),fixed=FALSE) || grepl("d<-d\\\\[,core,drop=FALSE\\\\]",b_mvmr))\n'''
    test.write_text(t)
