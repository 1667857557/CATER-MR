from pathlib import Path
import re

ROOT = Path('.')
core_path = ROOT / 'R' / 'CATER_MR_core.R'
lean_path = ROOT / 'R' / 'CATER_MR_lean.R'
entry_path = ROOT / 'CATER_MR.R'
workflow_path = ROOT / '.github' / 'workflows' / 'r-smoke.yml'

core = core_path.read_text()

def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f'{label}: expected exactly one match, found {n}')
    return text.replace(old, new, 1)

# Version and top-level description.
core = replace_once(core, '# CATER-MR v0.8\n', '# CATER-MR v0.9\n', 'version header')
core = replace_once(core, '.CATER_VERSION <- "0.8.2"', '.CATER_VERSION <- "0.9.0"', 'version constant')
core = core.replace('ignored by CATER-MR v0.8:', 'ignored by CATER-MR v0.9:')
core = core.replace("for the v0.8 contract", "for the v0.9 contract")

# Internal qualification helpers are inserted into the original core, not loaded as overrides.
helpers = r'''
.cater_preld_tf_anchor <- function(tf,snp,target_row,eqtl_dir,qtl_n) {
  d<-.cater_read_gene(tf,eqtl_dir,qtl_n)
  if(is.null(d)) return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=NA_real_,p_tf=NA_real_,status="QTL_NOT_AVAILABLE",stringsAsFactors=FALSE))
  d<-d[d$snp==snp,,drop=FALSE]
  if(!nrow(d)) return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=NA_real_,p_tf=NA_real_,status="SNP_NOT_AVAILABLE",stringsAsFactors=FALSE))
  same<-d$a1[1]==target_row$a1[1]&&d$a2[1]==target_row$a2[1]
  swap<-d$a1[1]==target_row$a2[1]&&d$a2[1]==target_row$a1[1]
  if(!(same||swap)) return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=d$se[1],p_tf=d$p[1],status="ALLELE_MISMATCH",stringsAsFactors=FALSE))
  beta<-if(swap)-d$beta[1] else d$beta[1]
  p<-d$p[1];if(!is.finite(p))p<-2*stats::pnorm(-abs(beta/d$se[1]))
  data.frame(tf=tf,beta_tf=beta,se_tf=d$se[1],p_tf=p,status="OK",stringsAsFactors=FALSE)
}

.cater_qualify_trans_candidates <- function(target,qtl,map,eqtl_dir,qtl_n,trans_hits=NULL,
                                             tf_anchor_p=5e-8,max_reported_trans_targets=1L,
                                             trans_set=c("core","extended")) {
  trans_set<-match.arg(trans_set)
  if(length(tf_anchor_p)!=1L||!is.finite(tf_anchor_p)||tf_anchor_p<=0||tf_anchor_p>=1).cater_stop("tf_anchor_p must be in (0,1)")
  if(length(max_reported_trans_targets)!=1L||is.na(max_reported_trans_targets)||max_reported_trans_targets<1||
     (is.finite(max_reported_trans_targets)&&max_reported_trans_targets!=as.integer(max_reported_trans_targets)))
    .cater_stop("max_reported_trans_targets must be a positive integer or Inf")
  empty_counts<-list(n_trans_raw=0L,n_trans_anchor_pass=0L,n_trans_unique_parent=0L,n_trans_core=0L,n_trans_extended=0L,n_parent_tf_core=0L)
  if(!nrow(qtl)||!nrow(map))return(list(qtl=qtl,map=map,table=data.frame(),anchors=data.frame(),counts=empty_counts))
  ti<-which(map$source=="trans")
  if(!length(ti))return(list(qtl=qtl,map=map,table=data.frame(),anchors=data.frame(),counts=empty_counts))
  rows<-vector("list",length(ti));anchor_rows<-list()
  for(ii in seq_along(ti)){
    mi<-ti[ii];s<-map$snp[mi];qi<-match(s,qtl$snp);tr<-qtl[qi,,drop=FALSE]
    parents<-.cater_parent_tokens(map$parent_tf[mi])[[1L]]
    aa<-if(length(parents))do.call(rbind,lapply(parents,function(tf).cater_preld_tf_anchor(tf,s,tr,eqtl_dir,qtl_n)))else data.frame()
    if(nrow(aa)){
      aa$target<-target;aa$snp<-s;anchor_rows[[length(anchor_rows)+1L]]<-aa
      pass<-aa$status=="OK"&is.finite(aa$p_tf)&aa$p_tf<tf_anchor_p;qtf<-unique(aa$tf[pass])
    }else qtf<-character()
    nqual<-length(qtf);nrep<-NA_integer_
    if(!is.null(trans_hits)&&nrow(trans_hits))nrep<-length(unique(trans_hits$gene[trans_hits$snp==s]))
    specificity<-is.infinite(max_reported_trans_targets)||is.na(nrep)||nrep<=max_reported_trans_targets
    extended<-nqual==1L;core_ok<-extended&&specificity
    reason<-if(nqual==0L)"TRANS_NO_TF_CIS_ANCHOR"else if(nqual>1L)"TRANS_AMBIGUOUS_PARENT"else if(!specificity)"TRANS_REPORTED_HOTSPOT"else"TRANS_CORE"
    chosen<-if(nqual==1L)qtf[[1L]]else NA_character_;aone<-if(nqual==1L)aa[match(chosen,aa$tf),,drop=FALSE]else NULL
    rows[[ii]]<-data.frame(target=target,snp=s,raw_parent_tf=map$parent_tf[mi],qualified_parent=chosen,n_qualifying_parent=nqual,
      beta_tf=if(is.null(aone))NA_real_ else aone$beta_tf[1],se_tf=if(is.null(aone))NA_real_ else aone$se_tf[1],p_tf=if(is.null(aone))NA_real_ else aone$p_tf[1],
      n_reported_trans_targets=nrep,anchor_pass=nqual>=1L,unique_parent=extended,specificity_pass=specificity,
      eligible_extended=extended,eligible_core=core_ok,reason=reason,stringsAsFactors=FALSE)
  }
  tab<-do.call(rbind,rows);rownames(tab)<-NULL;anchors<-if(length(anchor_rows))do.call(rbind,anchor_rows)else data.frame()
  use<-if(trans_set=="core")tab$eligible_core else tab$eligible_extended;keep_snps<-tab$snp[use]
  cis_idx<-which(map$source=="cis");trans_idx<-ti[match(keep_snps,map$snp[ti],nomatch=0L)];trans_idx<-trans_idx[trans_idx>0L]
  keep_idx<-c(cis_idx,trans_idx);qout<-qtl[match(map$snp[keep_idx],qtl$snp),,drop=FALSE];mout<-map[keep_idx,,drop=FALSE]
  if(length(trans_idx)){
    dd<-tab[match(mout$snp[mout$source=="trans"],tab$snp),,drop=FALSE];mout$parent_tf[mout$source=="trans"]<-dd$qualified_parent
    prov<-if(trans_set=="core")"REPORTED_TRANS_CORE"else"REPORTED_TRANS_EXTENDED";mout$provenance[mout$source=="trans"]<-prov
    qout$provenance[qout$snp%in%mout$snp[mout$source=="trans"]]<-prov
  }
  counts<-list(n_trans_raw=nrow(tab),n_trans_anchor_pass=sum(tab$anchor_pass),n_trans_unique_parent=sum(tab$unique_parent),
    n_trans_core=sum(tab$eligible_core),n_trans_extended=sum(tab$eligible_extended),
    n_parent_tf_core=length(unique(tab$qualified_parent[tab$eligible_core&!is.na(tab$qualified_parent)])))
  list(qtl=qout,map=mout,table=tab,anchors=anchors,counts=counts)
}
'''.strip('\n')
marker = '.cater_target_candidates <- function'
if helpers.splitlines()[0] not in core:
    core = replace_once(core, marker, helpers + '\n\n' + marker, 'insert qualification helpers')

old_sig = '.cater_target_candidates <- function(target,parents,annotation,eqtl_dir,qtl_n,trans_hits,instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc) {'
new_sig = '.cater_target_candidates <- function(target,parents,annotation,eqtl_dir,qtl_n,trans_hits,instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc,tf_anchor_p=instrument_p,max_reported_trans_targets=1L,trans_set=c("core","extended")) {'
core = replace_once(core, old_sig, new_sig, 'target candidate signature')

old_tail = '''  if(!length(qs)) return(list(qtl=data.frame(),map=data.frame(),regions=regions))
  q<-do.call(rbind,qs);mp<-do.call(rbind,ms);dup<-duplicated(q$snp);q<-q[!dup,,drop=FALSE];mp<-mp[!dup,,drop=FALSE];rownames(q)<-NULL;rownames(mp)<-NULL
  list(qtl=q,map=mp,regions=regions)
}'''
new_tail = '''  empty_counts<-list(n_trans_raw=0L,n_trans_anchor_pass=0L,n_trans_unique_parent=0L,n_trans_core=0L,n_trans_extended=0L,n_parent_tf_core=0L)
  if(!length(qs)) return(list(qtl=data.frame(),map=data.frame(),regions=regions,qualification=data.frame(),anchors=data.frame(),qualification_counts=empty_counts))
  q<-do.call(rbind,qs);mp<-do.call(rbind,ms);dup<-duplicated(q$snp);q<-q[!dup,,drop=FALSE];mp<-mp[!dup,,drop=FALSE];rownames(q)<-NULL;rownames(mp)<-NULL
  qualified<-.cater_qualify_trans_candidates(target,q,mp,eqtl_dir,qtl_n,trans_hits,tf_anchor_p,max_reported_trans_targets,trans_set)
  list(qtl=qualified$qtl,map=qualified$map,regions=regions,qualification=qualified$table,anchors=qualified$anchors,qualification_counts=qualified$counts)
}'''
core = replace_once(core, old_tail, new_tail, 'target candidate qualification')

# Directly simplify the original primary-decision function.
start = core.index('.cater_primary_decision <- function')
end = core.index('\n\n.cater_validate_input_manifest', start)
primary = r'''.cater_primary_decision <- function(fits,net=NULL,target=NULL,has_trans=FALSE,sibling_screen_performed=FALSE,sibling_screen_complete=TRUE,n_active_siblings=0L,primary_policy=c("cis_anchor","screened_cater"),allow_trans_only_primary=FALSE,n_trans_tf_loci=0L,min_trans_tf_loci=3L,allow_network_primary=FALSE,sibling_screen_independent=FALSE,sibling_screen_testable=TRUE){
  primary_policy<-match.arg(primary_policy)
  out<-list(model=NA_character_,beta=NA_real_,se=NA_real_,p=NA_real_,status="NO_CIS_PRIMARY",evidence_status=if(has_trans)"TRANS_SENSITIVITY_ONLY"else"NO_CIS_INSTRUMENT")
  if(!is.null(fits$cis)&&identical(fits$cis$status,"OK")){
    out$model<-"cis";out$beta<-fits$cis$beta;out$se<-fits$cis$se;out$p<-fits$cis$p;out$status<-"OK_CIS_PRIMARY";out$evidence_status<-"CIS_ANCHOR_PRIMARY"
  }
  out
}'''
core = core[:start] + primary + core[end:]

# Public API: advanced sensitivity defaults off; Lean controls appended to preserve positional compatibility.
old_api = '''                     sibling_fdr=0.05,enable_sibling_screen=TRUE,enable_mvmr=TRUE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2,
                     exclude_mhc=TRUE,trans_eqtl=NULL,trans_gene_col=NULL,instrument_p=5e-8,
                     trans_reporting_p=5e-8,ld_clump_r2=.01,ld_clump_kb=10000L,ld_threads=1L,
                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,
                     allow_network_primary=FALSE,sibling_screen_independent=FALSE,
                     input_manifest=NULL,strict_input_contract=FALSE) {'''
new_api = '''                     sibling_fdr=0.05,enable_sibling_screen=FALSE,enable_mvmr=FALSE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2,
                     exclude_mhc=TRUE,trans_eqtl=NULL,trans_gene_col=NULL,instrument_p=5e-8,
                     trans_reporting_p=5e-8,ld_clump_r2=.01,ld_clump_kb=10000L,ld_threads=1L,
                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,
                     allow_network_primary=FALSE,sibling_screen_independent=FALSE,
                     input_manifest=NULL,strict_input_contract=FALSE,
                     tf_anchor_p=instrument_p,max_reported_trans_targets=1L,trans_set=c("core","extended")) {'''
core = replace_once(core, old_api, new_api, 'cater_mr API')

old_policy = '  primary_policy<-match.arg(primary_policy)\n'
new_policy = '''  primary_policy<-match.arg(primary_policy)
  trans_set<-match.arg(trans_set)
  if(primary_policy!="cis_anchor")warning("primary_policy is retained for compatibility but Lean v0.9 always uses cis as primary when available",call.=FALSE)
  if(isTRUE(allow_network_primary))warning("allow_network_primary is ignored by Lean v0.9; network MVMR is sensitivity-only",call.=FALSE)
'''
core = replace_once(core, old_policy, new_policy, 'Lean policy validation')

validation_anchor = '  if(length(instrument_p)!=1L||!is.finite(instrument_p)||instrument_p<=0||instrument_p>=1).cater_stop("instrument_p must be in (0,1)")\n'
validation_new = validation_anchor + '''  if(length(tf_anchor_p)!=1L||!is.finite(tf_anchor_p)||tf_anchor_p<=0||tf_anchor_p>=1).cater_stop("tf_anchor_p must be in (0,1)")
  if(length(max_reported_trans_targets)!=1L||!is.finite(max_reported_trans_targets)||max_reported_trans_targets<1||max_reported_trans_targets!=as.integer(max_reported_trans_targets)).cater_stop("max_reported_trans_targets must be a positive integer")
'''
core = replace_once(core, validation_anchor, validation_new, 'Lean argument validation')

old_contract = '  .cater_resolve_trans_contract(manifest,trans_eqtl,instrument_p,trans_reporting_p)\n'
new_contract = '  trans_mode<-.cater_resolve_trans_contract(manifest,trans_eqtl,instrument_p,trans_reporting_p)\n  trans_specificity_limit<-if(identical(trans_mode,"full_summary")) Inf else max_reported_trans_targets\n'
core = replace_once(core, old_contract, new_contract, 'trans mode capture')

# Replace summary schema directly in the original cater_mr implementation.
pattern = re.compile(r'^  empty_summary<-function\(x,status\).*$', re.M)
new_summary = '''  empty_summary<-function(x,status)data.frame(cell_type=cell_type,target=x,trait=trait,status=status,n_parent_tf=0,n_parent_tf_mhc_excluded=0,
    n_candidate_snp=0,n_selected_signal=0,n_cis_signal=0,n_trans_signal=0,n_trans_raw=0,n_trans_anchor_pass=0,n_trans_unique_parent=0,n_trans_core=0,n_trans_extended=0,n_parent_tf_core=0,
    n_ld_reference_missing=0,n_ld_allele_mismatch=0,n_ld_nonfinite=0,n_mr_iv=0,n_mr_cis_iv=0,n_mr_trans_iv=0,
    effective_F=NA,cis_effective_F=NA,trans_effective_F=NA,augmented_effective_F=NA,primary_effective_F=NA,
    trans_information_fraction=NA,augmented_precision_fraction=NA,augmented_se_reduction=NA,
    n_tf_anchor_tested=0,n_tf_anchor_fdr=NA_integer_,n_sibling_candidate=0,n_sibling_active=0,n_sibling_incomplete=0,n_sibling_qtl_missing=0,n_sibling_iv_missing=0,
    sibling_screen_performed=FALSE,sibling_screen_complete=NA,cis_trans_heterogeneity_p=NA,max_tf_locus_weight=NA,leave_one_tf_max_delta=NA,
    network_status=NA_character_,beta_network=NA,se_network=NA,p_network=NA,conditional_F_X=NA,conditional_F_converged=NA,mvmr_rank=NA,mvmr_condition=NA,mvmr_n_iv=NA_integer_,
    mvmr_union_requested=NA_integer_,mvmr_union_ld_available=NA_integer_,mvmr_union_complete_case=NA_integer_,mvmr_cross_effect_missing=NA_integer_,mvmr_max_r2=NA,mvmr_ld_gate_pass=NA,
    tf_anchor_p=tf_anchor_p,max_reported_trans_targets=trans_specificity_limit,trans_set_used=trans_set,
    primary_model=NA_character_,primary_beta=NA,primary_se=NA,primary_p=NA,stringsAsFactors=FALSE)'''
core, n = pattern.subn(new_summary, core, count=1)
if n != 1: raise RuntimeError(f'empty_summary replacement: {n}')

old_cand_line = '''    cand<-.cater_target_candidates(target,parents,ann,eqtl_dir,qtl_n,trans_hits,instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc);qtl<-cand$qtl;cmap<-cand$map;srow$n_candidate_snp<-nrow(qtl);if(!nrow(qtl)){srow$status<-"NO_SIGNIFICANT_CANDIDATE";summ[[target]]<-srow;next}'''
new_cand_line = '''    cand<-.cater_target_candidates(target,parents,ann,eqtl_dir,qtl_n,trans_hits,instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc,tf_anchor_p,trans_specificity_limit,trans_set)
    qtl<-cand$qtl;cmap<-cand$map;qc<-cand$qualification_counts%||%list(n_trans_raw=0L,n_trans_anchor_pass=0L,n_trans_unique_parent=0L,n_trans_core=0L,n_trans_extended=0L,n_parent_tf_core=0L)
    for(nm in names(qc))srow[[nm]]<-as.integer(qc[[nm]])
    if(!is.null(cand$qualification)&&nrow(cand$qualification))utils::write.table(cand$qualification,file.path(outdir,"mechanism",paste0(target,"_trans_qualification.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
    if(!is.null(cand$anchors)&&nrow(cand$anchors))utils::write.table(cand$anchors,file.path(outdir,"mechanism",paste0(target,"_trans_anchor_candidates.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
    srow$n_candidate_snp<-nrow(qtl);if(!nrow(qtl)){srow$status<-"NO_SIGNIFICANT_CANDIDATE";summ[[target]]<-srow;next}'''
core = replace_once(core, old_cand_line, new_cand_line, 'candidate call')

# Replace model fitting/output line: augmented is the canonical name; combined is retained only as a compatibility alias.
fit_pattern = re.compile(r'^    fits<-list\(\);for\(m in c\("cis","trans","combined"\)\).*$', re.M)
fit_new = '''    fits<-list()
    for(m in c("cis","trans","augmented")){
      d<-if(m=="augmented")h else h[h$source==m,,drop=FALSE];RR<-if(nrow(d)).cater_subset_ld(R,d$snp)else matrix(numeric(),0,0);f<-.cater_givw(d,RR);fits[[m]]<-f
      long[[paste(target,m,sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model=m,f,stringsAsFactors=FALSE)
    }
    fits$combined<-fits$augmented
    long[[paste(target,"combined",sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model="combined",fits$augmented,stringsAsFactors=FALSE)
    srow$cis_effective_F<-fits$cis$effective_F;srow$trans_effective_F<-fits$trans$effective_F;srow$augmented_effective_F<-fits$augmented$effective_F
    srow$primary_effective_F<-if(identical(fits$cis$status,"OK"))fits$cis$effective_F else NA_real_;srow$effective_F<-srow$primary_effective_F
    srow$trans_information_fraction<-.cater_trans_information_fraction(fits$cis,fits$augmented);srow$augmented_precision_fraction<-.cater_trans_mr_precision_fraction(fits$cis,fits$augmented)
    if(identical(fits$cis$status,"OK")&&identical(fits$augmented$status,"OK")&&is.finite(fits$cis$se)&&fits$cis$se>0&&is.finite(fits$augmented$se))srow$augmented_se_reduction<-1-fits$augmented$se/fits$cis$se
    ht<-.cater_cis_trans_het(h,R);srow$cis_trans_heterogeneity_p<-ht["p"];ldgn<-.cater_locus_diagnostics(h,R);srow$max_tf_locus_weight<-ldgn$max_weight;srow$leave_one_tf_max_delta<-ldgn$leave_one_max_delta
    if(nrow(ldgn$table))utils::write.table(ldgn$table,file.path(outdir,"diagnostics",paste0(target,"_loci.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)'''
core, n = fit_pattern.subn(fit_new, core, count=1)
if n != 1: raise RuntimeError(f'fit block replacement: {n}')

# Remove the old post-outcome TF-anchor/hotspot diagnostics; qualification now occurs before LD/outcome use.
core, n1 = re.subn(r'^    anchor<-.*\n', '', core, count=1, flags=re.M)
core, n2 = re.subn(r'^    hot<-.*\n', '', core, count=1, flags=re.M)
if n1 != 1 or n2 != 1: raise RuntimeError(f'posthoc removal anchor={n1} hotspot={n2}')

core_path.write_text(core)

# The entrypoint now loads only the directly modified core plus LD helpers.
entry = entry_path.read_text()
entry = entry.replace('source(file.path("R", "CATER_MR_lean.R"), local = FALSE)\n', '')
entry_path.write_text(entry)

# CI parses/tests the direct core; there is no override layer.
wf = workflow_path.read_text().replace('          Rscript -e "parse(file=\'R/CATER_MR_lean.R\')"\n', '')
workflow_path.write_text(wf)

# Update tests to assert direct implementation semantics.
v09 = ROOT / 'tests' / 'v09_lean_trans_smoke.R'
t = v09.read_text()
t = t.replace('old <- .CATER_LEAN_STATE$active; .CATER_LEAN_STATE$active <- TRUE\npd <- .cater_primary_decision(fits,NULL,"X",TRUE,TRUE,TRUE,0L,"screened_cater")\n.CATER_LEAN_STATE$active <- old\n',
              'pd <- .cater_primary_decision(fits,NULL,"X",TRUE,TRUE,TRUE,0L,"screened_cater")\n')
post_start = t.index('# Post-processing exposes model-specific F statistics')
post_end = t.index('unlink(td,recursive=TRUE,force=TRUE)', post_start)
replacement = '''# The public implementation is direct: no override state/file exists, and the original
# cater_mr body owns model-specific F and augmented precision reporting.
stopifnot(!exists(".CATER_LEAN_STATE",inherits=FALSE))
stopifnot(!file.exists("R/CATER_MR_lean.R"))
body_txt <- paste(deparse(body(cater_mr)),collapse="\\n")
stopifnot(grepl("cis_effective_F",body_txt,fixed=TRUE))
stopifnot(grepl("augmented_effective_F",body_txt,fixed=TRUE))
stopifnot(grepl("augmented_precision_fraction",body_txt,fixed=TRUE))
stopifnot(grepl("augmented_se_reduction",body_txt,fixed=TRUE))
fc <- mkfit(F=12); fa <- mkfit(se=.04,F=22)
stopifnot(isTRUE(all.equal(1-fa$se/fc$se,0.2)))

'''
t = t[:post_start] + replacement + t[post_end:]
t = t.replace('unlink(td,recursive=TRUE,force=TRUE); unlink(out,recursive=TRUE,force=TRUE)', 'unlink(td,recursive=TRUE,force=TRUE)')
v09.write_text(t)

# v0.8 regression file becomes an input-contract compatibility test under direct v0.9 semantics.
v08 = ROOT / 'tests' / 'v08_significant_trans_smoke.R'
t = v08.read_text()
# Add the parent-TF full cis file required by the new direct gate.
needle = 'utils::write.table(cis,con,sep="\\t",quote=FALSE,row.names=FALSE); close(con)\n'
insert = needle + '''tf1 <- data.frame(SNP="rsT",CHR=2,BP=5000,A1="A",A2="G",EAF=.2,BETA=.25,SE=.02,P=1e-25,N=500)
con <- gzfile(file.path(td,"TF1.txt.gz"),"wt")
utils::write.table(tf1,con,sep="\\t",quote=FALSE,row.names=FALSE); close(con)
'''
t = replace_once(t, needle, insert, 'v08 TF anchor fixture')
# Network primary is no longer possible.
t = t.replace('stopifnot(d_blocked$model=="cis",d_allowed$model=="network")', 'stopifnot(d_blocked$model=="cis",d_allowed$model=="cis")')
v08.write_text(t)

# Core smoke: primary is always cis when available and never trans/network when cis is absent.
smoke = ROOT / 'tests' / 'smoke.R'
t = smoke.read_text()
t = t.replace('stopifnot(is.na(dec0$model),is.na(dec0$p),dec0$status=="TRANS_PLEIOTROPY_UNRESOLVED_NO_CIS")', 'stopifnot(is.na(dec0$model),is.na(dec0$p),dec0$status=="NO_CIS_PRIMARY")')
t = t.replace('stopifnot(is.na(dec1$model),is.na(dec1$p),dec1$status=="SIBLING_SCREEN_INCOMPLETE_NO_CIS")', 'stopifnot(is.na(dec1$model),is.na(dec1$p),dec1$status=="NO_CIS_PRIMARY")')
t = t.replace('dec2$status=="SIBLING_SCREEN_INCOMPLETE_CIS_FALLBACK")', 'dec2$status=="OK_CIS_PRIMARY")')
t = t.replace('stopifnot(dec3$model=="cis",dec3$status=="TRANS_UNSCREENED_CIS_FALLBACK")', 'stopifnot(dec3$model=="cis",dec3$status=="OK_CIS_PRIMARY")')
smoke.write_text(t)

# Documentation/API test must reject a second implementation file.
doc = ROOT / 'tests' / 'documentation_cleanup_smoke.R'
t = doc.read_text()
t = t.replace('"R/CATER_MR_core.R","R/CATER_MR_lean.R","R/LD_helpers.R")', '"R/CATER_MR_core.R","R/LD_helpers.R")')
t = t.replace('stopifnot(all(file.exists(required)))\n', 'stopifnot(all(file.exists(required)))\nstopifnot(!file.exists("R/CATER_MR_lean.R"))\n')
doc.write_text(t)

# Delete the override source itself.
if lean_path.exists(): lean_path.unlink()

print('Direct core refactor applied successfully')
