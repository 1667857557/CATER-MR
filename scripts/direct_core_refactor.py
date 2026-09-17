from pathlib import Path
import re

ROOT = Path('.')
core_path = ROOT / 'R' / 'CATER_MR_core.R'
core = core_path.read_text()


def function_span(text, name):
    marker = name + ' <- function'
    start = text.find(marker)
    if start < 0:
        raise RuntimeError(f'function not found: {name}')
    brace = text.find('{', start)
    if brace < 0:
        raise RuntimeError(f'opening brace not found: {name}')
    depth = 0
    quote = None
    escape = False
    comment = False
    for i in range(brace, len(text)):
        ch = text[i]
        if comment:
            if ch == '\n':
                comment = False
            continue
        if quote is not None:
            if escape:
                escape = False
            elif ch == '\\':
                escape = True
            elif ch == quote:
                quote = None
            continue
        if ch in ('"', "'"):
            quote = ch
        elif ch == '#':
            comment = True
        elif ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return start, i + 1
    raise RuntimeError(f'unclosed function: {name}')


def replace_function(text, name, replacement):
    a, b = function_span(text, name)
    return text[:a] + replacement.rstrip() + text[b:]


def one(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f'{label}: expected one match, found {n}')
    return text.replace(old, new, 1)


# Single implementation: version lives in the original core.
core = one(core, '# CATER-MR v0.8\n', '# CATER-MR v0.9\n', 'header')
core = one(core, '.CATER_VERSION <- "0.8.2"', '.CATER_VERSION <- "0.9.0"', 'version')
core = core.replace('ignored by CATER-MR v0.8:', 'ignored by CATER-MR v0.9:')
core = core.replace('for the v0.8 contract', 'for the v0.9 contract')

# Outcome-independent trans qualification helpers. Parent cis files are cached once
# per run, and reported-target counts are precomputed once rather than rescanning.
helpers = r'''
.cater_cached_read_gene <- function(gene,eqtl_dir,qtl_n,cache=NULL) {
  if(is.null(cache)) return(.cater_read_gene(gene,eqtl_dir,qtl_n))
  key<-as.character(gene)
  if(exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache,inherits=FALSE))
  d<-.cater_read_gene(gene,eqtl_dir,qtl_n);assign(key,d,envir=cache);d
}
.cater_preld_tf_anchor <- function(tf,snp,target_row,eqtl_dir,qtl_n,qtl_cache=NULL) {
  d<-.cater_cached_read_gene(tf,eqtl_dir,qtl_n,qtl_cache)
  if(is.null(d)) return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=NA_real_,p_tf=NA_real_,status="QTL_NOT_AVAILABLE",stringsAsFactors=FALSE))
  d<-d[d$snp==snp,,drop=FALSE]
  if(!nrow(d)) return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=NA_real_,p_tf=NA_real_,status="SNP_NOT_AVAILABLE",stringsAsFactors=FALSE))
  same<-d$a1[1]==target_row$a1[1]&&d$a2[1]==target_row$a2[1]
  swap<-d$a1[1]==target_row$a2[1]&&d$a2[1]==target_row$a1[1]
  if(!(same||swap)) return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=d$se[1],p_tf=d$p[1],status="ALLELE_MISMATCH",stringsAsFactors=FALSE))
  data.frame(tf=tf,beta_tf=if(swap)-d$beta[1]else d$beta[1],se_tf=d$se[1],p_tf=d$p[1],status="OK",stringsAsFactors=FALSE)
}
.cater_precompute_trans_target_counts <- function(trans_hits,trans_mode,instrument_p,trans_reporting_p) {
  if(is.null(trans_hits)||!nrow(trans_hits)) return(integer())
  threshold<-if(identical(trans_mode,"full_summary")) instrument_p else trans_reporting_p
  if(length(threshold)!=1L||!is.finite(threshold)||threshold<=0||threshold>1) .cater_stop("Invalid trans target-count threshold")
  d<-trans_hits[is.finite(trans_hits$p)&trans_hits$p<threshold,c("snp","gene"),drop=FALSE]
  if(!nrow(d)) return(integer())
  d<-unique(d);tt<-table(d$snp);out<-as.integer(tt);names(out)<-names(tt);out
}
.cater_qualify_trans_candidates <- function(target,qtl,map,eqtl_dir,qtl_n,trans_target_counts=integer(),
                                             tf_anchor_p=5e-8,max_reported_trans_targets=1L,
                                             trans_set=c("core","extended"),qtl_cache=NULL,
                                             trans_mode="significant_only") {
  trans_set<-match.arg(trans_set)
  if(length(tf_anchor_p)!=1L||!is.finite(tf_anchor_p)||tf_anchor_p<=0||tf_anchor_p>=1) .cater_stop("tf_anchor_p must be in (0,1)")
  if(length(max_reported_trans_targets)!=1L||!is.finite(max_reported_trans_targets)||max_reported_trans_targets<1||max_reported_trans_targets!=as.integer(max_reported_trans_targets)) .cater_stop("max_reported_trans_targets must be a positive integer")
  zero<-list(n_trans_raw=0L,n_trans_anchor_pass=0L,n_trans_unique_parent=0L,n_trans_core=0L,n_trans_extended=0L,n_parent_tf_core=0L)
  if(!nrow(qtl)||!nrow(map)) return(list(qtl=qtl,map=map,table=data.frame(),anchors=data.frame(),counts=zero))
  ti<-which(map$source=="trans")
  if(!length(ti)) return(list(qtl=qtl,map=map,table=data.frame(),anchors=data.frame(),counts=zero))
  rows<-vector("list",length(ti));anchors<-list()
  for(ii in seq_along(ti)) {
    mi<-ti[ii];s<-map$snp[mi];tr<-qtl[match(s,qtl$snp),,drop=FALSE];parents<-.cater_parent_tokens(map$parent_tf[mi])[[1L]]
    aa<-if(length(parents))do.call(rbind,lapply(parents,function(tf).cater_preld_tf_anchor(tf,s,tr,eqtl_dir,qtl_n,qtl_cache)))else data.frame()
    if(nrow(aa)) {aa$target<-target;aa$snp<-s;anchors[[length(anchors)+1L]]<-aa;pass<-aa$status=="OK"&is.finite(aa$p_tf)&aa$p_tf<tf_anchor_p;qtf<-unique(aa$tf[pass])} else qtf<-character()
    nqual<-length(qtf)
    nrep<-if(s%in%names(trans_target_counts))as.integer(trans_target_counts[[s]])else if(identical(trans_mode,"legacy_full_summary"))NA_integer_ else 0L
    specificity<-is.na(nrep)||nrep<=max_reported_trans_targets;extended<-nqual==1L;core_ok<-extended&&specificity
    reason<-if(nqual==0L)"TRANS_NO_TF_CIS_ANCHOR"else if(nqual>1L)"TRANS_AMBIGUOUS_PARENT"else if(!specificity)"TRANS_REPORTED_HOTSPOT"else"TRANS_CORE"
    chosen<-if(nqual==1L)qtf[[1L]]else NA_character_;aone<-if(nqual==1L)aa[match(chosen,aa$tf),,drop=FALSE]else NULL
    rows[[ii]]<-data.frame(target=target,snp=s,raw_parent_tf=map$parent_tf[mi],qualified_parent=chosen,n_qualifying_parent=nqual,
      beta_tf=if(is.null(aone))NA_real_ else aone$beta_tf[1],se_tf=if(is.null(aone))NA_real_ else aone$se_tf[1],p_tf=if(is.null(aone))NA_real_ else aone$p_tf[1],
      n_reported_trans_targets=nrep,anchor_pass=nqual>=1L,unique_parent=extended,specificity_pass=specificity,eligible_extended=extended,eligible_core=core_ok,reason=reason,stringsAsFactors=FALSE)
  }
  tab<-do.call(rbind,rows);rownames(tab)<-NULL;anchor_tab<-if(length(anchors))do.call(rbind,anchors)else data.frame();use<-if(trans_set=="core")tab$eligible_core else tab$eligible_extended
  keep_snps<-tab$snp[use];cis_idx<-which(map$source=="cis");trans_idx<-ti[map$snp[ti]%in%keep_snps];keep_idx<-c(cis_idx,trans_idx)
  qout<-qtl[match(map$snp[keep_idx],qtl$snp),,drop=FALSE];mout<-map[keep_idx,,drop=FALSE]
  if(length(trans_idx)) {dd<-tab[match(mout$snp[mout$source=="trans"],tab$snp),,drop=FALSE];mout$parent_tf[mout$source=="trans"]<-dd$qualified_parent;prov<-if(trans_set=="core")"REPORTED_TRANS_CORE"else"REPORTED_TRANS_EXTENDED";mout$provenance[mout$source=="trans"]<-prov;qout$provenance[qout$snp%in%mout$snp[mout$source=="trans"]]<-prov}
  counts<-list(n_trans_raw=nrow(tab),n_trans_anchor_pass=sum(tab$anchor_pass),n_trans_unique_parent=sum(tab$unique_parent),n_trans_core=sum(tab$eligible_core),n_trans_extended=sum(tab$eligible_extended),n_parent_tf_core=length(unique(tab$qualified_parent[tab$eligible_core&!is.na(tab$qualified_parent)])))
  list(qtl=qout,map=mout,table=tab,anchors=anchor_tab,counts=counts)
}
'''.strip('\n')
marker='.cater_target_candidates <- function'
if '.cater_cached_read_gene <- function' not in core:
    core=one(core,marker,helpers+'\n\n'+marker,'insert qualification helpers')

# Directly change the existing primary-decision function: no wrapper/alias.
primary=r'''.cater_primary_decision <- function(fits,net=NULL,target=NULL,has_trans=FALSE,sibling_screen_performed=FALSE,sibling_screen_complete=TRUE,n_active_siblings=0L,primary_policy=c("cis_anchor","screened_cater"),allow_trans_only_primary=FALSE,n_trans_tf_loci=0L,min_trans_tf_loci=3L,allow_network_primary=FALSE,sibling_screen_independent=FALSE,sibling_screen_testable=TRUE){
  primary_policy<-match.arg(primary_policy)
  out<-list(model=NA_character_,beta=NA_real_,se=NA_real_,p=NA_real_,status="NO_CIS_PRIMARY",evidence_status=if(has_trans)"TRANS_SENSITIVITY_ONLY"else"NO_CIS_INSTRUMENT")
  if(!is.null(fits$cis)&&identical(fits$cis$status,"OK")){out$model<-"cis";out$beta<-fits$cis$beta;out$se<-fits$cis$se;out$p<-fits$cis$p;out$status<-"OK_CIS_PRIMARY";out$evidence_status<-"CIS_ANCHOR_PRIMARY"}
  out
}'''
core=replace_function(core,'.cater_primary_decision',primary)

# Modify only the original cater_mr() body/signature.
a,b=function_span(core,'cater_mr');mr=core[a:b]
mr=one(mr,'sibling_fdr=0.05,enable_sibling_screen=TRUE,enable_mvmr=TRUE,','sibling_fdr=0.05,enable_sibling_screen=FALSE,enable_mvmr=FALSE,','advanced defaults')
mr=one(mr,'input_manifest=NULL,strict_input_contract=FALSE) {','input_manifest=NULL,strict_input_contract=FALSE,\n                     tf_anchor_p=instrument_p,max_reported_trans_targets=1L,trans_set=c("core","extended")) {','append lean args')
mr=one(mr,'  primary_policy<-match.arg(primary_policy)\n','  primary_policy<-match.arg(primary_policy);trans_set<-match.arg(trans_set)\n  if(primary_policy!="cis_anchor")warning("primary_policy is compatibility-only in CATER-MR v0.9; cis remains primary",call.=FALSE)\n  if(isTRUE(allow_network_primary))warning("allow_network_primary is ignored in CATER-MR v0.9; network fits are sensitivity-only",call.=FALSE)\n','policy')
needle='  if(length(instrument_p)!=1L||!is.finite(instrument_p)||instrument_p<=0||instrument_p>=1).cater_stop("instrument_p must be in (0,1)")\n'
mr=one(mr,needle,needle+'  if(length(tf_anchor_p)!=1L||!is.finite(tf_anchor_p)||tf_anchor_p<=0||tf_anchor_p>=1).cater_stop("tf_anchor_p must be in (0,1)")\n  if(length(max_reported_trans_targets)!=1L||!is.finite(max_reported_trans_targets)||max_reported_trans_targets<1||max_reported_trans_targets!=as.integer(max_reported_trans_targets)).cater_stop("max_reported_trans_targets must be a positive integer")\n','lean validation')
mr=one(mr,'  .cater_resolve_trans_contract(manifest,trans_eqtl,instrument_p,trans_reporting_p)\n','  trans_mode<-.cater_resolve_trans_contract(manifest,trans_eqtl,instrument_p,trans_reporting_p)\n','trans mode')
std='  outcome<-.cater_standardize_sumstats(outcome,label="outcome",require_position=FALSE);trans_hits<-.cater_standardize_trans_hits(trans_eqtl,trans_gene_col,qtl_n)\n'
mr=one(mr,std,std+'  trans_target_counts<-.cater_precompute_trans_target_counts(trans_hits,trans_mode,instrument_p,trans_reporting_p);qtl_cache<-new.env(parent=emptyenv())\n','precompute/cache')

# Replace summary row schema.
mr,n=re.subn(r'^  empty_summary<-function\(x,status\).*$', '''  empty_summary<-function(x,status)data.frame(cell_type=cell_type,target=x,trait=trait,status=status,n_parent_tf=0,n_parent_tf_mhc_excluded=0,n_candidate_snp=0,n_trans_raw=0,n_trans_anchor_pass=0,n_trans_unique_parent=0,n_trans_core=0,n_trans_extended=0,n_parent_tf_core=0,n_selected_signal=0,n_cis_signal=0,n_trans_signal=0,n_ld_reference_missing=0,n_ld_allele_mismatch=0,n_ld_nonfinite=0,n_mr_iv=0,n_mr_cis_iv=0,n_mr_trans_iv=0,effective_F=NA_real_,cis_effective_F=NA_real_,trans_effective_F=NA_real_,augmented_effective_F=NA_real_,primary_effective_F=NA_real_,trans_information_fraction=NA_real_,augmented_precision_fraction=NA_real_,augmented_se_reduction=NA_real_,n_tf_anchor_tested=0,n_tf_anchor_fdr=NA_integer_,n_sibling_candidate=0,n_sibling_active=0,n_sibling_incomplete=0,n_sibling_qtl_missing=0,n_sibling_iv_missing=0,sibling_screen_performed=FALSE,sibling_screen_complete=NA,cis_trans_heterogeneity_p=NA_real_,max_tf_locus_weight=NA_real_,leave_one_tf_max_delta=NA_real_,network_status=NA_character_,beta_network=NA_real_,se_network=NA_real_,p_network=NA_real_,conditional_F_X=NA_real_,conditional_F_converged=NA,mvmr_rank=NA_real_,mvmr_condition=NA_real_,mvmr_n_iv=NA_integer_,mvmr_union_requested=NA_integer_,mvmr_union_ld_available=NA_integer_,mvmr_union_complete_case=NA_integer_,mvmr_cross_effect_missing=NA_integer_,mvmr_max_r2=NA_real_,mvmr_ld_gate_pass=NA,primary_model=NA_character_,primary_beta=NA_real_,primary_se=NA_real_,primary_p=NA_real_,stringsAsFactors=FALSE)''',mr,count=1,flags=re.M)
if n!=1: raise RuntimeError('summary schema replacement failed')

old_cand='    cand<-.cater_target_candidates(target,parents,ann,eqtl_dir,qtl_n,trans_hits,instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc);qtl<-cand$qtl;cmap<-cand$map;srow$n_candidate_snp<-nrow(qtl);if(!nrow(qtl)){srow$status<-"NO_SIGNIFICANT_CANDIDATE";summ[[target]]<-srow;next}\n'
new_cand='''    cand<-.cater_target_candidates(target,parents,ann,eqtl_dir,qtl_n,trans_hits,instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc);raw_qtl<-cand$qtl;raw_map<-cand$map;srow$n_candidate_snp<-nrow(raw_qtl)
    qual<-.cater_qualify_trans_candidates(target,raw_qtl,raw_map,eqtl_dir,qtl_n,trans_target_counts,tf_anchor_p,max_reported_trans_targets,trans_set,qtl_cache,trans_mode);qtl<-qual$qtl;cmap<-qual$map
    for(nm in names(qual$counts))srow[[nm]]<-as.integer(qual$counts[[nm]])
    if(nrow(qual$table))utils::write.table(qual$table,file.path(outdir,"mechanism",paste0(target,"_trans_qualification.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
    if(nrow(qual$anchors)){utils::write.table(qual$anchors,file.path(outdir,"mechanism",paste0(target,"_trans_anchor_candidates.tsv")),sep="\\t",quote=FALSE,row.names=FALSE);srow$n_tf_anchor_tested<-sum(qual$anchors$status=="OK")}
    if(!nrow(qtl)){srow$status<-if(nrow(raw_qtl))"NO_CANDIDATE_AFTER_TRANS_QUALIFICATION"else"NO_SIGNIFICANT_CANDIDATE";summ[[target]]<-srow;next}
'''
mr=one(mr,old_cand,new_cand,'candidate qualification')

old_fit='    fits<-list();for(m in c("cis","trans","combined")){d<-if(m=="combined")h else h[h$source==m,,drop=FALSE];RR<-if(nrow(d)).cater_subset_ld(R,d$snp)else matrix(numeric(),0,0);f<-.cater_givw(d,RR);fits[[m]]<-f;long[[paste(target,m,sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model=m,f,stringsAsFactors=FALSE)};srow$effective_F<-fits$combined$effective_F;srow$trans_information_fraction<-.cater_trans_information_fraction(fits$cis,fits$combined);ht<-.cater_cis_trans_het(h,R);srow$cis_trans_heterogeneity_p<-ht["p"];ldgn<-.cater_locus_diagnostics(h,R);srow$max_tf_locus_weight<-ldgn$max_weight;srow$leave_one_tf_max_delta<-ldgn$leave_one_max_delta;if(nrow(ldgn$table))utils::write.table(ldgn$table,file.path(outdir,"diagnostics",paste0(target,"_loci.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)\n'
new_fit='''    fits<-list();for(m in c("cis","trans","combined")){d<-if(m=="combined")h else h[h$source==m,,drop=FALSE];RR<-if(nrow(d)).cater_subset_ld(R,d$snp)else matrix(numeric(),0,0);f<-.cater_givw(d,RR);fits[[m]]<-f;long[[paste(target,m,sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model=m,f,stringsAsFactors=FALSE)}
    long[[paste(target,"augmented",sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model="augmented",fits$combined,stringsAsFactors=FALSE)
    srow$cis_effective_F<-fits$cis$effective_F;srow$trans_effective_F<-fits$trans$effective_F;srow$augmented_effective_F<-fits$combined$effective_F;srow$primary_effective_F<-if(identical(fits$cis$status,"OK"))fits$cis$effective_F else NA_real_;srow$effective_F<-srow$primary_effective_F
    srow$trans_information_fraction<-.cater_trans_information_fraction(fits$cis,fits$combined);srow$augmented_precision_fraction<-.cater_trans_mr_precision_fraction(fits$cis,fits$combined)
    if(identical(fits$cis$status,"OK")&&identical(fits$combined$status,"OK")&&is.finite(fits$cis$se)&&fits$cis$se>0&&is.finite(fits$combined$se))srow$augmented_se_reduction<-1-fits$combined$se/fits$cis$se
    ht<-.cater_cis_trans_het(h,R);srow$cis_trans_heterogeneity_p<-ht["p"];ldgn<-.cater_locus_diagnostics(h,R);srow$max_tf_locus_weight<-ldgn$max_weight;srow$leave_one_tf_max_delta<-ldgn$leave_one_max_delta;if(nrow(ldgn$table))utils::write.table(ldgn$table,file.path(outdir,"diagnostics",paste0(target,"_loci.tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
'''
mr=one(mr,old_fit,new_fit,'fit outputs')
# Old post-hoc anchor/hotspot scans are superseded by pre-LD qualification.
mr=re.sub(r'^    anchor<-\.cater_tf_anchor.*\n','',mr,count=1,flags=re.M)
mr=re.sub(r'^    hot<-\.cater_trans_hotspot.*\n','',mr,count=1,flags=re.M)
core=core[:a]+mr+core[b:]
core_path.write_text(core)

# Remove the override source and file: one implementation path only.
entry=Path('CATER_MR.R');entry.write_text(entry.read_text().replace('source(file.path("R", "CATER_MR_lean.R"), local = FALSE)\n',''))
lean=Path('R/CATER_MR_lean.R')
if lean.exists(): lean.unlink()

# CI parses the original core directly.
wf=Path('.github/workflows/r-smoke.yml');wf.write_text(wf.read_text().replace('          Rscript -e "parse(file=\'R/CATER_MR_lean.R\')"\n',''))

# Direct-core regression tests.
p=Path('tests/documentation_cleanup_smoke.R');s=p.read_text().replace('"R/CATER_MR_core.R","R/CATER_MR_lean.R","R/LD_helpers.R"','"R/CATER_MR_core.R","R/LD_helpers.R"')
s=s.replace('obsolete_files <- c("METHODOLOGY_V08.md","LD_DIAGNOSIS.md",','obsolete_files <- c("METHODOLOGY_V08.md","LD_DIAGNOSIS.md","R/CATER_MR_lean.R",');p.write_text(s)

p=Path('tests/v08_significant_trans_smoke.R');s=p.read_text().replace('# Network primary requires explicit policy + conditional-F acceptance + both public gates.','# Lean v0.9 primary remains cis even when an advanced network sensitivity is available.').replace('stopifnot(d_blocked$model=="cis",d_allowed$model=="network")','stopifnot(d_blocked$model=="cis",d_allowed$model=="cis")');p.write_text(s)

p=Path('tests/smoke.R');s=p.read_text();x=s.find('# Unresolved or incompletely screened trans estimates never populate primary_*.');y=s.find('# Deterministic local MVMR',x)
if x<0 or y<0: raise RuntimeError('smoke primary block not found')
new_primary='''# Lean v0.9 keeps cis as the only primary model; trans/network fits are sensitivity analyses.\nmkfit <- function(status,beta=.2,se=.05,p=.01,n_iv=2L,information=10)\n  data.frame(n_iv=n_iv,beta=beta,se=se,p=p,Q=NA,Q_p=NA,information=information,effective_F=NA,mean_F=NA,min_F=NA,status=status)\nfits0 <- list(cis=mkfit("NO_IV",n_iv=0L,information=NA),trans=mkfit("OK"),combined=mkfit("OK"))\ndec0 <- .cater_primary_decision(fits0,has_trans=TRUE,primary_policy="screened_cater")\nstopifnot(is.na(dec0$model),dec0$status=="NO_CIS_PRIMARY")\nfits1 <- fits0; fits1$cis <- mkfit("OK",beta=.11,p=.02)\ndec1 <- .cater_primary_decision(fits1,has_trans=TRUE,n_active_siblings=1L,primary_policy="screened_cater",allow_network_primary=TRUE)\nstopifnot(dec1$model=="cis",abs(dec1$beta-.11)<1e-12,dec1$status=="OK_CIS_PRIMARY")\n\n'''
p.write_text(s[:x]+new_primary+s[y:])

Path('tests/v09_lean_trans_smoke.R').write_text(r'''source("CATER_MR.R")
stopifnot(identical(.CATER_VERSION,"0.9.0"),!file.exists("R/CATER_MR_lean.R"))
fml<-formals(cater_mr);stopifnot(identical(fml$enable_sibling_screen,FALSE),identical(fml$enable_mvmr,FALSE),all(c("tf_anchor_p","max_reported_trans_targets","trans_set")%in%names(fml)))

td<-tempfile("cater_v09_");dir.create(td)
write_gene<-function(gene,dat){con<-gzfile(file.path(td,paste0(gene,".txt.gz")),"wt");utils::write.table(dat,con,sep="\t",quote=FALSE,row.names=FALSE);close(con)}
base<-data.frame(SNP=c("rsCore","rsHot","rsAmb"),CHR=2,BP=c(5000,5100,5200),A1="A",A2="G",EAF=.2,BETA=c(.20,.18,.16),SE=.02,P=c(1e-20,1e-18,1e-16),N=500,stringsAsFactors=FALSE)
write_gene("TF1",base);tf2<-base[3,,drop=FALSE];tf2$BETA<-.12;tf2$P<-1e-12;write_gene("TF2",tf2)
qtl<-.cater_standardize_sumstats(base,n_default=500,label="target trans");qtl$provenance<-"REPORTED_TRANS"
map<-data.frame(snp=qtl$snp,source="trans",parent_tf=c("TF1","TF1","TF1;TF2"),locus_id=c("TF:TF1","TF:TF1","TF:TF1+TF2"),provenance="REPORTED_TRANS",stringsAsFactors=FALSE)
trans<-data.frame(GENE=c("X","X","Z","X"),SNP=c("rsCore","rsHot","rsHot","rsAmb"),CHR=2,BP=c(5000,5100,5100,5200),A1="A",A2="G",EAF=.2,BETA=c(.20,.18,.11,.16),SE=.02,P=c(1e-20,1e-18,1e-12,1e-16),N=500,stringsAsFactors=FALSE)
th<-.cater_standardize_trans_hits(trans,qtl_n=500);counts<-.cater_precompute_trans_target_counts(th,"significant_only",5e-8,5e-8);cache<-new.env(parent=emptyenv())
core<-.cater_qualify_trans_candidates("X",qtl,map,td,500,counts,5e-8,1L,"core",cache,"significant_only")
stopifnot(identical(core$qtl$snp,"rsCore"),core$counts$n_trans_raw==3L,core$counts$n_trans_anchor_pass==3L,core$counts$n_trans_unique_parent==2L,core$counts$n_trans_core==1L,core$counts$n_trans_extended==2L,core$counts$n_parent_tf_core==1L)
stopifnot(core$table$reason[core$table$snp=="rsHot"]=="TRANS_REPORTED_HOTSPOT",core$table$reason[core$table$snp=="rsAmb"]=="TRANS_AMBIGUOUS_PARENT",setequal(ls(cache),c("TF1","TF2")))
extended<-.cater_qualify_trans_candidates("X",qtl,map,td,500,counts,5e-8,1L,"extended",cache,"significant_only");stopifnot(setequal(extended$qtl$snp,c("rsCore","rsHot")))

# Codex P1: full-summary rows below significance do not inflate the target burden.
full<-data.frame(GENE=c("X","Z","W"),SNP="rsCore",CHR=2,BP=5000,A1="A",A2="G",EAF=.2,BETA=c(.2,.01,.01),SE=.02,P=c(1e-20,.5,.8),N=500,stringsAsFactors=FALSE)
full_std<-.cater_standardize_trans_hits(full,qtl_n=500);full_counts<-.cater_precompute_trans_target_counts(full_std,"full_summary",5e-8,1)
stopifnot(identical(unname(full_counts["rsCore"]),1L))
full_core<-.cater_qualify_trans_candidates("X",qtl[1,,drop=FALSE],map[1,,drop=FALSE],td,500,full_counts,5e-8,1L,"core",cache,"full_summary");stopifnot(identical(full_core$qtl$snp,"rsCore"))

q_no<-qtl[1,,drop=FALSE];q_no$snp<-"rsNoAnchor";m_no<-data.frame(snp="rsNoAnchor",source="trans",parent_tf="TF1",locus_id="TF:TF1",provenance="REPORTED_TRANS",stringsAsFactors=FALSE)
z_no<-.cater_qualify_trans_candidates("X",q_no,m_no,td,500,integer(),5e-8,1L,"core",cache,"significant_only");stopifnot(nrow(z_no$qtl)==0L,z_no$table$reason=="TRANS_NO_TF_CIS_ANCHOR")

mkfit<-function(status="OK",b=.2,se=.05,p=.01,F=20)data.frame(n_iv=2L,beta=b,se=se,p=p,Q=NA,Q_p=NA,information=40,joint_wald_per_df=F,effective_F=F,precision_information=25,mean_F=F,min_F=F,status=status)
fits<-list(cis=mkfit(),trans=mkfit(b=.4),combined=mkfit(b=.3));pd<-.cater_primary_decision(fits,NULL,"X",TRUE,TRUE,TRUE,1L,"screened_cater",allow_network_primary=TRUE);stopifnot(pd$model=="cis",pd$status=="OK_CIS_PRIMARY")
unlink(td,recursive=TRUE,force=TRUE);cat("CATER-MR Lean v0.9 direct-core smoke tests passed\n")
''')

# Update docs: full-summary burden counts only significant associations.
for fn in ('README.md','EQTL_INPUT.md','MATHEMATICS.md'):
    p=Path(fn);s=p.read_text();s=s.replace('With `trans_data_mode="full_summary"`, the observed-target-count filter is disabled because the file may contain nonsignificant tested rows.','With `trans_data_mode="full_summary"`, target burden counts only SNP-gene associations with `P < instrument_p`; nonsignificant tested rows do not count as reported targets.')
    p.write_text(s)

core_path.write_text(core)
