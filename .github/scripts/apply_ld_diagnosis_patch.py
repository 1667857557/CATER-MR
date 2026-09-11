from pathlib import Path

core = Path('R/CATER_MR_core_v05.R')
s = core.read_text()

def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f'{label}: expected 1 match, found {n}')
    return text.replace(old, new, 1)

s = replace_once(s, '.CATER_VERSION <- "0.7.0"', '.CATER_VERSION <- "0.7.1"', 'version')

start = s.index('.cater_run_cojo <- function')
end = s.index('\n\n.cater_joint_ld <- function', start)
new_run = r'''.cater_run_cojo <- function(qtl,candidate_map,ld_bfile,manc_cojo_bin,cojo_p,
                            cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,
                            plink_bin="plink",enable_ld_diagnosis=TRUE,
                            ld_diag_loglr=2,ld_diag_abs_z=2) {
  if (length(ld_bfile)!=1L || !nzchar(ld_bfile)) .cater_stop("Current CATER-MR expects one COJO LD cohort")
  if (!all(file.exists(paste0(ld_bfile,c(".bed",".bim",".fam"))))) .cater_stop("ld_bfile must be a PLINK prefix: %s",ld_bfile)
  exe <- Sys.which(manc_cojo_bin); if (!nzchar(exe)) .cater_stop("Cannot find Manc-COJO executable '%s'",manc_cojo_bin)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE); .cater_remove_prefix_outputs(prefix)
  diag_res <- list(candidate_map=candidate_map,diagnostics=data.frame(),removed=character())
  if (isTRUE(enable_ld_diagnosis)) {
    diag_res <- .cater_ld_diagnosis_filter(qtl,candidate_map,ld_bfile,plink_bin,cojo_threads,prefix,
      loglr_cutoff=ld_diag_loglr,abs_z_cutoff=ld_diag_abs_z,verbose=verbose)
    candidate_map <- diag_res$candidate_map
    if (nrow(diag_res$diagnostics)) utils::write.table(diag_res$diagnostics,paste0(prefix,".ld_diagnosis.tsv"),
      sep="\t",quote=FALSE,row.names=FALSE)
  }
  finish <- function(selected=data.frame(),selected_step=selected,joint_removed=character(),ld=matrix(numeric(),0,0))
    list(selected=selected,selected_step=selected_step,joint_removed=joint_removed,ld=ld,
         ld_diagnosis=diag_res$diagnostics,ld_diagnosis_removed=diag_res$removed)
  if (!nrow(candidate_map)) return(finish())
  ma <- paste0(prefix,".sumstat"); extract <- paste0(prefix,".candidate.snplist")
  .cater_write_cojo_ma(qtl,ma); writeLines(unique(candidate_map$snp),extract)
  sp <- paste0(prefix,".select")
  args <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",extract,"--cojo-slct",
            "--cojo-p",format(cojo_p,scientific=TRUE),"--cojo-wind",as.character(as.integer(cojo_wind_kb)),
            "--cojo-collinear",as.character(cojo_collinear),"--thread-num",as.character(as.integer(cojo_threads)),"--out",sp)
  st <- system2(exe,args=.cater_quote_args(args),stdout=paste0(sp,".stdout"),stderr=paste0(sp,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO selection failed for %s",prefix)
  selected_step <- .cater_read_jma_ref(paste0(sp,".jma.cojo"))
  if (!nrow(selected_step)) return(finish(selected=selected_step,selected_step=selected_step))
  if (nrow(selected_step)==1L) return(finish(selected=selected_step,selected_step=selected_step,
    ld=matrix(1,1,1,dimnames=list(selected_step$snp,selected_step$snp))))
  jp <- paste0(prefix,".joint"); .cater_remove_prefix_outputs(jp)
  args2 <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",paste0(sp,".jma.cojo"),"2","header",
             "--cojo-joint","--thread-num",as.character(as.integer(cojo_threads)),"--output-all","--out",jp)
  st <- system2(exe,args=.cater_quote_args(args2),stdout=paste0(jp,".stdout"),stderr=paste0(jp,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO joint/LD failed for %s",prefix)
  joint_ref <- .cater_read_jma_ref(paste0(jp,".jma.cojo"))
  if (!nrow(joint_ref)) .cater_stop("Manc-COJO joint retained no SNPs for %s",prefix)
  if (length(setdiff(joint_ref$snp,selected_step$snp))) .cater_stop("Manc-COJO joint output contains unselected SNPs for %s",prefix)
  ld <- .cater_read_manc_ldr(paste0(jp,".ldr.cojo"),joint_ref$snp)
  finish(selected=joint_ref,selected_step=selected_step,joint_removed=setdiff(selected_step$snp,joint_ref$snp),ld=ld)
}'''
s = s[:start] + new_run + s[end:]

start = s.index('.cater_run_sibling_cis <- function')
end = s.index('\n\n.cater_build_mvmr <- function', start)
new_sib = r'''.cater_run_sibling_cis <- function(gene,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
                                   cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,
                                   plink_bin="plink",enable_ld_diagnosis=TRUE,
                                   ld_diag_loglr=2,ld_diag_abs_z=2) {
  q<-.cater_read_gene(gene,eqtl_dir,qtl_n);if(is.null(q)) return(NULL)
  reg<-.cater_make_regions(gene,character(),annotation,cis_window,cis_window)
  if(is.null(reg)) return(NULL)
  cmap<-.cater_candidate_map(q,reg);if(!nrow(cmap)) return(list(qtl=q,selected=data.frame()))
  co<-.cater_run_cojo(q,cmap,ld_bfile,manc_cojo_bin,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,
    plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z)
  list(qtl=q,selected=co$selected)
}'''
s = s[:start] + new_sib + s[end:]

old = '''.cater_build_mvmr <- function(target,target_qtl,target_sel,active_siblings,annotation,eqtl_dir,qtl_n,
                              ld_bfile,manc_cojo_bin,cis_window,cojo_p,cojo_wind_kb,cojo_collinear,
                              cojo_threads,outdir,outcome,drop_palindromic,exposure_corr,min_cond_F,
                              mvmr_max_r2,mvmr_max_condition,verbose) {'''
new = '''.cater_build_mvmr <- function(target,target_qtl,target_sel,active_siblings,annotation,eqtl_dir,qtl_n,
                              ld_bfile,manc_cojo_bin,cis_window,cojo_p,cojo_wind_kb,cojo_collinear,
                              cojo_threads,outdir,outcome,drop_palindromic,exposure_corr,min_cond_F,
                              mvmr_max_r2,mvmr_max_condition,verbose,plink_bin="plink",
                              enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2) {'''
s = replace_once(s,old,new,'mvmr signature')
old = '''    x<-tryCatch(.cater_run_sibling_cis(z,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
      cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,pr,verbose),error=function(e)NULL)'''
new = '''    x<-tryCatch(.cater_run_sibling_cis(z,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
      cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,pr,verbose,
      plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,
      ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z),error=function(e)NULL)'''
s = replace_once(s,old,new,'sibling call')

old = '''                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater")) {'''
new = '''                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2) {'''
s = replace_once(s,old,new,'public signature')
old = '''  primary_policy<-match.arg(primary_policy)
  if(!dir.exists(eqtl_dir)) .cater_stop("eqtl_dir does not exist: %s",eqtl_dir)'''
new = '''  primary_policy<-match.arg(primary_policy)
  if(length(plink_bin)!=1L||is.na(plink_bin)||!nzchar(as.character(plink_bin))) .cater_stop("plink_bin must be a non-empty scalar")
  if(length(enable_ld_diagnosis)!=1L||is.na(enable_ld_diagnosis)) .cater_stop("enable_ld_diagnosis must be TRUE or FALSE")
  if(length(ld_diag_loglr)!=1L||!is.finite(ld_diag_loglr)||ld_diag_loglr<0) .cater_stop("ld_diag_loglr must be a finite non-negative scalar")
  if(length(ld_diag_abs_z)!=1L||!is.finite(ld_diag_abs_z)||ld_diag_abs_z<0) .cater_stop("ld_diag_abs_z must be a finite non-negative scalar")
  if(!dir.exists(eqtl_dir)) .cater_stop("eqtl_dir does not exist: %s",eqtl_dir)'''
s = replace_once(s,old,new,'public validation')

old = '''    n_parent_tf=0,n_candidate_snp=0,n_cojo_signal=0,n_cis_signal=0,n_trans_signal=0,
    n_ld_aligned_signal=0,n_mr_iv=0,n_mr_cis_iv=0,n_mr_trans_iv=0,'''
new = '''    n_parent_tf=0,n_candidate_snp=0,n_ld_diag_removed=0,n_ld_diag_susie=0,
    n_ld_reference_missing=0,n_ld_allele_mismatch=0,n_ld_nonfinite=0,ld_diag_lambda_max=NA_real_,
    n_cojo_signal=0,n_cis_signal=0,n_trans_signal=0,n_ld_aligned_signal=0,
    n_mr_iv=0,n_mr_cis_iv=0,n_mr_trans_iv=0,'''
s = replace_once(s,old,new,'summary fields')

old = '''    co<-tryCatch(.cater_run_cojo(qtl,cmap,ld_bfile,manc_cojo_bin,cojo_p,cojo_wind_kb,
      cojo_collinear,cojo_threads,file.path(outdir,"cojo",target),verbose),error=function(e)e)
    if(inherits(co,"error")){warning(sprintf("[%s] %s",target,conditionMessage(co)));srow$status<-"COJO_FAILED";summ[[target]]<-srow;next}
    if(!nrow(co$selected)){srow$status<-"NO_COJO_SIGNAL";summ[[target]]<-srow;next}'''
new = '''    co<-tryCatch(.cater_run_cojo(qtl,cmap,ld_bfile,manc_cojo_bin,cojo_p,cojo_wind_kb,
      cojo_collinear,cojo_threads,file.path(outdir,"cojo",target),verbose,
      plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,
      ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z),error=function(e)e)
    if(inherits(co,"error")){
      warning(sprintf("[%s] %s",target,conditionMessage(co)))
      srow$status<-if(grepl("SuSiE-RSS|LD diagnosis|PLINK",conditionMessage(co))) "LD_DIAGNOSIS_FAILED" else "COJO_FAILED"
      summ[[target]]<-srow;next
    }
    if(nrow(co$ld_diagnosis)){
      dd<-co$ld_diagnosis
      srow$n_ld_diag_removed<-sum(dd$remove,na.rm=TRUE)
      srow$n_ld_diag_susie<-sum(dd$status=="SUSIE_LD_INCONSISTENT",na.rm=TRUE)
      srow$n_ld_reference_missing<-sum(dd$status=="LD_REFERENCE_MISSING",na.rm=TRUE)
      srow$n_ld_allele_mismatch<-sum(dd$status=="LD_ALLELE_MISMATCH",na.rm=TRUE)
      srow$n_ld_nonfinite<-sum(dd$status=="LD_NONFINITE",na.rm=TRUE)
      if(any(is.finite(dd$lambda))) srow$ld_diag_lambda_max<-max(dd$lambda[is.finite(dd$lambda)])
    }
    if(!nrow(co$selected)){
      srow$status<-if(nrow(co$ld_diagnosis)&&all(co$ld_diagnosis$remove)) "NO_CANDIDATE_AFTER_LD_DIAGNOSIS" else "NO_COJO_SIGNAL"
      summ[[target]]<-srow;next
    }'''
s = replace_once(s,old,new,'main cojo call')

old = '''      net<-.cater_build_mvmr(target,qtl,co$selected,sib$active,ann,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
        cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,outdir,outcome,drop_palindromic,
        exposure_corr,min_cond_F,mvmr_max_r2,mvmr_max_condition,verbose)'''
new = '''      net<-.cater_build_mvmr(target,qtl,co$selected,sib$active,ann,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
        cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,outdir,outcome,drop_palindromic,
        exposure_corr,min_cond_F,mvmr_max_r2,mvmr_max_condition,verbose,
        plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,
        ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z)'''
s = replace_once(s,old,new,'main mvmr call')

core.write_text(s)

test = Path('tests/v06_hardening_smoke.R')
t = test.read_text()
old_test = '''# 12. Backward compatibility: new options must not shift established trailing positional arguments.
formal_names <- names(formals(cater_mr))
stopifnot(identical(tail(formal_names, 4L), c("outdir","drop_palindromic","verbose","primary_policy")))'''
new_test = '''# 12. Backward compatibility: appending new options must not shift established positional arguments.
formal_names <- names(formals(cater_mr))
legacy_tail <- c("outdir","drop_palindromic","verbose","primary_policy")
i_outdir <- match("outdir",formal_names)
stopifnot(identical(formal_names[i_outdir:(i_outdir+3L)],legacy_tail))
stopifnot(all(match(c("plink_bin","enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z"),formal_names) > match("primary_policy",formal_names)))'''
t = replace_once(t,old_test,new_test,'test 12')
marker = 'cat("CATER-MR direct-core evidence-safety tests passed\\n")'
if marker not in t:
    raise RuntimeError('test end marker not found')
extra = r'''

# 14. Pre-COJO z/LD alignment handles direct swaps and strand complements.
qz <- data.frame(snp=paste0("rs",1:5),a1=c("A","A","A","A","A"),a2=c("G","G","C","C","C"),
                 beta=rep(.2,5),se=rep(.1,5),stringsAsFactors=FALSE)
rz <- data.frame(snp=paste0("rs",1:5),ld_a1=c("A","G","T","G","A"),ld_a2=c("G","A","G","T","G"),stringsAsFactors=FALSE)
az <- .cater_align_z_to_plink(qz,rz)
stopifnot(identical(az$alignment,c("same","swap","strand_same","strand_swap","mismatch")))
stopifnot(isTRUE(all.equal(az$z[1:4],c(2,-2,2,-2),tolerance=1e-12)),is.na(az$z[5]),!az$allele_match[5])

# 15. mapgen/SuSiE-RSS detection rule is strict: logLR > 2 and |z| > 2.
cd <- data.frame(z=c(2.1,2,5,-3),logLR=c(2.1,3,1,2.2))
stopifnot(identical(.cater_ld_outlier_index(cd,2,2),c(1L,4L)))

# 16. The public API exposes LD diagnosis controls while retaining legacy argument positions.
stopifnot(all(c("plink_bin","enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z") %in% names(formals(cater_mr))))
stopifnot(identical(formals(cater_mr)$enable_ld_diagnosis,TRUE))
'''
t = t.replace(marker, extra + '\n' + marker, 1)
test.write_text(t)
