from pathlib import Path

p=Path('R/CATER_MR_core_v05.R')
s=p.read_text()
s=s.replace('.CATER_VERSION <- "0.7.2"','.CATER_VERSION <- "0.7.3"',1)

anchor='''  out\n}\n\n\n.cater_standardize_grn <- function(grn) {\n'''
helper='''  out\n}\n\n# Conservative GRCh38/hg38 extended MHC exclusion used for genetic analyses.\n# The interval is inclusive and intentionally broader than the core GRC MHC\n# alternate-locus placement to avoid the long-range LD structure around MHC.\n.CATER_MHC_HG38 <- c(chr="6", start=25000000, end=36000000)\n\n.cater_is_mhc_hg38 <- function(chr,pos) {\n  ch <- sub("^chr","",as.character(chr),ignore.case=TRUE)\n  pp <- suppressWarnings(as.numeric(pos))\n  !is.na(ch) & ch==.CATER_MHC_HG38[["chr"]] & is.finite(pp) &\n    pp>=as.numeric(.CATER_MHC_HG38[["start"]]) & pp<=as.numeric(.CATER_MHC_HG38[["end"]])\n}\n\n.cater_filter_analysis_variants <- function(qtl, drop_palindromic=TRUE, exclude_mhc=TRUE) {\n  if (!nrow(qtl)) return(qtl)\n  mhc <- if (isTRUE(exclude_mhc)) .cater_is_mhc_hg38(qtl$chr,qtl$pos) else rep(FALSE,nrow(qtl))\n  pal <- if (isTRUE(drop_palindromic)) .cater_is_palindromic(qtl$a1,qtl$a2) else rep(FALSE,nrow(qtl))\n  out <- qtl[!(mhc|pal),,drop=FALSE]\n  rownames(out) <- NULL\n  attr(out,"analysis_qc") <- list(n_mhc_excluded=sum(mhc),n_palindromic_excluded=sum(pal & !mhc))\n  out\n}\n\n\n.cater_standardize_grn <- function(grn) {\n'''
assert anchor in s
s=s.replace(anchor,helper,1)

# Public API: append one flag only, preserving all existing positional argument meanings.
old='''                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",\n                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2) {\n'''
new='''                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",\n                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2,\n                     exclude_mhc=TRUE) {\n'''
assert old in s
s=s.replace(old,new,1)
old='''  if(!is.logical(enable_ld_diagnosis)||length(enable_ld_diagnosis)!=1L||is.na(enable_ld_diagnosis)) .cater_stop("enable_ld_diagnosis must be TRUE or FALSE")\n'''
new=old+'''  if(!is.logical(exclude_mhc)||length(exclude_mhc)!=1L||is.na(exclude_mhc)) .cater_stop("exclude_mhc must be TRUE or FALSE")\n  if(!is.logical(drop_palindromic)||length(drop_palindromic)!=1L||is.na(drop_palindromic)) .cater_stop("drop_palindromic must be TRUE or FALSE")\n'''
s=s.replace(old,new,1)

# Summary audit fields.
old='''    n_parent_tf=0,n_candidate_snp=0,n_ld_diag_removed=0,n_ld_diag_susie=0,\n'''
new='''    n_parent_tf=0,n_parent_tf_mhc_excluded=0,n_candidate_snp=0,n_mhc_excluded=0,n_palindromic_excluded=0,\n    n_ld_diag_removed=0,n_ld_diag_susie=0,\n'''
assert old in s
s=s.replace(old,new,1)

# Main target: reject MHC targets, remove MHC parents, and prefilter variants before candidate mapping.
old='''    parents<-unique(grn$TF[grn$Target==target]);srow$n_parent_tf<-length(parents)\n    qtl<-.cater_read_gene(target,eqtl_dir,qtl_n)\n    if(is.null(qtl)){srow$status<-"NO_EQTL_FILE";summ[[target]]<-srow;next}\n    regions<-.cater_make_regions(target,parents,ann,cis_window,tf_window)\n    if(is.null(regions)){srow$status<-"NO_TARGET_ANNOTATION";summ[[target]]<-srow;next}\n    cmap<-.cater_candidate_map(qtl,regions);srow$n_candidate_snp<-nrow(cmap)\n'''
new='''    parents<-unique(grn$TF[grn$Target==target])\n    ta<-ann[ann$symbol==target,,drop=FALSE]\n    if(!nrow(ta)){srow$status<-"NO_TARGET_ANNOTATION";summ[[target]]<-srow;next}\n    if(isTRUE(exclude_mhc)&&.cater_is_mhc_hg38(ta$chr[1],ta$tss[1])){\n      srow$status<-"TARGET_IN_MHC_EXCLUDED";summ[[target]]<-srow;next\n    }\n    if(length(parents)&&isTRUE(exclude_mhc)){\n      pa<-ann[match(parents,ann$symbol),,drop=FALSE]\n      parent_mhc<-!is.na(pa$symbol)&.cater_is_mhc_hg38(pa$chr,pa$tss)\n      srow$n_parent_tf_mhc_excluded<-sum(parent_mhc)\n      parents<-parents[!parent_mhc]\n    }\n    srow$n_parent_tf<-length(parents)\n    qtl<-.cater_read_gene(target,eqtl_dir,qtl_n)\n    if(is.null(qtl)){srow$status<-"NO_EQTL_FILE";summ[[target]]<-srow;next}\n    qtl<-.cater_filter_analysis_variants(qtl,drop_palindromic=drop_palindromic,exclude_mhc=exclude_mhc)\n    aqc<-attr(qtl,"analysis_qc"); if(!is.null(aqc)){srow$n_mhc_excluded<-aqc$n_mhc_excluded;srow$n_palindromic_excluded<-aqc$n_palindromic_excluded}\n    regions<-.cater_make_regions(target,parents,ann,cis_window,tf_window)\n    cmap<-.cater_candidate_map(qtl,regions);srow$n_candidate_snp<-nrow(cmap)\n'''
assert old in s
s=s.replace(old,new,1)

# Sibling cis runner accepts and applies the same genetic QC, and excludes MHC genes.
old='''.cater_run_sibling_cis <- function(gene,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,\n                                   cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,\n                                   plink_bin="plink",enable_ld_diagnosis=TRUE,\n                                   ld_diag_loglr=2,ld_diag_abs_z=2) {\n  q<-.cater_read_gene(gene,eqtl_dir,qtl_n);if(is.null(q)) return(NULL)\n  reg<-.cater_make_regions(gene,character(),annotation,cis_window,cis_window)\n'''
new='''.cater_run_sibling_cis <- function(gene,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,\n                                   cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,\n                                   plink_bin="plink",enable_ld_diagnosis=TRUE,\n                                   ld_diag_loglr=2,ld_diag_abs_z=2,drop_palindromic=TRUE,exclude_mhc=TRUE) {\n  ga<-annotation[annotation$symbol==gene,,drop=FALSE]\n  if(isTRUE(exclude_mhc)&&nrow(ga)&&.cater_is_mhc_hg38(ga$chr[1],ga$tss[1])) return(NULL)\n  q<-.cater_read_gene(gene,eqtl_dir,qtl_n);if(is.null(q)) return(NULL)\n  q<-.cater_filter_analysis_variants(q,drop_palindromic=drop_palindromic,exclude_mhc=exclude_mhc)\n  reg<-.cater_make_regions(gene,character(),annotation,cis_window,cis_window)\n'''
assert old in s
s=s.replace(old,new,1)

# MVMR builder propagates the QC settings to sibling cis selection.
old='''                              mvmr_max_r2,mvmr_max_condition,verbose,plink_bin="plink",\n                              enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2) {\n'''
new='''                              mvmr_max_r2,mvmr_max_condition,verbose,plink_bin="plink",\n                              enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2,exclude_mhc=TRUE) {\n'''
assert old in s
s=s.replace(old,new,1)
old='''      plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,\n      ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z),error=function(e)NULL)\n'''
new='''      plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,\n      ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z,\n      drop_palindromic=drop_palindromic,exclude_mhc=exclude_mhc),error=function(e)NULL)\n'''
# only replace first occurrence after sibling runner call inside MVMR; locate via function region
idx=s.index('.cater_build_mvmr <- function')
pre=s[:idx]; tail=s[idx:]
assert old in tail
tail=tail.replace(old,new,1); s=pre+tail

# Main MVMR call propagates exclude_mhc.
old='''        plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,\n        ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z)\n'''
new='''        plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,\n        ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z,exclude_mhc=exclude_mhc)\n'''
# replace last/main occurrence
pos=s.rfind(old); assert pos!=-1
s=s[:pos]+new+s[pos+len(old):]
p.write_text(s)

# Documentation.
p=Path('LD_DIAGNOSIS.md'); d=p.read_text()
d=d.replace('CATER-MR performs an LD/summary-statistic consistency gate before every Manc-COJO selection step by default.',
'''CATER-MR performs an LD/summary-statistic consistency gate before every Manc-COJO selection step by default. Before locus construction, the default genetic QC also removes palindromic A/T and C/G variants (`drop_palindromic=TRUE`) and excludes the extended MHC using GRCh38/hg38 coordinates `chr6:25,000,000-36,000,000` (`exclude_mhc=TRUE`). Targets whose TSS lies inside this MHC interval are skipped, and parent TFs whose TSS lies inside it do not contribute trans loci.''',1)
d=d.replace('Palindromic A/T and C/G variants are not assigned a signed z score from allele labels alone; without frequency-based strand disambiguation they are excluded from SuSiE-RSS diagnosis, retained for downstream COJO, and reported as `NOT_DIAGNOSABLE_PALINDROMIC`.',
'''With the default `drop_palindromic=TRUE`, palindromic A/T and C/G variants are removed before candidate-locus LD diagnosis and COJO. If that option is explicitly disabled, the LD-diagnosis helper still refuses to infer their signed z score from allele labels alone and reports them as `NOT_DIAGNOSABLE_PALINDROMIC`.''',1)
d=d.replace('Set `enable_ld_diagnosis=FALSE` only to reproduce the legacy COJO path without this QC gate.',
'''Set `enable_ld_diagnosis=FALSE` only to reproduce the legacy COJO path without this QC gate. `exclude_mhc=TRUE` is independent of the SuSiE gate and remains enabled by default; its coordinates are fixed to the conservative GRCh38/hg38 extended-MHC interval chr6:25-36 Mb.''',1)
p.write_text(d)

# Regression tests.
p=Path('tests/v06_hardening_smoke.R'); t=p.read_text()
footer='cat("CATER-MR direct-core evidence-safety tests passed\\n")'
assert footer in t
extra=r'''# 23. Default analysis QC uses GRCh38/hg38 extended MHC chr6:25-36 Mb and drops palindromes.
stopifnot(identical(.CATER_MHC_HG38,c(chr="6",start="25000000",end="36000000")) ||
          (as.character(.CATER_MHC_HG38[["chr"]])=="6" &&
           as.numeric(.CATER_MHC_HG38[["start"]])==25000000 &&
           as.numeric(.CATER_MHC_HG38[["end"]])==36000000))
stopifnot(identical(.cater_is_mhc_hg38(c("6","chr6","6","6"),
                                       c(25000000,36000000,24999999,36000001)),
                    c(TRUE,TRUE,FALSE,FALSE)))
q_qc <- data.frame(snp=c("mhc","pal","ok"),chr=c("6","1","1"),pos=c(30000000,100,200),
                   a1=c("A","A","A"),a2=c("C","T","G"),beta=1,se=1,p=.1,eaf=.2,n=100,
                   stringsAsFactors=FALSE)
q_keep <- .cater_filter_analysis_variants(q_qc,drop_palindromic=TRUE,exclude_mhc=TRUE)
stopifnot(identical(q_keep$snp,"ok"),attr(q_keep,"analysis_qc")$n_mhc_excluded==1L,
          attr(q_keep,"analysis_qc")$n_palindromic_excluded==1L)
q_keep_pal <- .cater_filter_analysis_variants(q_qc,drop_palindromic=FALSE,exclude_mhc=TRUE)
stopifnot(identical(q_keep_pal$snp,c("pal","ok")))

'''
t=t.replace(footer,extra+footer,1); p.write_text(t)
