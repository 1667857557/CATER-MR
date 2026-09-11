source("CATER_MR.R")

# 1. Manc-COJO LD parsing must fail closed when a requested SNP is absent.
f <- tempfile(fileext=".ldr.cojo")
writeLines(c("# Chromosome 1","SNP\trs1","rs1\t1"),f)
stopifnot(inherits(try(.cater_read_manc_ldr(f,c("rs1","rs2")),silent=TRUE),"try-error"))
writeLines(c("# Chromosome 1","SNP\trs1","rs1\t1","# Chromosome 2","SNP\trs2","rs2\t1"),f)
R <- .cater_read_manc_ldr(f,c("rs1","rs2"))
stopifnot(R["rs1","rs2"]==0,R["rs1","rs1"]==1,R["rs2","rs2"]==1)
unlink(f)

# 2. Generic annotation requires an explicit TSS and rejects ambiguous gene-level TSS.
tx <- data.frame(symbol="X",chr="1",txStart=100,stringsAsFactors=FALSE)
stopifnot(inherits(try(.cater_standardize_annotation(tx),silent=TRUE),"try-error"))
amb <- data.frame(symbol=c("X","X"),chr=c("1","1"),tss=c(100,900))
stopifnot(inherits(try(.cater_standardize_annotation(amb),silent=TRUE),"try-error"))

# 3. MVMR numerical condition diagnostic is invariant to harmless exposure-unit rescaling.
B <- rbind(c(.20,.00),c(.10,.15),c(.00,.20),c(.12,.05))
colnames(B) <- c("X","Z")
theta <- c(X=.5,Z=.6); by <- as.numeric(B %*% theta)
seB <- matrix(.02,nrow(B),ncol(B),dimnames=dimnames(B))
Cexp <- diag(2); dimnames(Cexp) <- list(c("X","Z"),c("X","Z"))
m1 <- .cater_mvmr_fit(B,seB,by,rep(.02,nrow(B)),diag(nrow(B)),c("X","Z"),Cexp)
B2 <- B; SE2 <- seB; B2[,"Z"] <- B2[,"Z"]*1e4; SE2[,"Z"] <- SE2[,"Z"]*1e4
m2 <- .cater_mvmr_fit(B2,SE2,by,rep(.02,nrow(B)),diag(nrow(B)),c("X","Z"),Cexp)
stopifnot(m1$status=="OK",m2$status=="OK",abs(m1$condition-m2$condition)<1e-8)

# 4. Leave-one-TF removes every trans IV annotated to that TF, including shared-parent SNPs.
h <- data.frame(snp=c("c","a","b","d"),source=c("cis","trans","trans","trans"),
                parent_tf=c("","T1","T1;T2","T2"),locus_id=c("cis","TF:T1","TF:T1;T2","TF:T2"),
                bx=c(.2,.1,.08,.09),bx_se=.02,by=c(.1,.05,.04,.045),by_se=.03)
Rh <- diag(4); dimnames(Rh) <- list(h$snp,h$snp)
ldg <- .cater_locus_diagnostics(h,Rh,Rh)
t1 <- ldg$table[ldg$table$tf=="T1",]
stopifnot(t1$n_removed_iv==2L,all(c("a","b") %in% strsplit(t1$removed_snps,";",fixed=TRUE)[[1]]))

# 5. Total-effect safeguard excludes indirect descendants from sibling adjustment candidates.
grn <- data.frame(TF=c("T","T","X","M"),Target=c("X","Z","M","Z"))
td <- tempfile("cater_desc_");dir.create(td)
zss <- data.frame(SNP="g1",CHR=1,BP=100,A1="A",A2="G",EAF=.2,BETA=.2,SE=.02,P=1e-20,N=500)
con<-gzfile(file.path(td,"Z.txt.gz"),"wt");utils::write.table(zss,con,sep="\t",quote=FALSE,row.names=FALSE);close(con)
hs<-data.frame(snp="g1",source="trans",parent_tf="T",locus_id="TF:T",bx=.2,bx_se=.02,by=.1,by_se=.02)
ref<-data.frame(snp="g1",ld_a1="A",ld_a2="G");L<-matrix(1,1,1,dimnames=list("g1","g1"))
sib <- .cater_sibling_screen("X",hs,grn,ref,L,td,500,.05,preserve_total_effect=TRUE)
stopifnot(sib$n_candidate==0L,sib$status=="NO_TESTABLE_SIBLINGS",!sib$complete,!sib$testable)
unlink(td,recursive=TRUE)

# 6. Empty candidate sibling set is not interpreted as validated exclusion restriction evidence.
grn2 <- data.frame(TF="T",Target="X")
sib2 <- .cater_sibling_screen("X",hs,grn2,ref,L,tempdir(),500,.05)
stopifnot(sib2$status=="NO_TESTABLE_SIBLINGS",!sib2$complete,!sib2$testable)

# 7. Trans-only cannot become primary by default; explicit research opt-in requires multi-locus support.
mkfit <- function(status,b=.2,se=.05,p=.01,n=2L,info=10,prec=12) data.frame(n_iv=n,beta=b,se=se,p=p,Q=NA,Q_p=NA,information=info,joint_wald_per_df=NA,effective_F=NA,precision_information=prec,mean_F=NA,min_F=NA,status=status)
fits <- list(cis=mkfit("NO_IV",n=0L,info=NA,prec=NA),trans=mkfit("OK"),combined=mkfit("OK"))
d0 <- .cater_primary_decision(fits,has_trans=TRUE,sibling_screen_performed=TRUE,sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE)
stopifnot(is.na(d0$model),d0$status=="TRANS_ONLY_SENSITIVITY")
d1 <- .cater_primary_decision(fits,has_trans=TRUE,sibling_screen_performed=TRUE,sibling_screen_complete=TRUE,n_active_siblings=0L,primary_policy="screened_cater",allow_trans_only_primary=TRUE,n_trans_tf_loci=3L,min_trans_tf_loci=3L,sibling_screen_testable=TRUE)
stopifnot(d1$model=="combined",d1$status=="OK_CATER_TRANS_ONLY_EXPLORATORY")

# 8. Exposure-side signal fraction and first-order MR precision fraction are distinct diagnostics.
g2 <- data.frame(snp=c("c","t"),source=c("cis","trans"),parent_tf=c("","T"),locus_id=c("cis","TF:T"),bx=c(.1,.1),bx_se=c(.01,.001),by=c(.05,.05),by_se=c(.001,.1))
R2<-diag(2);dimnames(R2)<-list(g2$snp,g2$snp)
fc <- .cater_givw(g2[1,,drop=FALSE],matrix(1,1,1),matrix(1,1,1));fa <- .cater_givw(g2,R2,R2)
fe <- .cater_trans_exposure_signal_fraction(fc,fa);fp <- .cater_trans_mr_precision_fraction(fc,fa)
stopifnot(fe>.98,fp<.001)

# 9. Signed LD reorientation follows allele flips on rows and columns.
Ro<-matrix(c(1,.4,.4,1),2,2,dimnames=list(c("a","b"),c("a","b")))
from<-data.frame(snp=c("a","b"),ld_a1=c("A","C"),ld_a2=c("G","T"))
to<-data.frame(snp=c("a","b"),ld_a1=c("G","C"),ld_a2=c("A","T"))
Rr<-.cater_reorient_ld(Ro,from,to);stopifnot(abs(Rr["a","b"]+.4)<1e-12)

# 10. Strict biological/statistical input contract rejects cells-as-N, truncated QTL and build mismatch.
good <- list(grn_build="GRCh38",eqtl_build="GRCh38",ld_build="GRCh38",outcome_build="GRCh38",
             eqtl_n_unit="donors",eqtl_full_summary=TRUE,eqtl_ancestry="EUR",outcome_ancestry="EUR",ld_ancestry="EUR",trans_qtl_qc="source-study QC")
stopifnot(is.list(.cater_validate_input_manifest(good,TRUE)))
badN<-good;badN$eqtl_n_unit<-"cells";stopifnot(inherits(try(.cater_validate_input_manifest(badN,TRUE),silent=TRUE),"try-error"))
badF<-good;badF$eqtl_full_summary<-FALSE;stopifnot(inherits(try(.cater_validate_input_manifest(badF,TRUE),silent=TRUE),"try-error"))
badB<-good;badB$ld_build<-"GRCh37";stopifnot(inherits(try(.cater_validate_input_manifest(badB,TRUE),silent=TRUE),"try-error"))

# 11. Direct-core Occam semantics: no qtl-outcome-overlap field and cis is the default primary anchor.
manifest_occam <- .cater_validate_input_manifest(good, TRUE)
stopifnot(!"qtl_outcome_overlap" %in% names(manifest_occam))
fits_occam <- list(cis=mkfit("OK",b=.11,se=.04,p=.01),trans=mkfit("OK"),combined=mkfit("OK",b=.20,se=.03,p=.001))
d_occam <- .cater_primary_decision(fits_occam,has_trans=TRUE,sibling_screen_performed=TRUE,
                                   sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE)
stopifnot(d_occam$model=="cis",d_occam$status=="OK_CIS_ANCHOR_PRIMARY")
stopifnot("primary_policy" %in% names(formals(cater_mr)))

# 12. Backward compatibility: appending new options must not shift established positional arguments.
formal_names <- names(formals(cater_mr))
legacy_tail <- c("outdir","drop_palindromic","verbose","primary_policy")
i_outdir <- match("outdir",formal_names)
stopifnot(identical(formal_names[i_outdir:(i_outdir+3L)],legacy_tail))
stopifnot(all(match(c("plink_bin","enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z"),formal_names) > match("primary_policy",formal_names)))

# 13. screened_cater promotes a valid combined fit, but a failed augmented fit must fall back to cis.
fits_combined_fail <- list(cis=mkfit("OK",b=.11,se=.04,p=.01),
                           trans=mkfit("OK"),combined=mkfit("LD_SINGULAR"))
d_combined_fail <- .cater_primary_decision(
  fits_combined_fail,has_trans=TRUE,sibling_screen_performed=TRUE,
  sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE,
  primary_policy="screened_cater")
stopifnot(d_combined_fail$model=="cis",
          d_combined_fail$status=="CATER_COMBINED_FAILED_CIS_FALLBACK",
          d_combined_fail$evidence_status=="COMBINED_FIT_FAILED")

# 14. Pre-COJO z/LD alignment accepts only direct A1/A2 same/swap coding.
qz <- data.frame(snp=paste0("rs",1:5),a1=c("A","A","A","A","A"),a2=c("G","G","C","C","C"),
                 beta=rep(.2,5),se=rep(.1,5),stringsAsFactors=FALSE)
rz <- data.frame(snp=paste0("rs",1:5),ld_a1=c("A","G","T","G","A"),ld_a2=c("G","A","G","T","G"),stringsAsFactors=FALSE)
az <- .cater_align_z_to_plink(qz,rz)
stopifnot(identical(az$alignment,c("same","swap","strand_same","strand_swap","mismatch")))
stopifnot(isTRUE(all.equal(az$z[1:2],c(2,-2),tolerance=1e-12)))
stopifnot(all(is.na(az$z[3:5])),identical(az$allele_match,c(TRUE,TRUE,FALSE,FALSE,FALSE)))

# 15. mapgen/SuSiE-RSS detection rule is strict: logLR > 2 and |z| > 2.
cd <- data.frame(z=c(2.1,2,5,-3),logLR=c(2.1,3,1,2.2))
stopifnot(identical(.cater_ld_outlier_index(cd,2,2),c(1L,4L)))

# 16. The public API exposes LD diagnosis controls while retaining legacy argument positions.
stopifnot(all(c("plink_bin","enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z") %in% names(formals(cater_mr))))
stopifnot(identical(formals(cater_mr)$enable_ld_diagnosis,TRUE))

# 17. Diagnostic LD validation leaves PSD/eigenvalue handling to susieR.
Rd <- matrix(c(1,.9,.9,.9,1,.6199,.9,.6199,1),3,3,byrow=TRUE)
stopifnot(identical(dim(.cater_validate_diag_ld(Rd)),c(3L,3L)))

# 18. Overlapping TF windows are merged by genomic interval before candidate assignment.
ann_locus <- data.frame(symbol=c("X","A","B","C"),chr=rep("1",4),
                        tss=c(1000,10000,11500,20000),stringsAsFactors=FALSE)
reg_locus <- .cater_make_regions("X",c("A","B","C"),ann_locus,cis_window=50,tf_window=1000)
tr_locus <- reg_locus[reg_locus$type=="trans",,drop=FALSE]
stopifnot(tr_locus$locus_id[tr_locus$gene=="A"]=="TF:A+B",
          tr_locus$locus_id[tr_locus$gene=="B"]=="TF:A+B",
          tr_locus$locus_id[tr_locus$gene=="C"]=="TF:C")
q_locus <- data.frame(snp=paste0("g",1:5),chr=rep("1",5),
                      pos=c(9500,10750,12000,12500,20000),stringsAsFactors=FALSE)
cm_locus <- .cater_candidate_map(q_locus,reg_locus)
stopifnot(identical(cm_locus$parent_tf,c("A","A;B","B","B","C")),
          identical(cm_locus$locus_id,c("TF:A+B","TF:A+B","TF:A+B","TF:A+B","TF:C")))
gg <- .cater_ld_diagnosis_groups(q_locus,cm_locus)
stopifnot(length(unique(gg[1:4]))==1L,gg[5]!=gg[1])

# 19. A TF window connected to the target cis window is absorbed into the full cis locus.
# A directly overlaps target cis; B overlaps A, so the whole connected component is cis.
ann_cis_tf <- data.frame(symbol=c("X","A","B","C"),chr=rep("1",4),
                         tss=c(10000,11500,13000,20000),stringsAsFactors=FALSE)
reg_cis_tf <- .cater_make_regions("X",c("A","B","C"),ann_cis_tf,cis_window=1000,tf_window=1000)
stopifnot(reg_cis_tf$locus_id[reg_cis_tf$gene=="X"]=="cis",
          reg_cis_tf$locus_id[reg_cis_tf$gene=="A"]=="cis",
          reg_cis_tf$locus_id[reg_cis_tf$gene=="B"]=="cis",
          reg_cis_tf$locus_id[reg_cis_tf$gene=="C"]=="TF:C")
q_cis_tf <- data.frame(snp=paste0("ct",1:4),chr=rep("1",4),
                       pos=c(9500,11500,13000,20000),stringsAsFactors=FALSE)
cm_cis_tf <- .cater_candidate_map(q_cis_tf,reg_cis_tf)
stopifnot(identical(cm_cis_tf$source,c("cis","cis","cis","trans")),
          identical(cm_cis_tf$locus_id,c("cis","cis","cis","TF:C")),
          identical(cm_cis_tf$parent_tf,c("","","","C")))
gg_cis_tf <- .cater_ld_diagnosis_groups(q_cis_tf,cm_cis_tf)
stopifnot(length(unique(gg_cis_tf[1:3]))==1L,gg_cis_tf[4]!=gg_cis_tf[1])

# 20. Palindromic A/T and C/G variants are never signed from allele labels alone.
qp <- data.frame(snp=c("p1","p2","n1"),a1=c("A","C","A"),a2=c("T","G","G"),
                 beta=c(.2,.2,.2),se=c(.1,.1,.1),stringsAsFactors=FALSE)
rp <- data.frame(snp=c("p1","p2","n1"),ld_a1=c("T","C","G"),ld_a2=c("A","G","A"),
                 stringsAsFactors=FALSE)
ap <- .cater_align_z_to_plink(qp,rp)
stopifnot(identical(ap$palindromic,c(TRUE,TRUE,FALSE)),
          identical(ap$alignment,c("palindromic_unresolved","palindromic_unresolved","swap")),
          all(is.na(ap$z[1:2])),isTRUE(all.equal(ap$z[3],-2,tolerance=1e-12)),
          identical(ap$allele_match,c(FALSE,FALSE,TRUE)))

# 21. Every variant whose original LD row contains a non-finite entry is removed together.
Rnf <- matrix(c(1,NaN,0.2,NaN,1,0.3,0.2,0.3,1),3,3,byrow=TRUE,
              dimnames=list(c("a","b","c"),c("a","b","c")))
fnf <- .cater_filter_nonfinite_ld_rows(Rnf)
stopifnot(identical(sort(fnf$nonfinite),c("a","b")),
          identical(dim(fnf$R),c(1L,1L)),rownames(fnf$R)=="c")
Rnf2 <- matrix(c(1,NaN,NaN,1),2,2,byrow=TRUE,
               dimnames=list(c("x","y"),c("x","y")))
fnf2 <- .cater_filter_nonfinite_ld_rows(Rnf2)
stopifnot(setequal(fnf2$nonfinite,c("x","y")),nrow(fnf2$R)==0L)

# 22. The diagnosis switch must be a scalar logical; numeric/string truthy values fail closed.
for (bad_flag in list(1,"TRUE")) {
  err <- tryCatch({
    cater_mr(grn=data.frame(TF="A",Target="B"),eqtl_dir=tempdir(),outcome=data.frame(),
             ld_bfile="unused",enable_ld_diagnosis=bad_flag)
    NULL
  },error=function(e)e)
  stopifnot(inherits(err,"error"),grepl("enable_ld_diagnosis must be TRUE or FALSE",conditionMessage(err),fixed=TRUE))
}

# 23. Default analysis QC uses GRCh38/hg38 extended MHC chr6:25-36 Mb and drops palindromes.
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

cat("CATER-MR direct-core evidence-safety tests passed\n")