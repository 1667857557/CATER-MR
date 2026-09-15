source("CATER_MR.R")

# Annotation contract.
tx <- data.frame(symbol="X",chr="1",txStart=100,stringsAsFactors=FALSE)
stopifnot(inherits(try(.cater_standardize_annotation(tx),silent=TRUE),"try-error"))
amb <- data.frame(symbol=c("X","X"),chr=c("1","1"),tss=c(100,900))
stopifnot(inherits(try(.cater_standardize_annotation(amb),silent=TRUE),"try-error"))

# MVMR numerical condition diagnostic is invariant to exposure-unit rescaling.
B <- rbind(c(.20,.00),c(.10,.15),c(.00,.20),c(.12,.05)); colnames(B) <- c("X","Z")
theta <- c(X=.5,Z=.6); by <- as.numeric(B %*% theta)
seB <- matrix(.02,nrow(B),ncol(B),dimnames=dimnames(B))
Cexp <- diag(2); dimnames(Cexp) <- list(c("X","Z"),c("X","Z"))
m1 <- .cater_mvmr_fit(B,seB,by,rep(.02,nrow(B)),diag(nrow(B)),c("X","Z"),Cexp)
B2 <- B; SE2 <- seB; B2[,"Z"] <- B2[,"Z"]*1e4; SE2[,"Z"] <- SE2[,"Z"]*1e4
m2 <- .cater_mvmr_fit(B2,SE2,by,rep(.02,nrow(B)),diag(nrow(B)),c("X","Z"),Cexp)
stopifnot(m1$status=="OK",m2$status=="OK",abs(m1$condition-m2$condition)<1e-8)

# Leave-one-TF removes every IV annotated to the TF.
h <- data.frame(snp=c("c","a","b","d"),source=c("cis","trans","trans","trans"),
                parent_tf=c("","T1","T1;T2","T2"),locus_id=c("cis","TF:T1","TF:T1+T2","TF:T2"),
                bx=c(.2,.1,.08,.09),bx_se=.02,by=c(.1,.05,.04,.045),by_se=.03)
Rh <- diag(4); dimnames(Rh) <- list(h$snp,h$snp)
ldg <- .cater_locus_diagnostics(h,Rh)
t1 <- ldg$table[ldg$table$tf=="T1",]
stopifnot(t1$n_removed_iv==2L,all(c("a","b") %in% strsplit(t1$removed_snps,";",fixed=TRUE)[[1]]))

# Total-effect safeguard excludes descendants from sibling candidates.
grn <- data.frame(TF=c("T","T","X","M"),Target=c("X","Z","M","Z"))
td <- tempfile("cater_desc_"); dir.create(td)
zss <- data.frame(SNP="g1",CHR=1,BP=100,A1="A",A2="G",EAF=.2,BETA=.2,SE=.02,P=1e-20,N=500)
con <- gzfile(file.path(td,"Z.txt.gz"),"wt"); utils::write.table(zss,con,sep="\t",quote=FALSE,row.names=FALSE); close(con)
hs <- data.frame(snp="g1",source="trans",parent_tf="T",locus_id="TF:T",bx=.2,bx_se=.02,by=.1,by_se=.02)
ref <- data.frame(snp="g1",ld_a1="A",ld_a2="G"); L <- matrix(1,1,1,dimnames=list("g1","g1"))
sib <- .cater_sibling_screen("X",hs,grn,ref,L,td,500,.05,preserve_total_effect=TRUE)
stopifnot(sib$n_candidate==0L,sib$status=="NO_TESTABLE_SIBLINGS",!sib$complete,!sib$testable)
unlink(td,recursive=TRUE)

# Primary semantics.
mkfit <- function(status,b=.2,se=.05,p=.01,n=2L,info=10,prec=12) data.frame(n_iv=n,beta=b,se=se,p=p,Q=NA,Q_p=NA,information=info,joint_wald_per_df=NA,effective_F=NA,precision_information=prec,mean_F=NA,min_F=NA,status=status)
fits <- list(cis=mkfit("NO_IV",n=0L,info=NA,prec=NA),trans=mkfit("OK"),combined=mkfit("OK"))
d0 <- .cater_primary_decision(fits,has_trans=TRUE,sibling_screen_performed=TRUE,sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE)
stopifnot(is.na(d0$model),d0$status=="TRANS_ONLY_SENSITIVITY")
fits2 <- list(cis=mkfit("OK",b=.11,se=.04,p=.01),trans=mkfit("OK"),combined=mkfit("LD_SINGULAR"))
d2 <- .cater_primary_decision(fits2,has_trans=TRUE,sibling_screen_performed=TRUE,sibling_screen_complete=TRUE,n_active_siblings=0L,sibling_screen_testable=TRUE,primary_policy="screened_cater")
stopifnot(d2$model=="cis",d2$status=="CATER_COMBINED_FAILED_CIS_FALLBACK")

# Signed LD orientation.
Ro <- matrix(c(1,.4,.4,1),2,2,dimnames=list(c("a","b"),c("a","b")))
from <- data.frame(snp=c("a","b"),ld_a1=c("A","C"),ld_a2=c("G","T"))
to <- data.frame(snp=c("a","b"),ld_a1=c("G","C"),ld_a2=c("A","T"))
Rr <- .cater_reorient_ld(Ro,from,to); stopifnot(abs(Rr["a","b"]+.4)<1e-12)

# v0.8 strict input contract.
good <- list(grn_build="GRCh38",eqtl_build="GRCh38",ld_build="GRCh38",eqtl_n_unit="donors",
             eqtl_ancestry="EUR",ld_ancestry="EUR",cis_full_summary=TRUE,trans_data_mode="significant_only")
stopifnot(is.list(.cater_validate_input_manifest(good,TRUE)))
badN <- good; badN$eqtl_n_unit <- "cells"
stopifnot(inherits(try(.cater_validate_input_manifest(badN,TRUE),silent=TRUE),"try-error"))
badB <- good; badB$ld_build <- "GRCh37"
stopifnot(inherits(try(.cater_validate_input_manifest(badB,TRUE),silent=TRUE),"try-error"))

# Physical-locus merging and cis absorption.
ann_locus <- data.frame(symbol=c("X","A","B","C"),chr=rep("1",4),tss=c(1000,10000,11500,20000))
reg_locus <- .cater_make_regions("X",c("A","B","C"),ann_locus,cis_window=50,tf_window=1000)
stopifnot(reg_locus$locus_id[reg_locus$gene=="A"]=="TF:A+B",reg_locus$locus_id[reg_locus$gene=="B"]=="TF:A+B")
ann_cis <- data.frame(symbol=c("X","A","B","C"),chr=rep("1",4),tss=c(10000,11500,13000,20000))
reg_cis <- .cater_make_regions("X",c("A","B","C"),ann_cis,cis_window=1000,tf_window=1000)
stopifnot(all(reg_cis$locus_id[reg_cis$gene %in% c("X","A","B")]=="cis"),reg_cis$locus_id[reg_cis$gene=="C"]=="TF:C")

# Default QC.
q_qc <- data.frame(snp=c("mhc","pal","ok"),chr=c("6","1","1"),pos=c(30000000,100,200),
                   a1=c("A","A","A"),a2=c("C","T","G"),beta=1,se=1,p=.1,eaf=.2,n=100)
q_keep <- .cater_filter_analysis_variants(q_qc,drop_palindromic=TRUE,exclude_mhc=TRUE)
stopifnot(identical(q_keep$snp,"ok"))

# Obsolete COJO/LD-diagnosis API is absent.
fml <- names(formals(cater_mr))
stopifnot("ld_threads" %in% fml)
stopifnot(!any(c("manc_cojo_bin","cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads",
                 "enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z") %in% fml))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))

cat("CATER-MR core hardening tests passed\n")
