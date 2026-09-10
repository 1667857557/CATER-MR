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

cat("CATER-MR direct-core evidence-safety tests passed\n")
