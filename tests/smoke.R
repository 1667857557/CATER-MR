source("CATER_MR.R")
source("SCMORE_GRN.R")

# Basic summary-stat normalization.
ss <- data.frame(
  SNP=c("rs1","rs2"), CHR=c(1,1), BP=c(100,200),
  A1=c("A","C"), A2=c("G","T"), EAF=c(.2,.3),
  BETA=c(.2,-.1), SE=c(.05,.04), N=c(500,500)
)
s <- .cater_standardize_sumstats(ss)
stopifnot(nrow(s)==2L, all(is.finite(s$p)))

# Duplicate/multiallelic IDs are removed entirely.
dup_ss <- rbind(ss[1,], transform(ss[1,], A1="T", A2="G"), ss[2,])
dup_std <- .cater_standardize_sumstats(dup_ss)
stopifnot(identical(dup_std$snp,"rs2"))

# One-hop region assignment.
ann <- data.frame(symbol=c("X","TF1","TF2"), chr=c("1","1","2"), tss=c(1000,5000,10000))
regions <- .cater_make_regions("X",c("TF1","TF2"),ann,100,100)
q <- data.frame(snp=c("c","t1","t2"),chr=c("1","1","2"),pos=c(1000,5000,10000),
                a1="A",a2="G",beta=.1,se=.02,p=1e-9,eaf=.2,n=500)
cm <- .cater_candidate_map(q,regions)
stopifnot(identical(cm$source,c("cis","trans","trans")))

# Manc-COJO chromosome-block LD parser.
f <- tempfile(fileext=".ldr.cojo")
writeLines(c("# Chromosome 1","SNP\trs1\trs2","rs1\t1\t0.25","rs2\t0.25\t1",
             "# Chromosome 2","SNP\trs3","rs3\t1"),f)
ld <- .cater_read_manc_ldr(f,c("rs1","rs2","rs3"))
stopifnot(abs(ld["rs1","rs2"]-.25)<1e-12,ld["rs1","rs3"]==0)
unlink(f)

# Allele-to-LD orientation flips beta and EAF when required.
d <- data.frame(snp=c("a","b"),chr="1",pos=1:2,a1=c("A","G"),a2=c("G","A"),
                beta=c(.2,.3),se=.02,p=1e-10,eaf=c(.2,.3),n=500)
ref <- data.frame(snp=c("a","b"),ld_a1=c("A","A"),ld_a2=c("G","G"))
da <- .cater_align_to_ld(d,ref)
stopifnot(abs(da$beta[1]-.2)<1e-12,abs(da$beta[2]+.3)<1e-12,abs(da$eaf[2]-.7)<1e-12)

# GIVW with identity LD equals ordinary IVW; effective_F equals mean single-SNP F.
g <- data.frame(snp=c("a","b"),source=c("cis","trans"),parent_tf=c("","TF1"),
                locus_id=c("cis","TF:TF1"),bx=c(.2,.1),bx_se=c(.02,.02),
                by=c(.1,.06),by_se=c(.02,.03))
fit <- .cater_givw(g,diag(2))
w <- 1/g$by_se^2
manual <- sum(w*g$bx*g$by)/sum(w*g$bx^2)
stopifnot(abs(fit$beta-manual)<1e-12)
stopifnot(abs(fit$effective_F-mean((g$bx/g$bx_se)^2))<1e-10)

# Same-IV omnibus reduces to squared Wald z for one SNP.
om <- .cater_omnibus(.12,.03,matrix(1,1,1))
stopifnot(abs(om["Q"]-16)<1e-10,om["df"]==1)

# Cis/trans heterogeneity is zero when both ratios are identical.
h <- g; h$by <- .5*h$bx
R2 <- diag(2); dimnames(R2) <- list(h$snp,h$snp)
ht <- .cater_cis_trans_het(h,R2)
stopifnot(abs(ht["z"])<1e-10,abs(ht["p"]-1)<1e-10)

# Trans-only analyses have 100% trans information rather than NA.
empty_cis <- .cater_givw(g[FALSE,,drop=FALSE],matrix(numeric(),0,0))
trans_fit <- .cater_givw(g[g$source=="trans",,drop=FALSE],matrix(1,1,1))
stopifnot(abs(.cater_trans_information_fraction(empty_cis,trans_fit)-1)<1e-12)

# Conditional F uses covariance-weighted nuisance regression, cross-SNP LD, and residual df m-p+1.
Bc <- cbind(X=c(.20,.10,.40,.05), Z=c(.10,.20,.10,.30))
SEc <- cbind(X=c(.01,.10,.01,.20), Z=rep(.05,4))
Cexp <- diag(2); dimnames(Cexp) <- list(c("X","Z"),c("X","Z"))
Rc <- diag(4)
cf_ind <- .cater_conditional_f(Bc,SEc,Rc,c("X","Z"),Cexp)
ols_delta <- as.numeric(stats::lm.fit(Bc[,"Z",drop=FALSE],Bc[,"X"])$coefficients)
cf_delta <- attr(cf_ind,"delta")$X
stopifnot(isTRUE(attr(cf_ind,"converged")["X"]),attr(cf_ind,"df")["X"]==3L)
stopifnot(abs(cf_delta-ols_delta)>0.1)
r <- Bc[,"X"]-Bc[,"Z"]*cf_delta
qv <- c(1,-cf_delta)
V <- .cater_residual_exposure_cov(SEc,Rc,qv,Cexp)
Qmanual <- as.numeric(crossprod(r,solve(V,r)))
stopifnot(abs(cf_ind["X"]-Qmanual/3)<1e-8)
Rc_ld <- matrix(c(1,.2,0,0,.2,1,.1,0,0,.1,1,.3,0,0,.3,1),4,4,byrow=TRUE)
cf_ld <- .cater_conditional_f(Bc,SEc,Rc_ld,c("X","Z"),Cexp)
stopifnot(isTRUE(attr(cf_ld,"converged")["X"]),abs(cf_ld["X"]-cf_ind["X"])>1e-4)

# Partial same-IV sibling coverage is explicitly incomplete even when the tested subset is significant.
td <- tempfile("cater_sib_"); dir.create(td)
zss <- data.frame(SNP="g1",CHR=1,BP=100,A1="A",A2="G",EAF=.2,BETA=.20,SE=.02,P=1e-20,N=500)
con <- gzfile(file.path(td,"Z.txt.gz"),"wt")
utils::write.table(zss,con,sep="\t",quote=FALSE,row.names=FALSE); close(con)
hs <- data.frame(snp=c("g1","g2"),source="trans",parent_tf="TF1",locus_id="TF:TF1")
gs <- data.frame(TF=c("TF1","TF1"),Target=c("X","Z"))
rs <- data.frame(snp=c("g1","g2"),ld_a1="A",ld_a2="G")
Ls <- diag(2); dimnames(Ls) <- list(c("g1","g2"),c("g1","g2"))
sib <- .cater_sibling_screen("X",hs,gs,rs,Ls,td,500,0.05)
stopifnot(sib$n_candidate==1L,sib$n_incomplete==1L,sib$n_iv_missing==1L,!sib$complete)
stopifnot(sib$table$n_iv_requested==2L,sib$table$n_iv_tested==1L,
          sib$table$status=="PARTIAL_SNP_COVERAGE",isTRUE(sib$table$active))
stopifnot("co_perturbation_detected" %in% names(sib$table),
          sib$table$interpretation=="DETECTED_SIBLING_COPERTURBATION")
unlink(td,recursive=TRUE)

# Unresolved or incompletely screened trans estimates never populate primary_*.
mkfit <- function(status,beta=.2,se=.05,p=.01,n_iv=2L,information=10)
  data.frame(n_iv=n_iv,beta=beta,se=se,p=p,Q=NA,Q_p=NA,information=information,
             effective_F=NA,mean_F=NA,min_F=NA,status=status)
fits0 <- list(cis=mkfit("NO_IV",n_iv=0L,information=NA),trans=mkfit("OK"),combined=mkfit("OK"))
dec0 <- .cater_primary_decision(fits0,has_trans=TRUE,sibling_screen_performed=TRUE,
                                sibling_screen_complete=TRUE,n_active_siblings=1L)
stopifnot(is.na(dec0$model),is.na(dec0$p),dec0$status=="TRANS_PLEIOTROPY_UNRESOLVED_NO_CIS")
dec1 <- .cater_primary_decision(fits0,has_trans=TRUE,sibling_screen_performed=TRUE,
                                sibling_screen_complete=FALSE,n_active_siblings=0L)
stopifnot(is.na(dec1$model),is.na(dec1$p),dec1$status=="SIBLING_SCREEN_INCOMPLETE_NO_CIS")
fits1 <- fits0; fits1$cis <- mkfit("OK",beta=.11,p=.02)
dec2 <- .cater_primary_decision(fits1,has_trans=TRUE,sibling_screen_performed=TRUE,
                                sibling_screen_complete=FALSE,n_active_siblings=0L)
stopifnot(dec2$model=="cis",abs(dec2$beta-.11)<1e-12,
          dec2$status=="SIBLING_SCREEN_INCOMPLETE_CIS_FALLBACK")
dec3 <- .cater_primary_decision(fits1,has_trans=TRUE,sibling_screen_performed=FALSE,
                                sibling_screen_complete=FALSE,n_active_siblings=0L)
stopifnot(dec3$model=="cis",dec3$status=="TRANS_UNSCREENED_CIS_FALLBACK")

# Deterministic local MVMR recovers known coefficients and conditional-F diagnostics converge.
B <- rbind(c(.20,.00),c(.10,.15),c(.00,.20),c(.12,.05))
colnames(B) <- c("X","Z")
theta <- c(X=.5,Z=.6)
by <- as.numeric(B %*% theta)
seB <- matrix(.02,nrow(B),ncol(B),dimnames=dimnames(B))
Cexp <- diag(2); dimnames(Cexp) <- list(c("X","Z"),c("X","Z"))
mf <- .cater_mvmr_fit(B,seB,by,rep(.02,nrow(B)),diag(nrow(B)),c("X","Z"),Cexp)
stopifnot(mf$status=="OK",abs(mf$beta["X"]-.5)<1e-10,abs(mf$beta["Z"]-.6)<1e-10)
stopifnot(isTRUE(mf$conditional_F_converged["X"]),mf$conditional_F_df["X"]==3L)

# scMORE adapter forwards the audited createRegulon public defaults exactly.
sa <- .cater_scmore_create_args()
stopifnot(identical(sa$n_targets,5),identical(sa$peak2gene_method,"Signac"),
          identical(sa$infer_method,"glm"),identical(sa$tss_upstream,100000),
          identical(sa$tss_downstream,0),identical(sa$exclude_exon_regions,TRUE))
stopifnot(!"conserved_regions" %in% names(sa))

# Omitting conserved_regions lets the upstream function evaluate its own default;
# specifying it forwards the object unchanged.
mock_create <- function(single_cell,n_targets=5,peak2gene_method="Signac",infer_method="glm",
                        tss_upstream=100000,tss_downstream=0,exclude_exon_regions=TRUE,
                        conserved_regions="UPSTREAM_DEFAULT") {
  list(grn=data.frame(TF="TF1",Target="X",Regions="chr1-1-2",Pval=.01),
       tf_names="TF1",seen=list(single_cell=single_cell,n_targets=n_targets,
       peak2gene_method=peak2gene_method,infer_method=infer_method,tss_upstream=tss_upstream,
       tss_downstream=tss_downstream,exclude_exon_regions=exclude_exon_regions,
       conserved_regions=conserved_regions))
}
mock1 <- .cater_scmore_call_create_regulon("CELL_SUBSET",sa,create_fun=mock_create)
stopifnot(identical(mock1$seen$single_cell,"CELL_SUBSET"),
          identical(mock1$seen$conserved_regions,"UPSTREAM_DEFAULT"))
custom_regions <- structure(list(id=1L),class="mock_regions")
sa2 <- .cater_scmore_create_args(conserved_regions=custom_regions)
mock2 <- .cater_scmore_call_create_regulon("CELL_SUBSET",sa2,create_fun=mock_create)
stopifnot(identical(mock2$seen$conserved_regions,custom_regions))
.cater_scmore_validate_output(mock2,"MockCell")
stopifnot(identical(.cater_scmore_safe_name("CD8+ T / effector"),"CD8_T_effector"))

cat("CATER-MR v0.5 core-estimator + scMORE adapter regression tests passed\n")