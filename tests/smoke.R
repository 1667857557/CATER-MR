source("CATER_MR.R")

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

# Deterministic local MVMR recovers known coefficients.
B <- rbind(c(.20,.00),c(.10,.15),c(.00,.20),c(.12,.05))
colnames(B) <- c("X","Z")
theta <- c(X=.5,Z=.6)
by <- as.numeric(B %*% theta)
seB <- matrix(.02,nrow(B),ncol(B),dimnames=dimnames(B))
Cexp <- diag(2); dimnames(Cexp) <- list(c("X","Z"),c("X","Z"))
mf <- .cater_mvmr_fit(B,seB,by,rep(.02,nrow(B)),diag(nrow(B)),c("X","Z"),Cexp)
stopifnot(mf$status=="OK",abs(mf$beta["X"]-.5)<1e-10,abs(mf$beta["Z"]-.6)<1e-10)

cat("CATER-MR v0.4 math smoke tests passed\n")
