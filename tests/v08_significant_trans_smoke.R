source("CATER_MR.R")

# Significant-only trans input must preserve the same SNP when it is reported for
# multiple genes; uniqueness is per SNP-gene pair, not global across the catalog.
trans <- data.frame(
  GENE=c("X","Z","X"), SNP=c("rsT","rsT","rsWeak"), CHR=c(2,2,2), BP=c(5000,5000,5100),
  A1="A", A2="G", EAF=.2, BETA=c(.20,.15,.05), SE=c(.02,.02,.02),
  P=c(1e-20,1e-12,1e-4), N=500, stringsAsFactors=FALSE
)
tx <- .cater_standardize_trans_hits(trans,qtl_n=500)
stopifnot(nrow(tx)==3L,sum(tx$snp=="rsT")==2L,setequal(tx$gene,c("X","Z")))

# Full cis and significant-only trans are different observation mechanisms but
# share the same instrument-eligibility threshold.
td <- tempfile("cater_v08_"); dir.create(td)
cis <- data.frame(
  SNP=c("rsC","rsCweak"), CHR=1, BP=c(1000,1050), A1=c("A","C"), A2=c("G","T"),
  EAF=c(.2,.3), BETA=c(.20,.04), SE=c(.02,.02), P=c(1e-20,1e-3), N=500
)
con <- gzfile(file.path(td,"X.txt.gz"),"wt")
utils::write.table(cis,con,sep="\t",quote=FALSE,row.names=FALSE); close(con)
ann <- data.frame(symbol=c("X","TF1"),chr=c("1","2"),tss=c(1000,5000),stringsAsFactors=FALSE)
cand <- .cater_target_candidates("X","TF1",ann,td,500,tx,5e-8,100,200,TRUE,FALSE)
stopifnot(setequal(cand$qtl$snp,c("rsC","rsT")))
stopifnot(cand$map$source[match("rsC",cand$map$snp)]=="cis")
stopifnot(cand$map$source[match("rsT",cand$map$snp)]=="trans")
stopifnot(cand$map$provenance[match("rsT",cand$map$snp)]=="REPORTED_TRANS")
stopifnot(!"rsWeak" %in% cand$qtl$snp,!"rsCweak" %in% cand$qtl$snp)

# A reported trans hit outside every direct-parent TF locus is biologically
# ineligible and must not enter CATER merely because it is genome-wide significant.
outside <- trans[1,,drop=FALSE]; outside$SNP <- "rsOutside"; outside$BP <- 9000; outside$P <- 1e-30
outside_std <- .cater_standardize_trans_hits(outside,qtl_n=500)
cand2 <- .cater_target_candidates("X","TF1",ann,td,500,outside_std,5e-8,100,200,TRUE,FALSE)
stopifnot(identical(cand2$qtl$snp,"rsC"),all(cand2$map$source=="cis"))

# Hit-only absence is not evidence of zero effect: sibling screening remains
# incomplete unless every selected trans IV has a usable sibling association.
h <- data.frame(snp=c("g1","g2"),source="trans",parent_tf="TF1",locus_id="TF:TF1",
                provenance="REPORTED_TRANS",bx=c(.2,.2),bx_se=.02,by=.1,by_se=.02)
grn <- data.frame(TF=c("TF1","TF1"),Target=c("X","Z"))
ref <- data.frame(snp=c("g1","g2"),ld_a1="A",ld_a2="G")
R <- diag(2); dimnames(R) <- list(c("g1","g2"),c("g1","g2"))
zrep <- data.frame(GENE="Z",SNP="g1",CHR=2,BP=5000,A1="A",A2="G",EAF=.2,
                   BETA=.2,SE=.02,P=1e-20,N=500,stringsAsFactors=FALSE)
zrep <- .cater_standardize_trans_hits(zrep,qtl_n=500)
sib <- .cater_sibling_screen("X",h,grn,ref,R,td,500,.05,TRUE,zrep,NULL)
stopifnot(sib$n_iv_missing==1L,!sib$complete,isTRUE(sib$table$co_perturbation_detected))

# The v0.8 strict contract explicitly permits full cis + significance-censored trans.
manifest <- list(grn_build="GRCh38",eqtl_build="GRCh38",ld_build="GRCh38",
                 eqtl_n_unit="donors",eqtl_ancestry="EUR",ld_ancestry="EUR",
                 cis_full_summary=TRUE,trans_data_mode="significant_only",
                 trans_reporting_threshold=5e-8)
stopifnot(is.list(.cater_validate_input_manifest(manifest,TRUE)))

# The custom correlated conditional-F remains available as a diagnostic but is
# explicitly labelled experimental; it is not silently treated as validated.
B <- rbind(c(.20,.00),c(.10,.15),c(.00,.20),c(.12,.05)); colnames(B)<-c("X","Z")
SE <- matrix(.02,nrow(B),ncol(B),dimnames=dimnames(B)); C <- diag(2); dimnames(C)<-list(c("X","Z"),c("X","Z"))
cf <- .cater_conditional_f(B,SE,diag(nrow(B)),c("X","Z"),C)
stopifnot(identical(attr(cf,"method"),"EXPERIMENTAL_CORRELATED_IV_CONDITIONAL_F"))

# Public API exposes the v0.8 data contract while preserving the historical
# positional slots. The old implementations remain removed; these formals exist
# only to prevent silent positional rebinding in pre-v0.8 calls.
fml <- names(formals(cater_mr))
stopifnot(all(c("trans_eqtl","instrument_p","trans_reporting_p","ld_clump_r2",
                "cross_effect_lookup","accept_experimental_conditional_f") %in% fml))
legacy_prefix <- c("grn","eqtl_dir","outcome","gene_annotation","ld_bfile",
                   "manc_cojo_bin","targets","cell_type","trait","cis_window","tf_window",
                   "cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads","qtl_n")
stopifnot(identical(fml[seq_along(legacy_prefix)],legacy_prefix))
legacy_tail <- c("outdir","drop_palindromic","verbose","primary_policy")
i <- match("outdir",fml); stopifnot(identical(fml[i:(i+3L)],legacy_tail))
plink_i <- match("plink_bin",fml)
stopifnot(identical(fml[(plink_i+1L):(plink_i+3L)],c("enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z")))
stopifnot("ld_threads" %in% fml)
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))

unlink(td,recursive=TRUE,force=TRUE)
cat("CATER-MR v0.8 significant-trans architecture smoke tests passed\n")


# Codex review regression tests (PR #13)
# Significant-only manifests must never silently fall back to legacy full-summary trans.
m_sig <- list(trans_data_mode="significant_only")
stopifnot(inherits(try(.cater_resolve_trans_contract(m_sig,NULL,5e-8,5e-8),silent=TRUE),"try-error"))
stopifnot(identical(.cater_resolve_trans_contract(m_sig,data.frame(),5e-8,5e-8),"significant_only"))

# A source reporting threshold stricter than the requested instrument threshold is not identifiable.
stopifnot(inherits(try(.cater_resolve_trans_contract(m_sig,data.frame(),5e-8,1e-10),silent=TRUE),"try-error"))

# Network primary requires explicit policy + conditional-F acceptance + both public gates.
mk_review_fit <- function(status="OK",b=.2,se=.05,p=.01) data.frame(n_iv=2L,beta=b,se=se,p=p,Q=NA,Q_p=NA,information=10,joint_wald_per_df=5,effective_F=5,precision_information=10,mean_F=20,min_F=15,status=status)
fits_review <- list(cis=mk_review_fit(),trans=mk_review_fit(),combined=mk_review_fit())
net_review <- list(status="OK",primary_eligible=TRUE,beta=c(X=.3),se=c(X=.06),p=c(X=.001))
d_blocked <- .cater_primary_decision(fits_review,net_review,"X",TRUE,TRUE,TRUE,1L,"screened_cater",allow_network_primary=FALSE,sibling_screen_independent=TRUE,sibling_screen_testable=TRUE)
d_allowed <- .cater_primary_decision(fits_review,net_review,"X",TRUE,TRUE,TRUE,1L,"screened_cater",allow_network_primary=TRUE,sibling_screen_independent=TRUE,sibling_screen_testable=TRUE)
stopifnot(d_blocked$model=="cis",d_allowed$model=="network")

# The MVMR builder creates a per-target directory before calling LD selection.
b_mvmr <- paste(deparse(body(.cater_build_mvmr)),collapse="\n")
stopifnot(grepl("dir.create\\(mvmr_dir",b_mvmr),grepl("file.path\\(mvmr_dir, \"ld\"\\)",b_mvmr))

# Candidate frames are normalized to a common summary-stat schema before rbind.
stopifnot(grepl("d <- d\\[, core, drop = FALSE\\]",gsub(";", "; ", b_mvmr),fixed=FALSE) || grepl("d<-d\\[,core,drop=FALSE\\]",b_mvmr))
