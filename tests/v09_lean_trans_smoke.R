source("CATER_MR.R")

stopifnot(identical(.CATER_VERSION,"0.9.0"))
fml <- formals(cater_mr)
stopifnot(identical(fml$enable_sibling_screen,FALSE))
stopifnot(identical(fml$enable_mvmr,FALSE))
stopifnot(all(c("tf_anchor_p","max_reported_trans_targets","trans_set") %in% names(fml)))

# Build minimal full-cis files for two candidate parent TFs.
td <- tempfile("cater_v09_"); dir.create(td)
write_gene <- function(gene,dat) {
  con <- gzfile(file.path(td,paste0(gene,".txt.gz")),"wt")
  utils::write.table(dat,con,sep="\t",quote=FALSE,row.names=FALSE); close(con)
}
base <- data.frame(SNP=c("rsCore","rsHot","rsAmb"),CHR=2,BP=c(5000,5100,5200),
                   A1="A",A2="G",EAF=.2,BETA=c(.20,.18,.16),SE=.02,
                   P=c(1e-20,1e-18,1e-16),N=500,stringsAsFactors=FALSE)
write_gene("TF1",base)
tf2 <- base[3,,drop=FALSE]; tf2$BETA <- .12; tf2$P <- 1e-12
write_gene("TF2",tf2)

qtl <- .cater_standardize_sumstats(base,n_default=500,label="target trans")
qtl$provenance <- "REPORTED_TRANS"
map <- data.frame(snp=qtl$snp,source="trans",
                  parent_tf=c("TF1","TF1","TF1;TF2"),
                  locus_id=c("TF:TF1","TF:TF1","TF:TF1+TF2"),
                  provenance="REPORTED_TRANS",stringsAsFactors=FALSE)
trans <- data.frame(GENE=c("X","X","Z","X"),
                    SNP=c("rsCore","rsHot","rsHot","rsAmb"),
                    CHR=2,BP=c(5000,5100,5100,5200),A1="A",A2="G",EAF=.2,
                    BETA=c(.20,.18,.11,.16),SE=.02,P=c(1e-20,1e-18,1e-12,1e-16),N=500,
                    stringsAsFactors=FALSE)
th <- .cater_standardize_trans_hits(trans,qtl_n=500)

core <- .cater_qualify_trans_candidates("X",qtl,map,td,500,th,5e-8,1L,"core")
stopifnot(identical(core$qtl$snp,"rsCore"))
stopifnot(core$counts$n_trans_raw==3L)
stopifnot(core$counts$n_trans_anchor_pass==3L)
stopifnot(core$counts$n_trans_unique_parent==2L)
stopifnot(core$counts$n_trans_core==1L)
stopifnot(core$counts$n_trans_extended==2L)
stopifnot(core$counts$n_parent_tf_core==1L)
stopifnot(core$table$reason[core$table$snp=="rsHot"]=="TRANS_REPORTED_HOTSPOT")
stopifnot(core$table$reason[core$table$snp=="rsAmb"]=="TRANS_AMBIGUOUS_PARENT")
stopifnot(core$map$parent_tf[core$map$snp=="rsCore"]=="TF1")

extended <- .cater_qualify_trans_candidates("X",qtl,map,td,500,th,5e-8,1L,"extended")
stopifnot(setequal(extended$qtl$snp,c("rsCore","rsHot")))
stopifnot(!"rsAmb" %in% extended$qtl$snp)

# A positional TF-locus match without an observed TF cis-eQTL is not eligible.
q_no <- qtl[1,,drop=FALSE]; q_no$snp <- "rsNoAnchor"
m_no <- data.frame(snp="rsNoAnchor",source="trans",parent_tf="TF1",locus_id="TF:TF1",
                   provenance="REPORTED_TRANS",stringsAsFactors=FALSE)
z_no <- .cater_qualify_trans_candidates("X",q_no,m_no,td,500,NULL,5e-8,1L,"core")
stopifnot(nrow(z_no$qtl)==0L,z_no$table$reason=="TRANS_NO_TF_CIS_ANCHOR")

# Lean primary decision is cis-only while the workflow is active; the historical
# internal helper remains available outside the public workflow for compatibility.
mkfit <- function(status="OK",b=.2,se=.05,p=.01,F=20) data.frame(n_iv=2L,beta=b,se=se,p=p,Q=NA,Q_p=NA,
  information=40,joint_wald_per_df=F,effective_F=F,precision_information=25,mean_F=F,min_F=F,status=status)
fits <- list(cis=mkfit(),trans=mkfit(b=.4),combined=mkfit(b=.3))
old <- .CATER_LEAN_STATE$active; .CATER_LEAN_STATE$active <- TRUE
pd <- .cater_primary_decision(fits,NULL,"X",TRUE,TRUE,TRUE,0L,"screened_cater")
.CATER_LEAN_STATE$active <- old
stopifnot(pd$model=="cis",pd$beta==fits$cis$beta)

# Post-processing exposes model-specific F statistics and maps the legacy generic
# effective_F to the actual primary (cis), not to the augmented fit.
out <- tempfile("cater_v09_out_"); dir.create(out); dir.create(file.path(out,"mechanism"))
rr <- rbind(
  data.frame(cell_type="ct",target="X",trait="Y",model="cis",mkfit(F=12),q=NA),
  data.frame(cell_type="ct",target="X",trait="Y",model="trans",mkfit(F=30),q=NA),
  data.frame(cell_type="ct",target="X",trait="Y",model="combined",mkfit(se=.04,F=22),q=NA)
)
tg <- data.frame(target="X",primary_model="cis",primary_p=.01,effective_F=22,stringsAsFactors=FALSE)
.CATER_LEAN_STATE$diagnostics <- list()
pp <- .cater_lean_postprocess(list(results=rr,targets=tg),out,5e-8,1L,"core")
stopifnot(any(pp$results$model=="augmented"))
stopifnot(pp$targets$cis_effective_F==12)
stopifnot(pp$targets$augmented_effective_F==22)
stopifnot(pp$targets$primary_effective_F==12)
stopifnot(pp$targets$effective_F==12)
stopifnot(pp$targets$primary_model=="cis")
stopifnot(isTRUE(all.equal(pp$targets$augmented_se_reduction,0.2)))

unlink(td,recursive=TRUE,force=TRUE); unlink(out,recursive=TRUE,force=TRUE)
cat("CATER-MR Lean v0.9 trans qualification smoke tests passed\n")
