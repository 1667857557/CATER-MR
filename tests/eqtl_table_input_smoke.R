source("CATER_MR.R")
source("CATER_EQTL_INPUT.R")

combined <- data.frame(
  GENE = c("X","TF1","X","TF1"),
  SNP = c("rs1","rs1","rs2","rs2"),
  CHR = c(1,1,1,1),
  BP = c(100,100,200,200),
  A1 = c("A","A","C","C"),
  A2 = c("G","G","T","T"),
  EAF = c(.2,.2,.3,.3),
  BETA = c(.20,.10,.15,.05),
  SE = c(.02,.03,.03,.04),
  P = c(1e-20,1e-4,1e-8,.2),
  N = 500,
  stringsAsFactors = FALSE
)

# GENE is detected automatically and materialization preserves the exact
# per-gene full-summary statistics used by the existing core parser.
cache <- cater_prepare_eqtl_table(combined, verbose=FALSE)
stopifnot(cache$n_genes == 2L, cache$n_rows == 4L, cache$gene_col == "GENE")
stopifnot(all(file.exists(file.path(cache$eqtl_dir,c("X.txt.gz","TF1.txt.gz")))))
qx_cache <- .cater_read_gene("X",cache$eqtl_dir,NULL)
qx_direct <- .cater_standardize_sumstats(
  combined[combined$GENE=="X",setdiff(names(combined),"GENE"),drop=FALSE],
  label="X eQTL"
)
stopifnot(identical(qx_cache,qx_direct))

# The same contract works from a single compressed merged file.
f <- tempfile(fileext=".txt.gz")
con <- gzfile(f,"wt")
utils::write.table(combined,con,sep="\t",quote=FALSE,row.names=FALSE)
close(con)
cache2 <- cater_prepare_eqtl_table(f,verbose=FALSE)
qtf_cache <- .cater_read_gene("TF1",cache2$eqtl_dir,NULL)
qtf_direct <- .cater_standardize_sumstats(
  combined[combined$GENE=="TF1",setdiff(names(combined),"GENE"),drop=FALSE],
  label="TF1 eQTL"
)
stopifnot(identical(qtf_cache,qtf_direct))

# Common alternate gene-column names are accepted case-insensitively.
alt <- combined; names(alt)[names(alt)=="GENE"] <- "gene_symbol"
cache3 <- cater_prepare_eqtl_table(alt,verbose=FALSE)
stopifnot(cache3$gene_col == "gene_symbol",cache3$n_genes==2L)

# A custom gene column can be supplied explicitly.
custom <- combined; names(custom)[names(custom)=="GENE"] <- "EXPOSURE"
cache4 <- cater_prepare_eqtl_table(custom,gene_col="EXPOSURE",verbose=FALSE)
stopifnot(cache4$gene_col == "EXPOSURE")

# Missing gene identifiers and path-unsafe identifiers fail explicitly.
bad <- combined; bad$GENE[1] <- NA_character_
stopifnot(inherits(try(cater_prepare_eqtl_table(bad,verbose=FALSE),silent=TRUE),"try-error"))
bad2 <- combined; bad2$GENE[1] <- "A/B"
stopifnot(inherits(try(cater_prepare_eqtl_table(bad2,verbose=FALSE),silent=TRUE),"try-error"))

# Wrapper forwards the materialized directory to the unchanged core estimator.
core_cater_mr <- cater_mr
seen <- NULL
cater_mr <- function(grn,eqtl_dir,outcome,...) {
  seen <<- list(grn=grn,eqtl_dir=eqtl_dir,outcome=outcome,dots=list(...))
  stopifnot(dir.exists(eqtl_dir),file.exists(file.path(eqtl_dir,"X.txt.gz")))
  list(ok=TRUE)
}
res <- cater_mr_from_table(
  grn=data.frame(TF="TF1",Target="X"),
  eqtl_table=combined,
  outcome=data.frame(SNP="rs1",A1="A",A2="G",BETA=.1,SE=.02),
  materialize_verbose=FALSE,
  trait="Disease"
)
stopifnot(isTRUE(res$ok),seen$dots$trait=="Disease")
cater_mr <- core_cater_mr

unlink(cache$eqtl_dir,recursive=TRUE,force=TRUE)
unlink(cache2$eqtl_dir,recursive=TRUE,force=TRUE)
unlink(cache3$eqtl_dir,recursive=TRUE,force=TRUE)
unlink(cache4$eqtl_dir,recursive=TRUE,force=TRUE)
unlink(f,force=TRUE)

cat("Combined GENE-table eQTL input smoke tests passed\n")
