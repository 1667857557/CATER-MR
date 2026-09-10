# CATER-MR v0.6 evidence-safe hardening layer.
#
# This module intentionally overrides only functions whose behavior changed after
# the independent + literature + systems audit. The v0.5 core remains frozen in
# R/CATER_MR_core_v05.R so each hardening change is reviewable and testable.

.CATER_VERSION <- "0.6.0"

.cater_atomic_write_table <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(.Machine$integer.max, 1L))
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)
  utils::write.table(x, tmp, ...)
  if (!file.rename(tmp, path)) .cater_stop("Could not atomically replace output: %s", path)
  invisible(path)
}

.cater_standardize_sumstats <- function(x, n_default = NULL, label = "summary data",
                                        require_position = TRUE) {
  if (!is.data.frame(x)) .cater_stop("%s must be a data.frame", label)
  n_input <- nrow(x)
  snp <- .cater_pick_col(x, c("SNP","rsid","variant_id","variant","ID"), label = label)
  chr <- .cater_pick_col(x, c("CHR","chrom","chromosome"), required = require_position, label = label)
  pos <- .cater_pick_col(x, c("BP","POS","position","base_pair_location"), required = require_position, label = label)
  a1 <- .cater_pick_col(x, c("A1","EA","effect_allele","ALT"), label = label)
  a2 <- .cater_pick_col(x, c("A2","NEA","other_allele","non_effect_allele","REF"), label = label)
  b <- .cater_pick_col(x, c("b","beta","BETA","effect","estimate"), label = label)
  se <- .cater_pick_col(x, c("se","SE","stderr","standard_error"), label = label)
  p <- .cater_pick_col(x, c("p","P","pval","p_value","pvalue"), required = FALSE)
  eaf <- .cater_pick_col(x, c("freq","EAF","eaf","effect_allele_frequency","AF"), required = FALSE)
  n <- .cater_pick_col(x, c("N","n","samplesize","sample_size"), required = FALSE)

  out <- data.frame(
    snp = as.character(x[[snp]]),
    chr = if (is.null(chr)) NA_character_ else sub("^chr", "", as.character(x[[chr]]), ignore.case = TRUE),
    pos = if (is.null(pos)) NA_real_ else suppressWarnings(as.numeric(x[[pos]])),
    a1 = toupper(as.character(x[[a1]])),
    a2 = toupper(as.character(x[[a2]])),
    beta = suppressWarnings(as.numeric(x[[b]])),
    se = suppressWarnings(as.numeric(x[[se]])),
    stringsAsFactors = FALSE
  )
  out$p <- if (is.null(p)) 2 * stats::pnorm(-abs(out$beta / out$se)) else suppressWarnings(as.numeric(x[[p]]))
  out$eaf <- if (is.null(eaf)) NA_real_ else suppressWarnings(as.numeric(x[[eaf]]))
  out$n <- if (is.null(n)) {
    if (is.null(n_default)) NA_real_ else rep(as.numeric(n_default), nrow(x))
  } else suppressWarnings(as.numeric(x[[n]]))

  valid_snp <- !is.na(out$snp) & nzchar(out$snp)
  valid_pos <- !require_position | is.finite(out$pos)
  valid_beta_se <- is.finite(out$beta) & is.finite(out$se) & out$se > 0
  valid_p <- is.finite(out$p) & out$p >= 0 & out$p <= 1
  valid_allele <- out$a1 %in% c("A","C","G","T") & out$a2 %in% c("A","C","G","T") & out$a1 != out$a2
  keep <- valid_snp & valid_pos & valid_beta_se & valid_p & valid_allele
  out <- out[keep,,drop = FALSE]
  dup <- duplicated(out$snp) | duplicated(out$snp, fromLast = TRUE)
  n_dup <- sum(dup)
  out <- out[!dup,,drop = FALSE]
  rownames(out) <- NULL
  attr(out, "qc") <- list(
    n_input = n_input,
    n_valid_pre_duplicate = sum(keep),
    n_duplicate_rows_removed = n_dup,
    n_output = nrow(out),
    n_invalid_snp = sum(!valid_snp),
    n_invalid_position = sum(!valid_pos),
    n_invalid_beta_se = sum(!valid_beta_se),
    n_invalid_p = sum(!valid_p),
    n_invalid_allele = sum(!valid_allele)
  )
  out
}

.cater_validate_unique_annotation <- function(out, label = "gene annotation") {
  out <- out[!is.na(out$symbol) & nzchar(out$symbol) & !is.na(out$chr) &
               nzchar(out$chr) & is.finite(out$tss),,drop = FALSE]
  if (!nrow(out)) return(out)
  key <- paste(out$chr, format(out$tss, scientific = FALSE, trim = TRUE), sep = ":")
  amb <- vapply(split(key, out$symbol), function(z) length(unique(z)) > 1L, logical(1))
  if (any(amb)) {
    bad <- names(amb)[amb]
    .cater_stop("%s has multiple distinct chr/TSS values for %d symbols (e.g. %s). Supply one explicit gene-level TSS per symbol.",
                label, length(bad), paste(utils::head(bad, 10L), collapse = ", "))
  }
  out[!duplicated(out$symbol),,drop = FALSE]
}

.cater_annotation_from_grn <- function(grn_raw) {
  if (is.list(grn_raw) && !is.data.frame(grn_raw) && !is.null(grn_raw$grn)) grn_raw <- grn_raw$grn
  if (!is.data.frame(grn_raw)) return(NULL)
  tf <- .cater_pick_col(grn_raw,c("TF","tf","regulator"),required=FALSE)
  tg <- .cater_pick_col(grn_raw,c("Target","target","gene"),required=FALSE)
  tf_chr <- .cater_pick_col(grn_raw,c("TF_chr","tf_chr"),required=FALSE)
  tf_tss <- .cater_pick_col(grn_raw,c("TF_tss","tf_tss"),required=FALSE)
  tg_chr <- .cater_pick_col(grn_raw,c("Target_chr","target_chr"),required=FALSE)
  tg_tss <- .cater_pick_col(grn_raw,c("Target_tss","target_tss"),required=FALSE)
  if (any(vapply(list(tf,tg,tf_chr,tf_tss,tg_chr,tg_tss), is.null, logical(1)))) return(NULL)
  a <- rbind(
    data.frame(symbol=as.character(grn_raw[[tf]]), chr=as.character(grn_raw[[tf_chr]]), tss=suppressWarnings(as.numeric(grn_raw[[tf_tss]]))),
    data.frame(symbol=as.character(grn_raw[[tg]]), chr=as.character(grn_raw[[tg_chr]]), tss=suppressWarnings(as.numeric(grn_raw[[tg_tss]])))
  )
  a$chr <- sub("^chr", "", a$chr, ignore.case = TRUE)
  .cater_validate_unique_annotation(a, "GRN coordinate columns")
}

.cater_standardize_annotation <- function(annotation) {
  if (is.null(annotation)) return(NULL)
  if (!is.data.frame(annotation)) .cater_stop("gene_annotation must be a data.frame")
  sym <- .cater_pick_col(annotation,c("symbol","gene","gene_symbol","SYMBOL"),label="gene annotation symbol")
  chr <- .cater_pick_col(annotation,c("chr","CHR","chrom","chromosome"),label="gene annotation chromosome")
  # Deliberately reject txStart as a TSS alias. txStart is wrong on the negative
  # strand unless a transcript-aware conversion has already been performed.
  tss <- .cater_pick_col(annotation,c("tss","TSS"),label="gene annotation TSS")
  out <- data.frame(
    symbol = as.character(annotation[[sym]]),
    chr = sub("^chr", "", as.character(annotation[[chr]]), ignore.case = TRUE),
    tss = suppressWarnings(as.numeric(annotation[[tss]])),
    stringsAsFactors = FALSE
  )
  .cater_validate_unique_annotation(out, "gene_annotation")
}

.cater_validate_ld <- function(ld, label = "LD", tol = 1e-7) {
  ld <- as.matrix(ld)
  if (!is.numeric(ld) || nrow(ld) != ncol(ld)) .cater_stop("%s must be a square numeric matrix", label)
  if (any(!is.finite(ld))) .cater_stop("%s contains non-finite values", label)
  if (max(abs(ld - t(ld))) > tol) .cater_stop("%s is not symmetric", label)
  if (any(abs(diag(ld) - 1) > tol)) .cater_stop("%s diagonal is not one", label)
  if (any(abs(ld) > 1 + tol)) .cater_stop("%s contains correlations outside [-1,1]", label)
  ev <- eigen((ld + t(ld))/2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -tol * max(1, max(abs(ev)))) .cater_stop("%s is not positive semidefinite", label)
  ld
}

.cater_read_manc_ldr <- function(path, selected) {
  selected <- unique(as.character(selected))
  if (!length(selected)) return(matrix(numeric(), 0, 0))
  if (!file.exists(path)) .cater_stop("Manc-COJO LD output not found: %s", path)
  lines <- trimws(readLines(path, warn = FALSE)); lines <- lines[nzchar(lines)]
  blocks <- list(); i <- 1L
  while (i <= length(lines)) {
    if (startsWith(lines[i], "#")) { i <- i + 1L; next }
    if (!grepl("^SNP(\\s|$)", lines[i])) { i <- i + 1L; next }
    header <- strsplit(lines[i], "\\s+")[[1L]]; cols <- header[-1L]; i <- i + 1L
    rows <- list()
    while (i <= length(lines) && !startsWith(lines[i], "#") && !grepl("^SNP(\\s|$)", lines[i])) {
      tok <- strsplit(lines[i], "\\s+")[[1L]]
      if (length(tok) >= 2L) rows[[length(rows)+1L]] <- tok
      i <- i + 1L
    }
    if (!length(rows) || !length(cols)) next
    rn <- vapply(rows, `[`, character(1), 1L)
    if (anyDuplicated(rn) || anyDuplicated(cols)) .cater_stop("Duplicate SNP in Manc-COJO LD block: %s", path)
    vals <- do.call(rbind, lapply(rows, function(z) suppressWarnings(as.numeric(z[-1L]))))
    if (ncol(vals) != length(cols) || length(rn) != nrow(vals) || !setequal(rn, cols))
      .cater_stop("Malformed/non-square Manc-COJO LD block in %s", path)
    rownames(vals) <- rn; colnames(vals) <- cols; vals <- vals[rn,rn,drop=FALSE]
    blocks[[length(blocks)+1L]] <- .cater_validate_ld(vals, paste0("Manc-COJO LD block ", length(blocks)+1L))
  }
  if (!length(blocks)) .cater_stop("No LD blocks found in Manc-COJO output: %s", path)
  observed <- unlist(lapply(blocks, rownames), use.names = FALSE)
  if (anyDuplicated(observed)) .cater_stop("A SNP appears in multiple Manc-COJO LD blocks: %s", path)
  missing <- setdiff(selected, observed)
  if (length(missing)) .cater_stop("Manc-COJO LD output is missing %d requested SNP(s): %s",
                                   length(missing), paste(utils::head(missing,10L), collapse=", "))
  ld <- matrix(0, length(selected), length(selected), dimnames = list(selected, selected))
  for (b in blocks) {
    common <- intersect(rownames(b), selected)
    if (length(common)) ld[common,common] <- b[common,common,drop=FALSE]
  }
  .cater_validate_ld(ld, "assembled Manc-COJO LD")
}

.cater_read_jma_ref <- function(path) {
  if (!file.exists(path)) return(data.frame())
  x <- utils::read.table(path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  if (!nrow(x)) return(data.frame())
  snp <- .cater_pick_col(x,"SNP",label="Manc-COJO SNP")
  a1 <- .cater_pick_col(x,"A1",label="Manc-COJO A1")
  a2 <- .cater_pick_col(x,"A2",label="Manc-COJO A2")
  out <- data.frame(snp=as.character(x[[snp]]), ld_a1=toupper(as.character(x[[a1]])),
                    ld_a2=toupper(as.character(x[[a2]])), stringsAsFactors=FALSE)
  out <- out[!is.na(out$snp) & nzchar(out$snp),,drop=FALSE]
  if (anyDuplicated(out$snp)) .cater_stop("Duplicate SNP IDs in Manc-COJO .jma.cojo: %s", path)
  out
}

.cater_remove_prefix_outputs <- function(prefix) {
  z <- Sys.glob(paste0(prefix, "*")); if (length(z)) unlink(z, recursive=TRUE, force=TRUE)
  invisible(NULL)
}

.cater_quote_args <- function(args) {
  vapply(as.character(args), function(z) if (grepl("[[:space:]]",z)) shQuote(z) else z,
         character(1), USE.NAMES=FALSE)
}

.cater_run_cojo <- function(qtl,candidate_map,ld_bfile,manc_cojo_bin,cojo_p,
                            cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose) {
  if (length(ld_bfile)!=1L || !nzchar(ld_bfile)) .cater_stop("Current CATER-MR expects one COJO LD cohort")
  if (!all(file.exists(paste0(ld_bfile,c(".bed",".bim",".fam"))))) .cater_stop("ld_bfile must be a PLINK prefix: %s",ld_bfile)
  exe <- Sys.which(manc_cojo_bin); if (!nzchar(exe)) .cater_stop("Cannot find Manc-COJO executable '%s'",manc_cojo_bin)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE); .cater_remove_prefix_outputs(prefix)
  ma <- paste0(prefix,".sumstat"); extract <- paste0(prefix,".candidate.snplist")
  .cater_write_cojo_ma(qtl,ma); writeLines(unique(candidate_map$snp),extract)
  sp <- paste0(prefix,".select")
  args <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",extract,"--cojo-slct",
            "--cojo-p",format(cojo_p,scientific=TRUE),"--cojo-wind",as.character(as.integer(cojo_wind_kb)),
            "--cojo-collinear",as.character(cojo_collinear),"--thread-num",as.character(as.integer(cojo_threads)),"--out",sp)
  st <- system2(exe,args=.cater_quote_args(args),stdout=paste0(sp,".stdout"),stderr=paste0(sp,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO selection failed for %s",prefix)
  selected_step <- .cater_read_jma_ref(paste0(sp,".jma.cojo"))
  if (!nrow(selected_step)) return(list(selected=selected_step,selected_step=selected_step,joint_removed=character(),ld=matrix(numeric(),0,0)))
  if (nrow(selected_step)==1L) return(list(selected=selected_step,selected_step=selected_step,joint_removed=character(),ld=matrix(1,1,1,dimnames=list(selected_step$snp,selected_step$snp))))
  jp <- paste0(prefix,".joint"); .cater_remove_prefix_outputs(jp)
  args2 <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",paste0(sp,".jma.cojo"),"2","header",
             "--cojo-joint","--thread-num",as.character(as.integer(cojo_threads)),"--output-all","--out",jp)
  st <- system2(exe,args=.cater_quote_args(args2),stdout=paste0(jp,".stdout"),stderr=paste0(jp,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO joint/LD failed for %s",prefix)
  joint_ref <- .cater_read_jma_ref(paste0(jp,".jma.cojo"))
  if (!nrow(joint_ref)) .cater_stop("Manc-COJO joint retained no SNPs for %s",prefix)
  if (length(setdiff(joint_ref$snp,selected_step$snp))) .cater_stop("Manc-COJO joint output contains unselected SNPs for %s",prefix)
  ld <- .cater_read_manc_ldr(paste0(jp,".ldr.cojo"),joint_ref$snp)
  list(selected=joint_ref,selected_step=selected_step,joint_removed=setdiff(selected_step$snp,joint_ref$snp),ld=ld)
}

.cater_joint_ld <- function(qtl,snps,ld_bfile,manc_cojo_bin,cojo_threads,prefix,verbose) {
  requested <- unique(as.character(snps)); requested <- requested[requested %in% qtl$snp]
  if (!length(requested)) return(list(selected=data.frame(),requested=character(),removed=character(),ld=matrix(numeric(),0,0)))
  if (length(requested)==1L) {
    q <- qtl[match(requested,qtl$snp),,drop=FALSE]
    ref <- data.frame(snp=q$snp,ld_a1=q$a1,ld_a2=q$a2,stringsAsFactors=FALSE)
    return(list(selected=ref,requested=requested,removed=character(),ld=matrix(1,1,1,dimnames=list(requested,requested))))
  }
  exe <- Sys.which(manc_cojo_bin); if (!nzchar(exe)) .cater_stop("Cannot find Manc-COJO executable '%s'",manc_cojo_bin)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE); .cater_remove_prefix_outputs(prefix)
  ma <- paste0(prefix,".sumstat"); ex <- paste0(prefix,".snplist")
  .cater_write_cojo_ma(qtl,ma); writeLines(requested,ex)
  args <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",ex,"--cojo-joint",
            "--thread-num",as.character(as.integer(cojo_threads)),"--output-all","--out",prefix)
  st <- system2(exe,args=.cater_quote_args(args),stdout=paste0(prefix,".stdout"),stderr=paste0(prefix,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO joint LD failed for %s",prefix)
  ref <- .cater_read_jma_ref(paste0(prefix,".jma.cojo"))
  if (!nrow(ref)) return(list(selected=ref,requested=requested,removed=requested,ld=matrix(numeric(),0,0)))
  if (length(setdiff(ref$snp,requested))) .cater_stop("Manc-COJO joint output contains unrequested SNPs")
  ld <- .cater_read_manc_ldr(paste0(prefix,".ldr.cojo"),ref$snp)
  list(selected=ref,requested=requested,removed=setdiff(requested,ref$snp),ld=ld)
}

.cater_reorient_ld <- function(ld, from_ref, to_ref) {
  snps <- rownames(ld)
  fr <- from_ref[match(snps,from_ref$snp),,drop=FALSE]
  to <- to_ref[match(snps,to_ref$snp),,drop=FALSE]
  if (any(is.na(fr$snp)) || any(is.na(to$snp))) .cater_stop("Cannot reorient LD: allele reference missing SNPs")
  same <- fr$ld_a1==to$ld_a1 & fr$ld_a2==to$ld_a2
  swap <- fr$ld_a1==to$ld_a2 & fr$ld_a2==to$ld_a1
  if (any(!(same|swap))) .cater_stop("Cannot reorient LD because allele pairs differ")
  s <- ifelse(swap,-1,1); names(s) <- snps
  out <- ld * (s %o% s); dimnames(out) <- dimnames(ld)
  .cater_validate_ld(out,"reoriented LD")
}

.cater_strength <- function(dat,ld) {
  if(!nrow(dat)) return(list(information=NA_real_,joint_wald_per_df=NA_real_,effective_F=NA_real_,mean_F=NA_real_,min_F=NA_real_))
  f <- (dat$bx/dat$bx_se)^2
  Dx <- diag(dat$bx_se,nrow=nrow(dat)); Sx <- Dx%*%ld%*%Dx; inv <- .cater_inv(Sx)
  info <- if(is.null(inv)) NA_real_ else as.numeric(crossprod(dat$bx,inv%*%dat$bx))
  jw <- if(is.finite(info)) info/nrow(dat) else NA_real_
  list(information=info,joint_wald_per_df=jw,effective_F=jw,mean_F=mean(f),min_F=min(f))
}

.cater_givw <- function(dat,ld_y,ld_x=ld_y) {
  st <- .cater_strength(dat,ld_x); n <- nrow(dat)
  if(!n) return(data.frame(n_iv=0L,beta=NA,se=NA,p=NA,Q=NA,Q_p=NA,information=NA,joint_wald_per_df=NA,effective_F=NA,precision_information=NA,mean_F=NA,min_F=NA,status="NO_IV"))
  if(n==1L) {
    b <- dat$by/dat$bx; s <- abs(dat$by_se/dat$bx); prec <- (dat$bx/dat$by_se)^2
    return(data.frame(n_iv=1L,beta=b,se=s,p=2*stats::pnorm(-abs(b/s)),Q=NA,Q_p=NA,information=st$information,joint_wald_per_df=st$joint_wald_per_df,effective_F=st$effective_F,precision_information=prec,mean_F=st$mean_F,min_F=st$min_F,status="OK"))
  }
  Dy <- diag(dat$by_se,nrow=n); Oy <- Dy%*%ld_y%*%Dy; inv <- .cater_inv(Oy)
  if(is.null(inv)) return(data.frame(n_iv=n,beta=NA,se=NA,p=NA,Q=NA,Q_p=NA,information=st$information,joint_wald_per_df=st$joint_wald_per_df,effective_F=st$effective_F,precision_information=NA,mean_F=st$mean_F,min_F=st$min_F,status="LD_SINGULAR"))
  den <- as.numeric(crossprod(dat$bx,inv%*%dat$bx))
  if(!is.finite(den)||den<=0) return(data.frame(n_iv=n,beta=NA,se=NA,p=NA,Q=NA,Q_p=NA,information=st$information,joint_wald_per_df=st$joint_wald_per_df,effective_F=st$effective_F,precision_information=NA,mean_F=st$mean_F,min_F=st$min_F,status="INVALID_INFORMATION"))
  b <- as.numeric(crossprod(dat$bx,inv%*%dat$by))/den; s <- sqrt(1/den)
  r <- dat$by-b*dat$bx; q <- as.numeric(crossprod(r,inv%*%r))
  data.frame(n_iv=n,beta=b,se=s,p=2*stats::pnorm(-abs(b/s)),Q=q,Q_p=stats::pchisq(q,n-1,lower.tail=FALSE),information=st$information,joint_wald_per_df=st$joint_wald_per_df,effective_F=st$effective_F,precision_information=den,mean_F=st$mean_F,min_F=st$min_F,status="OK")
}

.cater_increment_fraction <- function(base,full,field) {
  a <- full[[field]]; b <- base[[field]]
  if(!is.finite(a)||a<=0) return(NA_real_)
  if(isTRUE(base$n_iv==0L)) b <- 0
  if(!is.finite(b)) return(NA_real_)
  d <- a-b
  if(d < -sqrt(.Machine$double.eps)*max(1,abs(a),abs(b))) return(NA_real_)
  max(0,min(1,d/a))
}

.cater_trans_exposure_signal_fraction <- function(cis_fit,combined_fit) .cater_increment_fraction(cis_fit,combined_fit,"information")
.cater_trans_information_fraction <- .cater_trans_exposure_signal_fraction
.cater_trans_mr_precision_fraction <- function(cis_fit,combined_fit) .cater_increment_fraction(cis_fit,combined_fit,"precision_information")

.cater_parent_tokens <- function(x) {
  lapply(strsplit(as.character(x),";",fixed=TRUE),function(v) unique(v[nzchar(v)]))
}

.cater_locus_diagnostics <- function(h,ld_y,ld_x=ld_y) {
  tr <- h[h$source=="trans",,drop=FALSE]
  if(!nrow(tr)) return(list(table=data.frame(),max_exposure_signal_weight=NA_real_,max_precision_weight=NA_real_,leave_one_max_delta=NA_real_,max_weight=NA_real_))
  full <- .cater_givw(h,ld_y,ld_x)
  tokens <- .cater_parent_tokens(tr$parent_tf); tfs <- sort(unique(unlist(tokens,use.names=FALSE)))
  tabs <- list(); deltas <- numeric()
  for(tf in tfs) {
    remove_tr <- vapply(tokens,function(v) tf%in%v,logical(1)); remove_snps <- tr$snp[remove_tr]
    keep <- !h$snp%in%remove_snps
    minus <- if(any(keep)) .cater_givw(h[keep,,drop=FALSE],.cater_subset_ld(ld_y,h$snp[keep]),.cater_subset_ld(ld_x,h$snp[keep])) else NULL
    delta <- if(!is.null(minus)&&full$status=="OK"&&minus$status=="OK") abs(full$beta-minus$beta) else NA_real_
    deltas <- c(deltas,delta)
    exp_inc <- if(!is.null(minus)&&is.finite(full$information)&&is.finite(minus$information)) full$information-minus$information else NA_real_
    pr_inc <- if(!is.null(minus)&&is.finite(full$precision_information)&&is.finite(minus$precision_information)) full$precision_information-minus$precision_information else NA_real_
    ew <- if(is.finite(exp_inc)&&is.finite(full$information)&&full$information>0) exp_inc/full$information else NA_real_
    pw <- if(is.finite(pr_inc)&&is.finite(full$precision_information)&&full$precision_information>0) pr_inc/full$precision_information else NA_real_
    tabs[[tf]] <- data.frame(tf=tf,n_removed_iv=length(remove_snps),removed_snps=paste(remove_snps,collapse=";"),exposure_signal_increment=exp_inc,exposure_signal_weight=ew,mr_precision_increment=pr_inc,mr_precision_weight=pw,leave_one_delta=delta,stringsAsFactors=FALSE)
  }
  z <- do.call(rbind,tabs); rownames(z) <- NULL
  maxew <- if(any(is.finite(z$exposure_signal_weight))) max(z$exposure_signal_weight,na.rm=TRUE) else NA_real_
  maxpw <- if(any(is.finite(z$mr_precision_weight))) max(z$mr_precision_weight,na.rm=TRUE) else NA_real_
  list(table=z,max_exposure_signal_weight=maxew,max_precision_weight=maxpw,leave_one_max_delta=if(any(is.finite(deltas)))max(deltas,na.rm=TRUE)else NA_real_,max_weight=maxew)
}

.cater_descendants <- function(grn,target) {
  seen <- character(); frontier <- target
  repeat {
    nxt <- unique(grn$Target[grn$TF %in% frontier]); nxt <- setdiff(nxt,c(target,seen))
    if(!length(nxt)) break
    seen <- unique(c(seen,nxt)); frontier <- nxt
  }
  seen
}

.cater_sibling_screen <- function(target,h,grn,cojo_ref,ld,eqtl_dir,qtl_n,sibling_fdr,preserve_total_effect=TRUE) {
  empty <- function(status="NO_TRANS_IV") list(table=data.frame(),active=character(),n_candidate=0L,n_incomplete=0L,n_qtl_missing=0L,n_iv_missing=0L,complete=FALSE,testable=FALSE,status=status)
  tr <- h[h$source=="trans",,drop=FALSE]; if(!nrow(tr)) return(empty("NO_TRANS_IV"))
  excluded <- if(isTRUE(preserve_total_effect)) .cater_descendants(grn,target) else character()
  tf_to_snps <- list()
  for(i in seq_len(nrow(tr))) for(tf in .cater_parent_tokens(tr$parent_tf[i])[[1L]]) tf_to_snps[[tf]] <- unique(c(tf_to_snps[[tf]],tr$snp[i]))
  sib_parents <- list()
  for(tf in names(tf_to_snps)) {
    sibs <- setdiff(unique(grn$Target[grn$TF==tf]),c(target,excluded))
    for(z in sibs) sib_parents[[z]] <- unique(c(sib_parents[[z]],tf))
  }
  if(!length(sib_parents)) return(empty("NO_TESTABLE_SIBLINGS"))
  out <- list(); qtl_missing <- 0L
  for(z in names(sib_parents)) {
    tfs <- sib_parents[[z]]; snps <- unique(unlist(tf_to_snps[tfs],use.names=FALSE)); requested <- length(snps)
    d <- .cater_exact_effects(z,snps,eqtl_dir,qtl_n,cojo_ref)
    if(is.null(d)) {
      qtl_missing <- qtl_missing+1L
      out[[z]] <- data.frame(sibling=z,parent_tf=paste(sort(tfs),collapse=";"),n_iv_requested=requested,n_iv_tested=0L,n_iv_missing=requested,Q=NA,df=NA,p=NA,status="QTL_NOT_AVAILABLE",complete=FALSE,testable=FALSE)
      next
    }
    tested <- nrow(d); missing_iv <- requested-tested
    if(!tested) {
      out[[z]] <- data.frame(sibling=z,parent_tf=paste(sort(tfs),collapse=";"),n_iv_requested=requested,n_iv_tested=0L,n_iv_missing=requested,Q=NA,df=NA,p=NA,status="SNP_NOT_AVAILABLE",complete=FALSE,testable=FALSE)
      next
    }
    R <- .cater_subset_ld(ld,d$snp); om <- .cater_omnibus(d$beta,d$se,R); test_ok <- is.finite(om["p"])
    status <- if(missing_iv>0L) { if(test_ok) "PARTIAL_SNP_COVERAGE" else "PARTIAL_TEST_FAILED" } else if(test_ok) "OK" else "TEST_FAILED"
    out[[z]] <- data.frame(sibling=z,parent_tf=paste(sort(tfs),collapse=";"),n_iv_requested=requested,n_iv_tested=tested,n_iv_missing=missing_iv,Q=om["Q"],df=om["df"],p=om["p"],status=status,complete=(missing_iv==0L&&test_ok),testable=test_ok)
  }
  tab <- do.call(rbind,out); rownames(tab) <- NULL; tab$q <- NA_real_
  ii <- which(is.finite(tab$p)); if(length(ii)) tab$q[ii] <- p.adjust(tab$p[ii],method="BH")
  tab$active <- is.finite(tab$q)&tab$q<sibling_fdr
  list(table=tab,active=tab$sibling[tab$active],n_candidate=nrow(tab),n_incomplete=sum(!tab$complete),n_qtl_missing=qtl_missing,n_iv_missing=sum(tab$n_iv_missing),complete=all(tab$complete),testable=any(tab$testable),status=if(all(tab$complete))"COMPLETE" else "INCOMPLETE")
}

.cater_conditional_f_v05 <- .cater_conditional_f
.cater_conditional_f <- function(B,SE,ld,exposure_names,exposure_corr=NULL,max_iter=200L,tol=1e-9) {
  ans <- .cater_conditional_f_v05(B,SE,ld,exposure_names,exposure_corr,max_iter,tol)
  attr(ans,"method") <- "EXPERIMENTAL_CORRELATED_IV_EXTENSION"
  ans
}

.cater_mvmr_fit <- function(B,seB,by,seY,ld,exposure_names,exposure_corr=NULL) {
  m <- nrow(B); p <- ncol(B)
  if(m<=p) return(list(status="MVMR_UNDERIDENTIFIED"))
  Dy <- diag(seY,nrow=m); Oy <- Dy%*%ld%*%Dy; invY <- .cater_inv(Oy)
  if(is.null(invY)) return(list(status="MVMR_LD_SINGULAR"))
  L <- tryCatch(t(chol(Oy)),error=function(e)NULL); if(is.null(L)) return(list(status="MVMR_LD_SINGULAR"))
  Bw <- forwardsolve(L,B); norms <- sqrt(colSums(Bw^2))
  if(any(!is.finite(norms)|norms<=0)) return(list(status="MVMR_RANK_DEFICIENT"))
  Bwd <- sweep(Bw,2,norms,"/"); rankB <- qr(Bwd)$rank; condB <- if(rankB<p) Inf else kappa(Bwd)
  if(rankB<p) return(list(status="MVMR_RANK_DEFICIENT",rank=rankB,condition=condB))
  info <- t(B)%*%invY%*%B; iv <- .cater_inv(info)
  if(is.null(iv)) return(list(status="MVMR_RANK_DEFICIENT",rank=rankB,condition=condB))
  th <- as.numeric(iv%*%t(B)%*%invY%*%by); names(th) <- exposure_names
  ses <- sqrt(diag(iv)); names(ses) <- exposure_names; ps <- 2*stats::pnorm(-abs(th/ses)); names(ps) <- exposure_names
  r <- by-as.numeric(B%*%th); Q <- as.numeric(t(r)%*%invY%*%r)
  cf <- .cater_conditional_f(B,seB,ld,exposure_names,exposure_corr)
  list(status="OK",beta=th,se=ses,p=ps,Q=Q,Q_p=stats::pchisq(Q,m-p,lower.tail=FALSE),rank=rankB,condition=condB,conditional_F=cf,conditional_F_Q=attr(cf,"Q"),conditional_F_df=attr(cf,"df"),conditional_F_converged=attr(cf,"converged"),conditional_strength_method=attr(cf,"method"),covariance_assumption=if(is.null(exposure_corr))"zero_within_SNP_exposure_covariance" else "user_exposure_corr")
}

.cater_validate_input_manifest <- function(input_manifest, strict=TRUE) {
  keys <- c("grn_build","eqtl_build","ld_build","outcome_build","eqtl_n_unit","eqtl_full_summary","eqtl_ancestry","outcome_ancestry","ld_ancestry","trans_qtl_qc","outcome_ld_ancestry","qtl_outcome_overlap")
  if(is.null(input_manifest)) {
    if(isTRUE(strict)) .cater_stop("strict_input_contract=TRUE requires input_manifest; see INPUT_CONTRACT.md")
    z <- as.list(setNames(rep(NA_character_,length(keys)),keys)); z$eqtl_full_summary <- NA; return(z)
  }
  if(!is.list(input_manifest)) .cater_stop("input_manifest must be a named list")
  miss <- setdiff(keys[1:10],names(input_manifest))
  if(length(miss)&&isTRUE(strict)) .cater_stop("input_manifest missing required fields: %s",paste(miss,collapse=", "))
  z <- as.list(setNames(rep(NA_character_,length(keys)),keys)); for(k in intersect(names(input_manifest),keys)) z[[k]] <- input_manifest[[k]]
  norm_build <- function(x) toupper(gsub("[^A-Z0-9]","",as.character(x)))
  bs <- vapply(z[c("grn_build","eqtl_build","ld_build")],norm_build,character(1)); bs <- bs[nzchar(bs)&!is.na(bs)]
  if(length(bs)>1L&&length(unique(bs))>1L) .cater_stop("GRN/eQTL/LD genome builds must match")
  if(isTRUE(strict)&&!identical(tolower(as.character(z$eqtl_n_unit)),"donors")) .cater_stop("eQTL N must denote genetically independent donors, not cells")
  if(isTRUE(strict)&&!isTRUE(z$eqtl_full_summary)) .cater_stop("CATER-MR requires full-summary eQTL data")
  if(isTRUE(strict)&&(is.na(z$trans_qtl_qc)||!nzchar(as.character(z$trans_qtl_qc)))) .cater_stop("trans_qtl_qc provenance is required")
  if(isTRUE(strict)&&!identical(tolower(as.character(z$eqtl_ancestry)),tolower(as.character(z$ld_ancestry)))) .cater_stop("COJO LD ancestry must match eQTL ancestry")
  z
}

.cater_primary_decision <- function(fits,net=NULL,target=NULL,has_trans=FALSE,
                                    sibling_screen_performed=FALSE,sibling_screen_complete=TRUE,
                                    n_active_siblings=0L,primary_policy=c("screened_cater","cis_anchor"),
                                    allow_trans_only_primary=FALSE,n_trans_tf_loci=0L,min_trans_tf_loci=3L,
                                    allow_network_primary=FALSE,sibling_screen_independent=FALSE,
                                    sibling_screen_testable=TRUE) {
  primary_policy <- match.arg(primary_policy)
  out <- list(model=NA_character_,beta=NA_real_,se=NA_real_,p=NA_real_,status="MR_FAILED",evidence_status="UNRESOLVED")
  use_fit <- function(model,fit,status,evidence) {
    out$model <<- model; out$beta <<- fit$beta; out$se <<- fit$se; out$p <<- fit$p
    out$status <<- status; out$evidence_status <<- evidence
  }
  has_cis <- identical(fits$cis$status,"OK")
  if(!has_trans) {
    if(has_cis) use_fit("cis",fits$cis,"OK_CIS_ONLY","CIS_ANCHOR")
    return(out)
  }
  if(!isTRUE(sibling_screen_performed)) {
    if(has_cis) {
      use_fit("cis",fits$cis,"TRANS_UNSCREENED_CIS_FALLBACK","TRANS_UNSCREENED")
    } else {
      out$status <- "TRANS_UNSCREENED_NO_CIS"
    }
    return(out)
  }
  if(!isTRUE(sibling_screen_complete)) {
    if(has_cis) {
      use_fit("cis",fits$cis,"SIBLING_SCREEN_INCOMPLETE_CIS_FALLBACK","SIBLING_SCREEN_INCOMPLETE")
    } else {
      out$status <- "SIBLING_SCREEN_INCOMPLETE_NO_CIS"
    }
    return(out)
  }
  if(n_active_siblings>0L) {
    eligible <- primary_policy=="screened_cater" && isTRUE(allow_network_primary) &&
      isTRUE(sibling_screen_independent) && !is.null(net) && identical(net$status,"OK") &&
      isTRUE(net$primary_eligible %||% net$statistical_eligibility)
    if(eligible) {
      idx <- if(!is.null(target)&&target%in%names(net$beta)) target else 1L
      use_fit("network",data.frame(beta=net$beta[idx],se=net$se[idx],p=net$p[idx]),"OK_NETWORK_ADJUSTED","EXPERIMENTAL_NETWORK_PRIMARY")
    } else if(has_cis) {
      use_fit("cis",fits$cis,"TRANS_PLEIOTROPY_UNRESOLVED","MEASURED_CO_PERTURBATION")
    } else {
      out$status <- "TRANS_PLEIOTROPY_UNRESOLVED_NO_CIS"
    }
    return(out)
  }
  if(!has_cis) {
    if(primary_policy=="screened_cater" && isTRUE(allow_trans_only_primary) &&
       n_trans_tf_loci>=min_trans_tf_loci && identical(fits$combined$status,"OK")) {
      use_fit("combined",fits$combined,"OK_CATER_TRANS_ONLY_EXPLORATORY","TRANS_ONLY_MULTI_LOCUS_EXPLORATORY")
    } else {
      out$status <- "TRANS_ONLY_SENSITIVITY"
      out$evidence_status <- if(isTRUE(sibling_screen_testable)) "NO_CIS_ANCHOR" else "NO_TESTABLE_SIBLINGS"
    }
    return(out)
  }
  if(primary_policy=="cis_anchor") {
    use_fit("cis",fits$cis,"OK_CIS_ANCHOR_PRIMARY","CATER_SENSITIVITY_UNCALIBRATED")
  } else if(identical(fits$combined$status,"OK")) {
    use_fit("combined",fits$combined,"OK_CATER","MEASURED_SIBLING_SCREEN_NEGATIVE")
  }
  out
}

# Null-coalescing helper used only for backward-compatible network eligibility.
`%||%` <- function(x,y) if(is.null(x)) y else x
