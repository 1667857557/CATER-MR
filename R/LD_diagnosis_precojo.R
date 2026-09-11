# Pre-COJO LD consistency diagnostics for CATER-MR.
# This file defines helpers only; core CATER-MR functions are modified directly
# to call these helpers rather than being overridden at runtime.

.cater_complement_allele <- function(x) {
  x <- toupper(as.character(x))
  out <- chartr("ACGT", "TGCA", x)
  out[!x %in% c("A","C","G","T")] <- NA_character_
  out
}

.cater_align_z_to_plink <- function(qtl, ref) {
  if (!nrow(ref)) return(data.frame())
  d <- qtl[match(ref$snp, qtl$snp),,drop=FALSE]
  if (any(is.na(d$snp))) .cater_stop("PLINK LD reference contains SNPs absent from eQTL summary data")
  qa1 <- toupper(d$a1); qa2 <- toupper(d$a2)
  ra1 <- toupper(ref$ld_a1); ra2 <- toupper(ref$ld_a2)
  same <- qa1==ra1 & qa2==ra2
  swap <- qa1==ra2 & qa2==ra1
  palindromic <- paste0(qa1,qa2) %in% c("AT","TA","CG","GC")
  ca1 <- .cater_complement_allele(qa1); ca2 <- .cater_complement_allele(qa2)
  comp_same <- !same & !swap & ca1==ra1 & ca2==ra2
  comp_swap <- !same & !swap & ca1==ra2 & ca2==ra1
  # A/T and C/G variants cannot be strand-resolved from allele labels alone.
  # Do not infer their sign for SuSiE-RSS without a frequency-based contract.
  direct_match <- (same|swap) & !palindromic
  sign <- ifelse(palindromic, NA_real_, ifelse(same, 1, ifelse(swap, -1, NA_real_)))
  alignment <- ifelse(palindromic, "palindromic_unresolved",
    ifelse(same, "same", ifelse(swap, "swap",
      ifelse(comp_same, "strand_same", ifelse(comp_swap, "strand_swap", "mismatch")))))
  z0 <- d$beta/d$se
  data.frame(snp=d$snp,qtl_a1=qa1,qtl_a2=qa2,ld_a1=ra1,ld_a2=ra2,
             z_raw=z0,z=z0*sign,alignment=alignment,allele_match=direct_match,
             palindromic=palindromic,stringsAsFactors=FALSE)
}

.cater_ld_outlier_index <- function(conditional_dist, loglr_cutoff=2, abs_z_cutoff=2) {
  if (!is.data.frame(conditional_dist)) conditional_dist <- as.data.frame(conditional_dist)
  if (!all(c("logLR","z") %in% names(conditional_dist)))
    .cater_stop("SuSiE-RSS conditional distribution must contain logLR and z")
  which(is.finite(conditional_dist$logLR) & conditional_dist$logLR > loglr_cutoff &
          is.finite(conditional_dist$z) & abs(conditional_dist$z) > abs_z_cutoff)
}

.cater_ld_diagnosis_groups <- function(q, candidate_map) {
  if (nrow(q)!=nrow(candidate_map)) .cater_stop("LD diagnosis grouping requires aligned qtl and candidate_map rows")
  if (!nrow(candidate_map)) return(character())
  if (any(is.na(candidate_map$locus_id) | !nzchar(candidate_map$locus_id)))
    .cater_stop("LD diagnosis requires one physical locus_id per candidate")
  paste(q$chr,candidate_map$locus_id,sep="|")
}

.cater_validate_diag_ld <- function(R, label="diagnostic LD", tol=1e-6) {
  R <- as.matrix(R)
  if (!is.numeric(R) || nrow(R)!=ncol(R)) .cater_stop("%s must be a square numeric matrix",label)
  if (any(!is.finite(R))) .cater_stop("%s contains non-finite values",label)
  if (max(abs(R-t(R)))>tol) .cater_stop("%s is not symmetric",label)
  if (any(abs(diag(R)-1)>tol)) .cater_stop("%s diagonal is not one",label)
  if (any(abs(R)>1+tol)) .cater_stop("%s contains correlations outside [-1,1]",label)
  # Do not impose an extra PSD eigenvalue gate here. susieR performs its own
  # RSS eigenvalue handling; an additional strict gate can reject rounded PLINK
  # correlation matrices before the reference diagnostic gets a chance to run.
  R <- (R+t(R))/2
  diag(R) <- 1
  R
}

.cater_read_plink_bim <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size==0) return(data.frame())
  x <- utils::read.table(path, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE)
  if (!nrow(x)) return(data.frame())
  if (ncol(x) < 6L) .cater_stop("Malformed PLINK BIM: %s", path)
  out <- data.frame(snp=as.character(x[[2L]]),ld_a1=toupper(as.character(x[[5L]])),
                    ld_a2=toupper(as.character(x[[6L]])),stringsAsFactors=FALSE)
  if (anyDuplicated(out$snp)) .cater_stop("Duplicate SNP IDs in PLINK BIM: %s", path)
  out
}

.cater_filter_nonfinite_ld_rows <- function(R) {
  R <- as.matrix(R)
  if (!nrow(R)) return(list(R=R,nonfinite=character()))
  bad <- which(rowSums(!is.finite(R)) > 0L)
  if (!length(bad)) return(list(R=R,nonfinite=character()))
  ids <- rownames(R)[bad]
  keep <- setdiff(seq_len(nrow(R)),bad)
  list(R=R[keep,keep,drop=FALSE],nonfinite=ids)
}

.cater_read_plink_square_ld <- function(path, snps) {
  if (!length(snps)) return(matrix(numeric(),0,0))
  if (length(snps)==1L) return(matrix(1,1,1,dimnames=list(snps,snps)))
  if (!file.exists(path)) .cater_stop("PLINK LD matrix not found: %s", path)
  z <- scan(path, what=double(), quiet=TRUE)
  k <- length(snps)
  if (length(z) != k*k) .cater_stop("PLINK LD matrix has %d values; expected %d for %d SNPs",length(z),k*k,k)
  matrix(z,nrow=k,ncol=k,byrow=TRUE,dimnames=list(snps,snps))
}

.cater_plink_ld <- function(snps, ld_bfile, plink_bin, threads, prefix) {
  requested <- unique(as.character(snps)); requested <- requested[nzchar(requested)]
  if (!length(requested)) return(list(R=matrix(numeric(),0,0),ref=data.frame(),ref_all=data.frame(),missing=character(),nonfinite=character()))
  exe <- if (length(plink_bin)==1L && file.exists(plink_bin)) normalizePath(plink_bin,mustWork=TRUE) else Sys.which(plink_bin)
  if (!nzchar(exe)) .cater_stop("Cannot find PLINK executable '%s' required for pre-COJO LD diagnosis",plink_bin)
  .cater_remove_prefix_outputs(prefix)
  extract <- paste0(prefix,".extract"); writeLines(requested,extract)
  refp <- paste0(prefix,".ref")
  args1 <- c("--bfile",ld_bfile,"--extract",extract,"--keep-allele-order","--make-bed",
             "--threads",as.character(as.integer(threads)),"--out",refp)
  st <- system2(exe,args=.cater_quote_args(args1),stdout=paste0(refp,".stdout"),stderr=paste0(refp,".stderr"))
  if (!identical(st,0L)) {
    logtxt <- if (file.exists(paste0(refp,".log"))) paste(readLines(paste0(refp,".log"),warn=FALSE),collapse="\n") else ""
    if (grepl("No variants remaining",logtxt,fixed=TRUE))
      return(list(R=matrix(numeric(),0,0),ref=data.frame(),ref_all=data.frame(),missing=requested,nonfinite=character()))
    .cater_stop("PLINK subset generation failed for LD diagnosis: %s",prefix)
  }
  ref_all <- .cater_read_plink_bim(paste0(refp,".bim"))
  if (!nrow(ref_all)) return(list(R=matrix(numeric(),0,0),ref=data.frame(),ref_all=ref_all,missing=requested,nonfinite=character()))
  missing <- setdiff(requested,ref_all$snp)
  if (nrow(ref_all)==1L) {
    R <- matrix(1,1,1,dimnames=list(ref_all$snp,ref_all$snp))
  } else {
    ldp <- paste0(prefix,".r")
    args2 <- c("--bfile",refp,"--keep-allele-order","--r","square",
               "--threads",as.character(as.integer(threads)),"--out",ldp)
    st2 <- system2(exe,args=.cater_quote_args(args2),stdout=paste0(ldp,".stdout"),stderr=paste0(ldp,".stderr"))
    if (!identical(st2,0L)) .cater_stop("PLINK signed LD calculation failed for LD diagnosis: %s",prefix)
    R <- .cater_read_plink_square_ld(paste0(ldp,".ld"),ref_all$snp)
  }
  finite_filter <- .cater_filter_nonfinite_ld_rows(R)
  R <- finite_filter$R
  nonfinite <- finite_filter$nonfinite
  ref <- ref_all[match(rownames(R),ref_all$snp),,drop=FALSE]
  if (nrow(R)) R <- .cater_validate_diag_ld(R,"pre-COJO PLINK signed LD")
  list(R=R,ref=ref,ref_all=ref_all,missing=missing,nonfinite=unique(nonfinite))
}

.cater_susie_ld_diagnosis <- function(z, R, n=NULL, loglr_cutoff=2, abs_z_cutoff=2) {
  if (!requireNamespace("susieR",quietly=TRUE))
    .cater_stop("Package 'susieR' is required when enable_ld_diagnosis=TRUE")
  z <- as.numeric(z); R <- .cater_validate_diag_ld(R,"SuSiE-RSS diagnostic LD")
  if (length(z)!=nrow(R) || any(!is.finite(z))) .cater_stop("SuSiE-RSS LD diagnosis requires one finite z score per LD variant")
  args <- list(z=z,R=R)
  if (!is.null(n) && length(n)==1L && is.finite(n) && n>0) args$n <- as.numeric(n)
  lambda <- tryCatch(do.call(susieR::estimate_s_rss,args),error=function(e)
    .cater_stop("SuSiE-RSS estimate_s_rss failed: %s",conditionMessage(e)))
  kres <- tryCatch(do.call(susieR::kriging_rss,c(args,list(s=lambda))),error=function(e)
    .cater_stop("SuSiE-RSS kriging_rss failed: %s",conditionMessage(e)))
  cd <- as.data.frame(kres$conditional_dist)
  if (nrow(cd)!=length(z)) .cater_stop("SuSiE-RSS returned an unexpected conditional distribution size")
  if ("z_std_diff" %in% names(cd)) cd$p_diff <- stats::pchisq(cd$z_std_diff^2,df=1,lower.tail=FALSE) else cd$p_diff <- NA_real_
  cd$detected <- FALSE
  ii <- .cater_ld_outlier_index(cd,loglr_cutoff,abs_z_cutoff)
  if (length(ii)) cd$detected[ii] <- TRUE
  list(lambda=as.numeric(lambda)[1L],conditional_dist=cd)
}

.cater_ld_diagnosis_filter <- function(qtl,candidate_map,ld_bfile,plink_bin,threads,prefix,
                                       loglr_cutoff=2,abs_z_cutoff=2,verbose=TRUE) {
  if (!nrow(candidate_map)) return(list(candidate_map=candidate_map,diagnostics=data.frame(),removed=character()))
  q <- qtl[match(candidate_map$snp,qtl$snp),,drop=FALSE]
  if (any(is.na(q$snp))) .cater_stop("Candidate map contains SNPs absent from eQTL summary data")
  groups <- .cater_ld_diagnosis_groups(q,candidate_map)
  tabs <- list(); ug <- unique(groups)
  for (gidx in seq_along(ug)) {
    ii <- which(groups==ug[gidx]); ids <- candidate_map$snp[ii]; qq <- q[ii,,drop=FALSE]
    cm <- candidate_map[ii,,drop=FALSE]
    .cater_msg(verbose,"LD diagnosis %d/%d: %s (%d SNPs)",gidx,length(ug),ug[gidx],length(ids))
    pl <- .cater_plink_ld(ids,ld_bfile,plink_bin,threads,paste0(prefix,".lddiag.",gidx))
    rows <- data.frame(snp=ids,chr=qq$chr,locus_id=cm$locus_id,source=cm$source,
      status=NA_character_,remove=FALSE,qtl_a1=qq$a1,qtl_a2=qq$a2,ld_a1=NA_character_,ld_a2=NA_character_,
      z_raw=qq$beta/qq$se,z=NA_real_,alignment=NA_character_,condmean=NA_real_,condvar=NA_real_,
      z_std_diff=NA_real_,logLR=NA_real_,p_diff=NA_real_,lambda=NA_real_,n_used=NA_real_,stringsAsFactors=FALSE)
    if (nrow(pl$ref_all)) {
      jj <- match(rows$snp,pl$ref_all$snp); has <- !is.na(jj)
      rows$ld_a1[has] <- pl$ref_all$ld_a1[jj[has]]; rows$ld_a2[has] <- pl$ref_all$ld_a2[jj[has]]
    }
    miss <- rows$snp %in% pl$missing
    rows$status[miss] <- "LD_REFERENCE_MISSING"; rows$remove[miss] <- TRUE
    nf <- rows$snp %in% pl$nonfinite
    rows$status[nf] <- "LD_NONFINITE"; rows$remove[nf] <- TRUE
    if (nrow(pl$ref)) {
      al <- .cater_align_z_to_plink(qq,pl$ref)
      jj <- match(al$snp,rows$snp)
      rows$z[jj] <- al$z; rows$alignment[jj] <- al$alignment
      pal <- jj[al$palindromic]
      if (length(pal)) {
        rows$status[pal] <- "NOT_DIAGNOSABLE_PALINDROMIC"
        rows$remove[pal] <- FALSE
      }
      badal <- jj[!al$allele_match & !al$palindromic]
      if (length(badal)) { rows$status[badal] <- "LD_ALLELE_MISMATCH"; rows$remove[badal] <- TRUE }
      good <- al$snp[al$allele_match]
      if (length(good)) {
        nvals <- qq$n[match(good,qq$snp)]; nvals <- nvals[is.finite(nvals)&nvals>0]
        n_used <- if(length(nvals)) stats::median(nvals) else NULL
        gj <- match(good,rows$snp); if(!is.null(n_used)) rows$n_used[gj] <- n_used
        if (length(good)==1L) {
          rows$status[gj] <- "NOT_DIAGNOSABLE_SINGLETON"
        } else {
          Rg <- pl$R[good,good,drop=FALSE]
          zg <- al$z[match(good,al$snp)]
          su <- .cater_susie_ld_diagnosis(zg,Rg,n=n_used,loglr_cutoff=loglr_cutoff,abs_z_cutoff=abs_z_cutoff)
          cd <- su$conditional_dist
          getcol <- function(nm) if(nm %in% names(cd)) suppressWarnings(as.numeric(cd[[nm]])) else rep(NA_real_,nrow(cd))
          rows$condmean[gj] <- getcol("condmean"); rows$condvar[gj] <- getcol("condvar")
          rows$z_std_diff[gj] <- getcol("z_std_diff"); rows$logLR[gj] <- getcol("logLR")
          rows$p_diff[gj] <- getcol("p_diff"); rows$lambda[gj] <- su$lambda
          det <- as.logical(cd$detected)
          rows$status[gj] <- ifelse(det,"SUSIE_LD_INCONSISTENT","OK")
          rows$remove[gj] <- det
        }
      }
    }
    if (any(is.na(rows$status))) .cater_stop("Internal error: incomplete pre-COJO LD diagnosis status for %s",ug[gidx])
    tabs[[gidx]] <- rows
  }
  diag <- do.call(rbind,tabs); rownames(diag) <- NULL
  diag <- diag[match(candidate_map$snp,diag$snp),,drop=FALSE]
  removed <- unique(diag$snp[diag$remove])
  list(candidate_map=candidate_map[!candidate_map$snp %in% removed,,drop=FALSE],diagnostics=diag,removed=removed)
}
