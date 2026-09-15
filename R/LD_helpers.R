# PLINK signed-LD helpers used by CATER-MR v0.8.

.cater_remove_prefix_outputs <- function(prefix) {
  z <- Sys.glob(paste0(prefix, "*"))
  if (length(z)) unlink(z, recursive=TRUE, force=TRUE)
  invisible(NULL)
}

.cater_quote_args <- function(args) {
  vapply(as.character(args), function(z) if (grepl("[[:space:]]", z)) shQuote(z) else z,
         character(1), USE.NAMES=FALSE)
}

.cater_validate_plink_ld <- function(R, label="PLINK LD", tol=1e-6) {
  R <- as.matrix(R)
  if (!is.numeric(R) || nrow(R) != ncol(R)) .cater_stop("%s must be a square numeric matrix", label)
  if (any(!is.finite(R))) .cater_stop("%s contains non-finite values", label)
  if (max(abs(R-t(R))) > tol) .cater_stop("%s is not symmetric", label)
  if (any(abs(diag(R)-1) > tol)) .cater_stop("%s diagonal is not one", label)
  if (any(abs(R) > 1+tol)) .cater_stop("%s contains correlations outside [-1,1]", label)
  R <- (R+t(R))/2
  diag(R) <- 1
  R
}

.cater_read_plink_bim <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) return(data.frame())
  x <- utils::read.table(path, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE)
  if (!nrow(x)) return(data.frame())
  if (ncol(x) < 6L) .cater_stop("Malformed PLINK BIM: %s", path)
  out <- data.frame(
    snp=as.character(x[[2L]]),
    ld_a1=toupper(as.character(x[[5L]])),
    ld_a2=toupper(as.character(x[[6L]])),
    stringsAsFactors=FALSE
  )
  if (anyDuplicated(out$snp)) .cater_stop("Duplicate SNP IDs in PLINK BIM: %s", path)
  out
}

.cater_filter_nonfinite_ld_rows <- function(R) {
  R <- as.matrix(R)
  if (!nrow(R)) return(list(R=R, nonfinite=character()))
  bad <- which(rowSums(!is.finite(R)) > 0L)
  if (!length(bad)) return(list(R=R, nonfinite=character()))
  ids <- rownames(R)[bad]
  keep <- setdiff(seq_len(nrow(R)), bad)
  list(R=R[keep, keep, drop=FALSE], nonfinite=ids)
}

.cater_read_plink_square_ld <- function(path, snps) {
  if (!length(snps)) return(matrix(numeric(), 0, 0))
  if (length(snps) == 1L) return(matrix(1, 1, 1, dimnames=list(snps, snps)))
  if (!file.exists(path)) .cater_stop("PLINK LD matrix not found: %s", path)
  z <- scan(path, what=double(), quiet=TRUE)
  k <- length(snps)
  if (length(z) != k*k) .cater_stop("PLINK LD matrix has %d values; expected %d for %d SNPs", length(z), k*k, k)
  matrix(z, nrow=k, ncol=k, byrow=TRUE, dimnames=list(snps, snps))
}

.cater_plink_ld <- function(snps, ld_bfile, plink_bin="plink", threads=1L, prefix=tempfile("cater_ld_")) {
  requested <- unique(as.character(snps))
  requested <- requested[nzchar(requested)]
  if (!length(requested)) {
    return(list(R=matrix(numeric(),0,0), ref=data.frame(), ref_all=data.frame(),
                missing=character(), nonfinite=character()))
  }

  exe <- if (length(plink_bin)==1L && file.exists(plink_bin)) normalizePath(plink_bin, mustWork=TRUE) else Sys.which(plink_bin)
  if (!nzchar(exe)) .cater_stop("Cannot find PLINK executable '%s' required for LD selection", plink_bin)

  .cater_remove_prefix_outputs(prefix)
  dir.create(dirname(prefix), recursive=TRUE, showWarnings=FALSE)
  extract <- paste0(prefix, ".extract")
  writeLines(requested, extract)

  refp <- paste0(prefix, ".ref")
  args1 <- c("--bfile", ld_bfile, "--extract", extract, "--keep-allele-order", "--make-bed",
             "--threads", as.character(as.integer(threads)), "--out", refp)
  st <- system2(exe, args=.cater_quote_args(args1), stdout=paste0(refp,".stdout"), stderr=paste0(refp,".stderr"))
  if (!identical(st, 0L)) {
    logtxt <- if (file.exists(paste0(refp,".log"))) paste(readLines(paste0(refp,".log"), warn=FALSE), collapse="\n") else ""
    if (grepl("No variants remaining", logtxt, fixed=TRUE)) {
      return(list(R=matrix(numeric(),0,0), ref=data.frame(), ref_all=data.frame(),
                  missing=requested, nonfinite=character()))
    }
    .cater_stop("PLINK subset generation failed for LD selection: %s", prefix)
  }

  ref_all <- .cater_read_plink_bim(paste0(refp, ".bim"))
  if (!nrow(ref_all)) {
    return(list(R=matrix(numeric(),0,0), ref=data.frame(), ref_all=ref_all,
                missing=requested, nonfinite=character()))
  }
  missing <- setdiff(requested, ref_all$snp)

  if (nrow(ref_all) == 1L) {
    R <- matrix(1, 1, 1, dimnames=list(ref_all$snp, ref_all$snp))
  } else {
    ldp <- paste0(prefix, ".r")
    args2 <- c("--bfile", refp, "--keep-allele-order", "--r", "square",
               "--threads", as.character(as.integer(threads)), "--out", ldp)
    st2 <- system2(exe, args=.cater_quote_args(args2), stdout=paste0(ldp,".stdout"), stderr=paste0(ldp,".stderr"))
    if (!identical(st2, 0L)) .cater_stop("PLINK signed LD calculation failed for LD selection: %s", prefix)
    R <- .cater_read_plink_square_ld(paste0(ldp, ".ld"), ref_all$snp)
  }

  finite_filter <- .cater_filter_nonfinite_ld_rows(R)
  R <- finite_filter$R
  nonfinite <- finite_filter$nonfinite
  ref <- ref_all[match(rownames(R), ref_all$snp),,drop=FALSE]
  if (nrow(R)) R <- .cater_validate_plink_ld(R, "PLINK signed LD")

  list(R=R, ref=ref, ref_all=ref_all, missing=missing, nonfinite=unique(nonfinite))
}
