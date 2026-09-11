from pathlib import Path
import re

# 1) LD diagnostic hardening: unresolved palindromic alleles are not signed for
# SuSiE-RSS, and all rows touched by non-finite LD are removed simultaneously.
p = Path("R/LD_diagnosis_precojo.R")
s = p.read_text()

old = '''  same <- qa1==ra1 & qa2==ra2
  swap <- qa1==ra2 & qa2==ra1
  ca1 <- .cater_complement_allele(qa1); ca2 <- .cater_complement_allele(qa2)
  comp_same <- !same & !swap & ca1==ra1 & ca2==ra2
  comp_swap <- !same & !swap & ca1==ra2 & ca2==ra1
  # Keep the same direct A1/A2 contract as the downstream COJO/MR alignment.
  # Strand complements are recorded for audit, but are not inferred as matches.
  sign <- ifelse(same, 1, ifelse(swap, -1, NA_real_))
  alignment <- ifelse(same, "same", ifelse(swap, "swap",
    ifelse(comp_same, "strand_same", ifelse(comp_swap, "strand_swap", "mismatch"))))
  z0 <- d$beta/d$se
  data.frame(snp=d$snp,qtl_a1=qa1,qtl_a2=qa2,ld_a1=ra1,ld_a2=ra2,
             z_raw=z0,z=z0*sign,alignment=alignment,allele_match=same|swap,
             stringsAsFactors=FALSE)
'''
new = '''  same <- qa1==ra1 & qa2==ra2
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
'''
assert old in s, "alignment block not found"
s = s.replace(old, new, 1)

needle = '''.cater_read_plink_square_ld <- function(path, snps) {
'''
helper = '''.cater_filter_nonfinite_ld_rows <- function(R) {
  R <- as.matrix(R)
  if (!nrow(R)) return(list(R=R,nonfinite=character()))
  bad <- which(rowSums(!is.finite(R)) > 0L)
  if (!length(bad)) return(list(R=R,nonfinite=character()))
  ids <- rownames(R)[bad]
  keep <- setdiff(seq_len(nrow(R)),bad)
  list(R=R[keep,keep,drop=FALSE],nonfinite=ids)
}

'''
assert needle in s, "square-LD function anchor not found"
s = s.replace(needle, helper + needle, 1)

old = '''  nonfinite <- character()
  while (nrow(R) && any(!is.finite(R))) {
    bad_diag <- which(!is.finite(diag(R)))
    bad <- if (length(bad_diag)) bad_diag[1L] else which.max(rowSums(!is.finite(R)))
    nonfinite <- c(nonfinite,rownames(R)[bad])
    keep <- setdiff(seq_len(nrow(R)),bad)
    R <- R[keep,keep,drop=FALSE]
  }
  ref <- ref_all[match(rownames(R),ref_all$snp),,drop=FALSE]
'''
new = '''  finite_filter <- .cater_filter_nonfinite_ld_rows(R)
  R <- finite_filter$R
  nonfinite <- finite_filter$nonfinite
  ref <- ref_all[match(rownames(R),ref_all$snp),,drop=FALSE]
'''
assert old in s, "greedy non-finite loop not found"
s = s.replace(old, new, 1)

old = '''      badal <- jj[!al$allele_match]
      if (length(badal)) { rows$status[badal] <- "LD_ALLELE_MISMATCH"; rows$remove[badal] <- TRUE }
      good <- al$snp[al$allele_match]
'''
new = '''      pal <- jj[al$palindromic]
      if (length(pal)) {
        rows$status[pal] <- "NOT_DIAGNOSABLE_PALINDROMIC"
        rows$remove[pal] <- FALSE
      }
      badal <- jj[!al$allele_match & !al$palindromic]
      if (length(badal)) { rows$status[badal] <- "LD_ALLELE_MISMATCH"; rows$remove[badal] <- TRUE }
      good <- al$snp[al$allele_match]
'''
assert old in s, "allele status block not found"
s = s.replace(old, new, 1)
p.write_text(s)

# 2) Strict logical validation.
p = Path("R/CATER_MR_core_v05.R")
s = p.read_text()
old = '  if(length(enable_ld_diagnosis)!=1L||is.na(enable_ld_diagnosis)) .cater_stop("enable_ld_diagnosis must be TRUE or FALSE")\n'
new = '  if(!is.logical(enable_ld_diagnosis)||length(enable_ld_diagnosis)!=1L||is.na(enable_ld_diagnosis)) .cater_stop("enable_ld_diagnosis must be TRUE or FALSE")\n'
assert old in s, "enable_ld_diagnosis validation not found"
s = s.replace(old, new, 1)
p.write_text(s)

# 3) Regression tests for all remaining Codex findings.
p = Path("tests/v06_hardening_smoke.R")
s = p.read_text()
anchor = 'cat("CATER-MR direct-core evidence-safety tests passed\\n")\n'
assert anchor in s, "test footer not found"
extra = r'''# 20. Palindromic A/T and C/G variants are never signed from allele labels alone.
qp <- data.frame(snp=c("p1","p2","n1"),a1=c("A","C","A"),a2=c("T","G","G"),
                 beta=c(.2,.2,.2),se=c(.1,.1,.1),stringsAsFactors=FALSE)
rp <- data.frame(snp=c("p1","p2","n1"),ld_a1=c("T","C","G"),ld_a2=c("A","G","A"),
                 stringsAsFactors=FALSE)
ap <- .cater_align_z_to_plink(qp,rp)
stopifnot(identical(ap$palindromic,c(TRUE,TRUE,FALSE)),
          identical(ap$alignment,c("palindromic_unresolved","palindromic_unresolved","swap")),
          all(is.na(ap$z[1:2])),isTRUE(all.equal(ap$z[3],-2,tolerance=1e-12)),
          identical(ap$allele_match,c(FALSE,FALSE,TRUE)))

# 21. Every variant whose original LD row contains a non-finite entry is removed together.
Rnf <- matrix(c(1,NaN,0.2,NaN,1,0.3,0.2,0.3,1),3,3,byrow=TRUE,
              dimnames=list(c("a","b","c"),c("a","b","c")))
fnf <- .cater_filter_nonfinite_ld_rows(Rnf)
stopifnot(identical(sort(fnf$nonfinite),c("a","b")),
          identical(dim(fnf$R),c(1L,1L)),rownames(fnf$R)=="c")
Rnf2 <- matrix(c(1,NaN,NaN,1),2,2,byrow=TRUE,
               dimnames=list(c("x","y"),c("x","y")))
fnf2 <- .cater_filter_nonfinite_ld_rows(Rnf2)
stopifnot(setequal(fnf2$nonfinite,c("x","y")),nrow(fnf2$R)==0L)

# 22. The diagnosis switch must be a scalar logical; numeric/string truthy values fail closed.
for (bad_flag in list(1,"TRUE")) {
  err <- tryCatch({
    cater_mr(grn=data.frame(TF="A",Target="B"),eqtl_dir=tempdir(),outcome=data.frame(),
             ld_bfile="unused",enable_ld_diagnosis=bad_flag)
    NULL
  },error=function(e)e)
  stopifnot(inherits(err,"error"),grepl("enable_ld_diagnosis must be TRUE or FALSE",conditionMessage(err),fixed=TRUE))
}

'''
s = s.replace(anchor, extra + anchor, 1)
p.write_text(s)

# 4) Documentation of the conservative palindromic handling and fail-closed LD rows.
p = Path("LD_DIAGNOSIS.md")
s = p.read_text()
needle = "Variants absent from the LD reference, variants with non-finite LD rows, and variants whose allele pair cannot be reconciled with the direct A1/A2 contract are also removed and explicitly reported. A one-variant locus cannot be conditionally diagnosed and is retained with status `NOT_DIAGNOSABLE_SINGLETON`.\n"
replacement = "Variants absent from the LD reference, variants with non-finite LD rows, and variants whose allele pair cannot be reconciled with the direct A1/A2 contract are also removed and explicitly reported. All variants whose original PLINK LD row contains any non-finite entry are removed simultaneously, so this QC is not variant-order dependent. Palindromic A/T and C/G variants are not assigned a signed z score from allele labels alone; without frequency-based strand disambiguation they are excluded from SuSiE-RSS diagnosis, retained for downstream COJO, and reported as `NOT_DIAGNOSABLE_PALINDROMIC`. A one-variant diagnosable locus is retained with status `NOT_DIAGNOSABLE_SINGLETON`.\n"
assert needle in s, "documentation paragraph not found"
s = s.replace(needle,replacement,1)
p.write_text(s)
