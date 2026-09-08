# CATER-MR v0.1
# Cis And Trans eQTLs guided by Regulatory networks for drug-target MR
# Minimal, target-centric implementation.

.cater_stop <- function(...) stop(sprintf(...), call. = FALSE)
.cater_msg <- function(verbose, ...) if (isTRUE(verbose)) message(sprintf(...))

.cater_pick_col <- function(x, candidates, required = TRUE, label = NULL) {
  nm <- names(x)
  hit <- match(tolower(candidates), tolower(nm), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(nm[hit[1L]])
  if (required) .cater_stop("Missing required column%s. Tried: %s",
                            if (is.null(label)) "" else paste0(" for ", label),
                            paste(candidates, collapse = ", "))
  NULL
}

.cater_read_table <- function(path) {
  if (!file.exists(path)) .cater_stop("File does not exist: %s", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  utils::read.delim(con, header = TRUE, sep = "\t", quote = "", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE)
}

.cater_standardize_sumstats <- function(x, n_default = NULL, label = "summary data", require_position = TRUE) {
  if (!is.data.frame(x)) .cater_stop("%s must be a data.frame", label)
  snp <- .cater_pick_col(x, c("SNP", "rsid", "variant_id", "variant", "ID"), label = label)
  chr <- .cater_pick_col(x, c("CHR", "chrom", "chromosome"), required = require_position, label = label)
  pos <- .cater_pick_col(x, c("BP", "POS", "position", "base_pair_location"), required = require_position, label = label)
  a1  <- .cater_pick_col(x, c("A1", "EA", "effect_allele", "ALT"), label = label)
  a2  <- .cater_pick_col(x, c("A2", "NEA", "other_allele", "non_effect_allele", "REF"), label = label)
  b   <- .cater_pick_col(x, c("b", "beta", "BETA", "effect", "estimate"), label = label)
  se  <- .cater_pick_col(x, c("se", "SE", "stderr", "standard_error"), label = label)
  p   <- .cater_pick_col(x, c("p", "P", "pval", "p_value", "pvalue"), required = FALSE)
  eaf <- .cater_pick_col(x, c("freq", "EAF", "eaf", "effect_allele_frequency", "AF"), required = FALSE)
  n   <- .cater_pick_col(x, c("N", "n", "samplesize", "sample_size"), required = FALSE)

  out <- data.frame(
    snp = as.character(x[[snp]]),
    chr = if (is.null(chr)) rep(NA_character_, nrow(x)) else sub("^chr", "", as.character(x[[chr]]), ignore.case = TRUE),
    pos = if (is.null(pos)) rep(NA_real_, nrow(x)) else as.numeric(x[[pos]]),
    a1 = toupper(as.character(x[[a1]])),
    a2 = toupper(as.character(x[[a2]])),
    beta = as.numeric(x[[b]]),
    se = as.numeric(x[[se]]),
    stringsAsFactors = FALSE
  )
  out$p <- if (is.null(p)) 2 * stats::pnorm(-abs(out$beta / out$se)) else as.numeric(x[[p]])
  out$eaf <- if (is.null(eaf)) NA_real_ else as.numeric(x[[eaf]])
  out$n <- if (is.null(n)) {
    if (is.null(n_default)) NA_real_ else rep(as.numeric(n_default), nrow(x))
  } else as.numeric(x[[n]])

  keep <- !is.na(out$snp) & nzchar(out$snp) &
    (!require_position | is.finite(out$pos)) &
    is.finite(out$beta) & is.finite(out$se) & out$se > 0 &
    is.finite(out$p) & out$p >= 0 & out$p <= 1 &
    out$a1 %in% c("A", "C", "G", "T") & out$a2 %in% c("A", "C", "G", "T")
  out <- out[keep, , drop = FALSE]
  out <- out[!duplicated(out$snp), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.cater_standardize_grn <- function(grn) {
  if (is.list(grn) && !is.data.frame(grn) && !is.null(grn$grn)) grn <- grn$grn
  if (!is.data.frame(grn)) .cater_stop("grn must be a data.frame or a list containing $grn")
  tf <- .cater_pick_col(grn, c("TF", "tf", "regulator"), label = "GRN TF")
  target <- .cater_pick_col(grn, c("Target", "target", "gene"), label = "GRN target")
  out <- data.frame(TF = as.character(grn[[tf]]), Target = as.character(grn[[target]]),
                    stringsAsFactors = FALSE)
  out <- out[nzchar(out$TF) & nzchar(out$Target) & !is.na(out$TF) & !is.na(out$Target), , drop = FALSE]
  unique(out)
}

.cater_annotation_from_grn <- function(grn_raw) {
  if (is.list(grn_raw) && !is.data.frame(grn_raw) && !is.null(grn_raw$grn)) grn_raw <- grn_raw$grn
  if (!is.data.frame(grn_raw)) return(NULL)
  tf_chr <- .cater_pick_col(grn_raw, c("TF_chr", "tf_chr"), required = FALSE)
  tf_tss <- .cater_pick_col(grn_raw, c("TF_tss", "tf_tss"), required = FALSE)
  tg_chr <- .cater_pick_col(grn_raw, c("Target_chr", "target_chr"), required = FALSE)
  tg_tss <- .cater_pick_col(grn_raw, c("Target_tss", "target_tss"), required = FALSE)
  tf <- .cater_pick_col(grn_raw, c("TF", "tf", "regulator"), required = FALSE)
  tg <- .cater_pick_col(grn_raw, c("Target", "target", "gene"), required = FALSE)
  if (any(vapply(list(tf_chr, tf_tss, tg_chr, tg_tss, tf, tg), is.null, logical(1)))) return(NULL)
  a <- rbind(
    data.frame(symbol = as.character(grn_raw[[tf]]), chr = as.character(grn_raw[[tf_chr]]), tss = as.numeric(grn_raw[[tf_tss]])),
    data.frame(symbol = as.character(grn_raw[[tg]]), chr = as.character(grn_raw[[tg_chr]]), tss = as.numeric(grn_raw[[tg_tss]]))
  )
  a$chr <- sub("^chr", "", a$chr, ignore.case = TRUE)
  a <- a[!is.na(a$symbol) & nzchar(a$symbol) & is.finite(a$tss), , drop = FALSE]
  a[!duplicated(a$symbol), , drop = FALSE]
}

.cater_standardize_annotation <- function(annotation) {
  if (is.null(annotation)) return(NULL)
  if (!is.data.frame(annotation)) .cater_stop("gene_annotation must be a data.frame")
  sym <- .cater_pick_col(annotation, c("symbol", "gene", "gene_symbol", "SYMBOL"), label = "gene annotation symbol")
  chr <- .cater_pick_col(annotation, c("chr", "CHR", "chrom", "chromosome"), label = "gene annotation chromosome")
  tss <- .cater_pick_col(annotation, c("tss", "TSS", "txStart"), label = "gene annotation TSS")
  out <- data.frame(symbol = as.character(annotation[[sym]]),
                    chr = sub("^chr", "", as.character(annotation[[chr]]), ignore.case = TRUE),
                    tss = as.numeric(annotation[[tss]]), stringsAsFactors = FALSE)
  out <- out[!is.na(out$symbol) & nzchar(out$symbol) & is.finite(out$tss), , drop = FALSE]
  out[!duplicated(out$symbol), , drop = FALSE]
}

.cater_make_regions <- function(target, parents, annotation, cis_window, tf_window) {
  ta <- annotation[annotation$symbol == target, , drop = FALSE]
  if (!nrow(ta)) return(NULL)
  reg <- data.frame(type = "cis", gene = target, chr = ta$chr[1],
                    start = max(1, ta$tss[1] - cis_window), end = ta$tss[1] + cis_window,
                    stringsAsFactors = FALSE)
  if (length(parents)) {
    pa <- annotation[match(parents, annotation$symbol), , drop = FALSE]
    pa$gene <- parents
    pa <- pa[!is.na(pa$chr) & is.finite(pa$tss), , drop = FALSE]
    if (nrow(pa)) {
      treg <- data.frame(type = "trans", gene = pa$gene, chr = pa$chr,
                         start = pmax(1, pa$tss - tf_window), end = pa$tss + tf_window,
                         stringsAsFactors = FALSE)
      reg <- rbind(reg, treg)
    }
  }
  rownames(reg) <- NULL
  reg
}

.cater_candidate_map <- function(qtl, regions, target) {
  if (is.null(regions) || !nrow(qtl)) return(data.frame())
  cis <- regions[regions$type == "cis", , drop = FALSE]
  is_cis <- qtl$chr == cis$chr[1] & qtl$pos >= cis$start[1] & qtl$pos <= cis$end[1]
  trans_regs <- regions[regions$type == "trans", , drop = FALSE]
  parent_hits <- vector("list", nrow(qtl))
  if (nrow(trans_regs)) {
    for (k in seq_len(nrow(trans_regs))) {
      hit <- qtl$chr == trans_regs$chr[k] & qtl$pos >= trans_regs$start[k] & qtl$pos <= trans_regs$end[k] & !is_cis
      if (any(hit)) parent_hits[hit] <- lapply(parent_hits[hit], function(z) c(z, trans_regs$gene[k]))
    }
  }
  is_trans <- lengths(parent_hits) > 0L
  keep <- is_cis | is_trans
  if (!any(keep)) return(data.frame())
  data.frame(
    snp = qtl$snp[keep],
    source = ifelse(is_cis[keep], "cis", "trans"),
    parent_tf = vapply(parent_hits[keep], function(z) if (length(z)) paste(unique(z), collapse = ";") else "", character(1)),
    stringsAsFactors = FALSE
  )
}

.cater_write_cojo_ma <- function(qtl, path) {
  ok <- is.finite(qtl$eaf) & qtl$eaf > 0 & qtl$eaf < 1 & is.finite(qtl$n) & qtl$n > 0
  if (!all(ok)) {
    .cater_stop("COJO requires EAF/freq and N for all retained QTL rows; %d rows are missing/invalid", sum(!ok))
  }
  ma <- data.frame(SNP = qtl$snp, A1 = qtl$a1, A2 = qtl$a2, freq = qtl$eaf,
                   b = qtl$beta, se = qtl$se, p = qtl$p, N = qtl$n,
                   stringsAsFactors = FALSE)
  utils::write.table(ma, path, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

.cater_run_cojo <- function(qtl, candidate_map, ld_bfile, gcta_bin, cojo_p,
                            cojo_wind_kb, cojo_collinear, prefix, verbose) {
  if (!length(ld_bfile) || !nzchar(ld_bfile)) .cater_stop("ld_bfile is required for COJO")
  if (!file.exists(paste0(ld_bfile, ".bed")) || !file.exists(paste0(ld_bfile, ".bim")) || !file.exists(paste0(ld_bfile, ".fam"))) {
    .cater_stop("ld_bfile must point to a PLINK bed/bim/fam prefix: %s", ld_bfile)
  }
  exe <- Sys.which(gcta_bin)
  if (!nzchar(exe)) .cater_stop("Cannot find GCTA executable '%s' in PATH", gcta_bin)
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  ma <- paste0(prefix, ".ma")
  extract <- paste0(prefix, ".candidate.snplist")
  .cater_write_cojo_ma(qtl, ma)
  writeLines(unique(candidate_map$snp), extract)
  args <- c("--bfile", ld_bfile,
            "--cojo-file", ma,
            "--extract", extract,
            "--cojo-slct",
            "--cojo-p", format(cojo_p, scientific = TRUE),
            "--cojo-wind", as.character(as.integer(cojo_wind_kb)),
            "--cojo-collinear", as.character(cojo_collinear),
            "--out", prefix)
  .cater_msg(verbose, "COJO: %s %s", exe, paste(args, collapse = " "))
  status <- system2(exe, args = args, stdout = paste0(prefix, ".stdout"), stderr = paste0(prefix, ".stderr"))
  if (!identical(status, 0L)) .cater_stop("GCTA-COJO failed for %s (status %s); see %s.stderr", prefix, status, prefix)
  jma_path <- paste0(prefix, ".jma")
  if (!file.exists(jma_path)) return(list(selected = data.frame(), ld = matrix(numeric(), 0, 0)))
  jma <- utils::read.table(jma_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  snp_col <- .cater_pick_col(jma, c("SNP"), label = "COJO .jma SNP")
  selected <- data.frame(snp = as.character(jma[[snp_col]]), stringsAsFactors = FALSE)
  ldr_path <- paste0(prefix, ".jma.ldr")
  if (file.exists(ldr_path) && nrow(selected) > 1L) {
    ldr <- utils::read.table(ldr_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    rn <- as.character(ldr[[1L]])
    mat <- as.matrix(ldr[, -1L, drop = FALSE])
    storage.mode(mat) <- "double"
    rownames(mat) <- rn
    colnames(mat) <- names(ldr)[-1L]
    ld <- mat
  } else if (nrow(selected) == 1L) {
    ld <- matrix(1, 1, 1, dimnames = list(selected$snp, selected$snp))
  } else ld <- matrix(numeric(), 0, 0)
  list(selected = selected, ld = ld)
}

.cater_is_palindromic <- function(a1, a2) {
  paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
}

.cater_harmonize <- function(exp, outcome, drop_palindromic = TRUE) {
  y <- outcome[match(exp$snp, outcome$snp), , drop = FALSE]
  keep <- !is.na(y$snp)
  exp <- exp[keep, , drop = FALSE]
  y <- y[keep, , drop = FALSE]
  if (!nrow(exp)) return(data.frame())
  same <- exp$a1 == y$a1 & exp$a2 == y$a2
  swap <- exp$a1 == y$a2 & exp$a2 == y$a1
  keep <- same | swap
  if (drop_palindromic) keep <- keep & !.cater_is_palindromic(exp$a1, exp$a2)
  exp <- exp[keep, , drop = FALSE]
  y <- y[keep, , drop = FALSE]
  swap <- swap[keep]
  if (!nrow(exp)) return(data.frame())
  data.frame(
    snp = exp$snp,
    source = exp$source,
    parent_tf = exp$parent_tf,
    bx = exp$beta,
    bx_se = exp$se,
    by = ifelse(swap, -y$beta, y$beta),
    by_se = y$se,
    stringsAsFactors = FALSE
  )
}

.cater_subset_ld <- function(ld, snps) {
  if (!length(snps)) return(matrix(numeric(), 0, 0))
  if (length(snps) == 1L) return(matrix(1, 1, 1, dimnames = list(snps, snps)))
  if (!all(snps %in% rownames(ld)) || !all(snps %in% colnames(ld))) {
    .cater_stop("COJO LD matrix does not contain all harmonized selected SNPs")
  }
  ld[snps, snps, drop = FALSE]
}

.cater_givw <- function(dat, ld) {
  n <- nrow(dat)
  if (!n) return(data.frame(n_iv = 0L, beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_, mean_F = NA_real_, min_F = NA_real_))
  bx <- dat$bx; by <- dat$by; sy <- dat$by_se
  if (n == 1L) {
    b <- by / bx
    s <- abs(sy / bx)
    return(data.frame(n_iv = 1L, beta = b, se = s, p = 2 * stats::pnorm(-abs(b / s)),
                      Q = NA_real_, Q_p = NA_real_, mean_F = (bx / dat$bx_se)^2, min_F = (bx / dat$bx_se)^2))
  }
  D <- diag(sy, nrow = n)
  omega <- D %*% ld %*% D
  inv <- tryCatch(solve(omega), error = function(e) NULL)
  if (is.null(inv)) {
    eps <- max(diag(omega), na.rm = TRUE) * 1e-10
    inv <- solve(omega + diag(eps, n))
  }
  den <- as.numeric(crossprod(bx, inv %*% bx))
  if (!is.finite(den) || den <= 0) .cater_stop("Non-positive GIVW denominator")
  b <- as.numeric(crossprod(bx, inv %*% by)) / den
  s <- sqrt(1 / den)
  resid <- by - b * bx
  q <- as.numeric(crossprod(resid, inv %*% resid))
  f <- (bx / dat$bx_se)^2
  data.frame(n_iv = n, beta = b, se = s, p = 2 * stats::pnorm(-abs(b / s)),
             Q = q, Q_p = stats::pchisq(q, df = n - 1L, lower.tail = FALSE),
             mean_F = mean(f, na.rm = TRUE), min_F = min(f, na.rm = TRUE))
}

#' Run minimal CATER-MR analysis
#'
#' @param grn Cell-type-specific GRN data.frame, or list with $grn. Must contain TF and Target.
#' @param eqtl_dir Directory containing one genome-wide eQTL summary file per gene, named SYMBOL.txt.gz.
#' @param outcome Outcome GWAS summary-statistics data.frame.
#' @param gene_annotation Optional data.frame with symbol, chr and tss. Not needed if GRN embeds TF_chr/TF_tss and Target_chr/Target_tss.
#' @param ld_bfile PLINK bed/bim/fam prefix from an ancestry-matched LD reference; ideally the QTL donor genotypes.
#' @param gcta_bin GCTA executable name/path, default gcta64.
#' @param targets Optional character vector. Default is all unique TF and Target symbols in the GRN.
#' @param cis_window Cis window in bp around target TSS.
#' @param tf_window One-hop TF-locus window in bp around TF TSS. Defaults to cis_window.
#' @param cojo_p COJO selection threshold. V0.1 uses one threshold for cis and trans.
#' @param cojo_wind_kb GCTA --cojo-wind.
#' @param cojo_collinear GCTA --cojo-collinear.
#' @param qtl_n Optional constant QTL donor N when summary files do not contain N.
#' @param outdir Output directory.
#' @param drop_palindromic Drop A/T and C/G instruments during harmonization.
#' @param verbose Print progress.
#'
#' @return data.frame with cis, trans and combined MR results per target.
cater_mr <- function(grn,
                     eqtl_dir,
                     outcome,
                     gene_annotation = NULL,
                     ld_bfile,
                     gcta_bin = "gcta64",
                     targets = NULL,
                     cis_window = 1e6,
                     tf_window = cis_window,
                     cojo_p = 5e-8,
                     cojo_wind_kb = 10000L,
                     cojo_collinear = 0.9,
                     qtl_n = NULL,
                     outdir = "CATER_MR_results",
                     drop_palindromic = TRUE,
                     verbose = TRUE) {
  if (!dir.exists(eqtl_dir)) .cater_stop("eqtl_dir does not exist: %s", eqtl_dir)
  grn0 <- grn
  grn <- .cater_standardize_grn(grn)
  annotation <- .cater_standardize_annotation(gene_annotation)
  if (is.null(annotation)) annotation <- .cater_annotation_from_grn(grn0)
  if (is.null(annotation)) {
    .cater_stop("Gene coordinates are required. Supply gene_annotation with symbol/chr/tss, or embed TF_chr/TF_tss and Target_chr/Target_tss in the GRN object.")
  }
  outcome <- .cater_standardize_sumstats(outcome, label = "outcome", require_position = FALSE)
  if (is.null(targets)) targets <- sort(unique(c(grn$TF, grn$Target)))
  targets <- unique(as.character(targets))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "cojo"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "instruments"), recursive = TRUE, showWarnings = FALSE)

  all_results <- list()
  for (target in targets) {
    .cater_msg(verbose, "[%s] starting", target)
    qtl_file <- file.path(eqtl_dir, paste0(target, ".txt.gz"))
    if (!file.exists(qtl_file)) {
      all_results[[target]] <- data.frame(target = target, model = NA_character_, n_iv = 0L,
                                          beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_,
                                          mean_F = NA_real_, min_F = NA_real_, status = "NO_EQTL_FILE")
      next
    }
    parents <- unique(grn$TF[grn$Target == target])
    regions <- .cater_make_regions(target, parents, annotation, cis_window, tf_window)
    if (is.null(regions)) {
      all_results[[target]] <- data.frame(target = target, model = NA_character_, n_iv = 0L,
                                          beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_,
                                          mean_F = NA_real_, min_F = NA_real_, status = "NO_TARGET_ANNOTATION")
      next
    }
    qtl <- .cater_standardize_sumstats(.cater_read_table(qtl_file), n_default = qtl_n, label = paste0(target, " eQTL"))
    cmap <- .cater_candidate_map(qtl, regions, target)
    if (!nrow(cmap)) {
      all_results[[target]] <- data.frame(target = target, model = NA_character_, n_iv = 0L,
                                          beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_,
                                          mean_F = NA_real_, min_F = NA_real_, status = "NO_CANDIDATE_SNP")
      next
    }
    prefix <- file.path(outdir, "cojo", target)
    cojo <- tryCatch(
      .cater_run_cojo(qtl, cmap, ld_bfile, gcta_bin, cojo_p, cojo_wind_kb, cojo_collinear, prefix, verbose),
      error = function(e) e
    )
    if (inherits(cojo, "error")) {
      warning(sprintf("[%s] %s", target, conditionMessage(cojo)))
      all_results[[target]] <- data.frame(target = target, model = NA_character_, n_iv = 0L,
                                          beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_,
                                          mean_F = NA_real_, min_F = NA_real_, status = "COJO_FAILED")
      next
    }
    if (!nrow(cojo$selected)) {
      all_results[[target]] <- data.frame(target = target, model = NA_character_, n_iv = 0L,
                                          beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_,
                                          mean_F = NA_real_, min_F = NA_real_, status = "NO_COJO_SIGNAL")
      next
    }
    sel <- qtl[match(cojo$selected$snp, qtl$snp), , drop = FALSE]
    meta <- cmap[match(sel$snp, cmap$snp), , drop = FALSE]
    sel$source <- meta$source
    sel$parent_tf <- meta$parent_tf
    h <- .cater_harmonize(sel, outcome, drop_palindromic = drop_palindromic)
    if (!nrow(h)) {
      all_results[[target]] <- data.frame(target = target, model = NA_character_, n_iv = 0L,
                                          beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Q_p = NA_real_,
                                          mean_F = NA_real_, min_F = NA_real_, status = "NO_HARMONIZED_IV")
      next
    }
    utils::write.table(h, file.path(outdir, "instruments", paste0(target, ".tsv")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    models <- c("cis", "trans", "combined")
    rows <- lapply(models, function(m) {
      d <- if (m == "combined") h else h[h$source == m, , drop = FALSE]
      ld <- if (nrow(d)) .cater_subset_ld(cojo$ld, d$snp) else matrix(numeric(), 0, 0)
      est <- .cater_givw(d, ld)
      data.frame(target = target, model = m, est, status = if (nrow(d)) "OK" else paste0("NO_", toupper(m), "_IV"),
                 stringsAsFactors = FALSE)
    })
    all_results[[target]] <- do.call(rbind, rows)
  }
  result <- do.call(rbind, all_results)
  rownames(result) <- NULL
  utils::write.table(result, file.path(outdir, "cater_mr_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  result
}
