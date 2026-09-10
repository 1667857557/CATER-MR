# CATER-MR combined eQTL input adapter
#
# This adapter lets the unchanged core estimator consume a single merged
# cell-type-specific full-summary eQTL table with one gene-identifying column
# (default candidates include GENE/gene/symbol/gene_symbol).
#
# Mathematical contract:
#   merged_table[GENE == X, ] == Q_X
# The adapter only changes storage/layout. It does not filter SNPs, alter
# betas/SEs/P values, or change Manc-COJO / MR calculations.

.cater_eqtl_input_stop <- function(...) stop(sprintf(...), call. = FALSE)
.cater_eqtl_input_msg <- function(verbose, ...) if (isTRUE(verbose)) message(sprintf(...))

.cater_require_core <- function() {
  need <- c("cater_mr", ".cater_read_table", ".cater_standardize_sumstats")
  miss <- need[!vapply(need, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(miss)) {
    .cater_eqtl_input_stop(
      "Source CATER_MR.R before CATER_EQTL_INPUT.R; missing core functions: %s",
      paste(miss, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.cater_eqtl_pick_gene_col <- function(x, gene_col = NULL) {
  if (!is.data.frame(x)) .cater_eqtl_input_stop("Combined eQTL input must be a data.frame after reading")
  nm <- names(x)
  if (!is.null(gene_col)) {
    if (length(gene_col) != 1L || !nzchar(gene_col) || !gene_col %in% nm) {
      .cater_eqtl_input_stop("gene_col must name one existing column; requested: %s", paste(gene_col, collapse = ","))
    }
    return(gene_col)
  }
  candidates <- c("GENE", "gene", "gene_name", "gene_symbol", "symbol", "SYMBOL")
  hit <- match(tolower(candidates), tolower(nm), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    .cater_eqtl_input_stop(
      "Combined eQTL table needs a gene column. Tried: %s. Supply gene_col explicitly for another name.",
      paste(candidates, collapse = ", ")
    )
  }
  nm[hit[1L]]
}

.cater_read_combined_eqtl <- function(eqtl_table) {
  if (is.data.frame(eqtl_table)) return(eqtl_table)
  if (!is.character(eqtl_table) || length(eqtl_table) != 1L || !nzchar(eqtl_table)) {
    .cater_eqtl_input_stop("eqtl_table must be a data.frame or one tab-delimited .txt/.txt.gz file path")
  }
  if (!file.exists(eqtl_table)) .cater_eqtl_input_stop("Combined eQTL file does not exist: %s", eqtl_table)
  .cater_require_core()
  .cater_read_table(eqtl_table)
}

.cater_validate_gene_names <- function(gene) {
  gene <- trimws(as.character(gene))
  bad <- is.na(gene) | !nzchar(gene)
  if (any(bad)) {
    .cater_eqtl_input_stop("Combined eQTL table has %d rows with missing/empty gene names", sum(bad))
  }
  # The core directory contract is <SYMBOL>.txt.gz; path separators would make
  # the same biological identifier map to a different filesystem location.
  bad_path <- grepl("[/\\\\]", gene)
  if (any(bad_path)) {
    ex <- unique(gene[bad_path])
    .cater_eqtl_input_stop(
      "Gene names cannot contain path separators under the per-gene cache contract (e.g. %s)",
      paste(utils::head(ex, 5L), collapse = ", ")
    )
  }
  gene
}

#' Materialize one merged eQTL table into the core per-gene storage contract
#'
#' @param eqtl_table data.frame or tab-delimited .txt/.txt.gz path. One column
#'   identifies the exposure gene and all other columns are ordinary CATER-MR
#'   full-summary eQTL columns.
#' @param gene_col Gene column name. If NULL, common names including GENE,
#'   gene, symbol and gene_symbol are detected case-insensitively.
#' @param outdir Destination directory. If NULL, a temporary directory is made.
#' @param overwrite Whether existing <GENE>.txt.gz cache files may be replaced.
#' @param validate_sumstats If TRUE, validate every gene subset using the exact
#'   core summary-stat parser before writing. This does not rewrite effect data.
#' @param qtl_n Optional donor N used only during validation when N is absent.
#' @param verbose Print progress.
#'
#' @return cater_eqtl_cache object containing eqtl_dir, gene/row summary and
#'   whether the cache is temporary.
cater_prepare_eqtl_table <- function(eqtl_table,
                                     gene_col = NULL,
                                     outdir = NULL,
                                     overwrite = FALSE,
                                     validate_sumstats = TRUE,
                                     qtl_n = NULL,
                                     verbose = TRUE) {
  .cater_require_core()
  x <- .cater_read_combined_eqtl(eqtl_table)
  gcol <- .cater_eqtl_pick_gene_col(x, gene_col)
  gene <- .cater_validate_gene_names(x[[gcol]])
  if (!nrow(x)) .cater_eqtl_input_stop("Combined eQTL table has zero rows")

  temporary <- is.null(outdir)
  if (temporary) outdir <- tempfile("cater_eqtl_gene_cache_")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(outdir)) .cater_eqtl_input_stop("Could not create eQTL cache directory: %s", outdir)

  # Sort once so each gene is written from one contiguous block. This avoids
  # split(data.frame, gene), which can duplicate a very large table in memory.
  ord <- order(gene, method = "radix")
  gene_ord <- gene[ord]
  rr <- rle(gene_ord)
  ends <- cumsum(rr$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  payload_cols <- setdiff(names(x), gcol)
  if (!length(payload_cols)) .cater_eqtl_input_stop("Combined eQTL table contains only the gene column")

  summary_rows <- vector("list", length(rr$values))
  for (i in seq_along(rr$values)) {
    g <- rr$values[[i]]
    idx <- ord[starts[[i]]:ends[[i]]]
    d <- x[idx, payload_cols, drop = FALSE]
    if (isTRUE(validate_sumstats)) {
      # Exact same parser used later by cater_mr(); validation catches malformed
      # per-gene groups early but the original columns/values are persisted.
      .cater_standardize_sumstats(d, n_default = qtl_n, label = paste0(g, " eQTL"))
    }
    f <- file.path(outdir, paste0(g, ".txt.gz"))
    if (file.exists(f) && !isTRUE(overwrite)) {
      .cater_eqtl_input_stop("Cache file already exists for gene '%s': %s (set overwrite=TRUE to replace)", g, f)
    }
    con <- gzfile(f, "wt")
    tryCatch(
      utils::write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE),
      finally = close(con)
    )
    summary_rows[[i]] <- data.frame(gene = g, n_rows = nrow(d), file = f, stringsAsFactors = FALSE)
    if (i %% 100L == 0L || i == length(rr$values)) {
      .cater_eqtl_input_msg(verbose, "Materialized %d/%d genes", i, length(rr$values))
    }
  }

  sm <- do.call(rbind, summary_rows)
  rownames(sm) <- NULL
  structure(
    list(
      eqtl_dir = normalizePath(outdir, winslash = "/", mustWork = TRUE),
      gene_col = gcol,
      n_rows = nrow(x),
      n_genes = nrow(sm),
      genes = sm,
      temporary = temporary,
      source = if (is.data.frame(eqtl_table)) "data.frame" else normalizePath(eqtl_table, winslash = "/", mustWork = TRUE)
    ),
    class = "cater_eqtl_cache"
  )
}

#' Run the unchanged CATER-MR estimator from one merged eQTL table
#'
#' This is a storage adapter around cater_mr(). It materializes per-gene files
#' and then calls the exact same core estimator, so Manc-COJO, GIVW, TF-anchor,
#' sibling screening and triggered local MVMR are identical to directory mode.
cater_mr_from_table <- function(grn,
                                eqtl_table,
                                outcome,
                                ...,
                                gene_col = NULL,
                                eqtl_cache_dir = NULL,
                                eqtl_cache_overwrite = FALSE,
                                keep_eqtl_cache = FALSE,
                                validate_sumstats = TRUE,
                                materialize_verbose = TRUE) {
  .cater_require_core()
  dots <- list(...)
  qtl_n <- if ("qtl_n" %in% names(dots)) dots$qtl_n else NULL
  prep <- cater_prepare_eqtl_table(
    eqtl_table = eqtl_table,
    gene_col = gene_col,
    outdir = eqtl_cache_dir,
    overwrite = eqtl_cache_overwrite,
    validate_sumstats = validate_sumstats,
    qtl_n = qtl_n,
    verbose = materialize_verbose
  )
  if (isTRUE(prep$temporary) && !isTRUE(keep_eqtl_cache)) {
    on.exit(unlink(prep$eqtl_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  do.call(
    cater_mr,
    c(list(grn = grn, eqtl_dir = prep$eqtl_dir, outcome = outcome), dots)
  )
}

print.cater_eqtl_cache <- function(x, ...) {
  cat("CATER-MR combined eQTL cache\n")
  cat("  source:", x$source, "\n")
  cat("  gene column:", x$gene_col, "\n")
  cat("  rows:", x$n_rows, "\n")
  cat("  genes:", x$n_genes, "\n")
  cat("  eqtl_dir:", x$eqtl_dir, "\n")
  invisible(x)
}
