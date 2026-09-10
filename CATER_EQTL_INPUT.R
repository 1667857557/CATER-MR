# CATER-MR combined eQTL input adapter v0.6
#
# Storage contract only: merged_table[GENE == X, ] == Q_X.
# The adapter never performs significance filtering and never changes MR/COJO
# mathematics. It can materialize only GRN-relevant genes to avoid needless disk
# amplification, while direct calls still default to all genes for compatibility.

.cater_eqtl_input_stop <- function(...) stop(sprintf(...), call. = FALSE)
.cater_eqtl_input_msg <- function(verbose, ...) if (isTRUE(verbose)) message(sprintf(...))

.cater_require_core <- function() {
  need <- c("cater_mr", ".cater_read_table", ".cater_standardize_sumstats", ".cater_standardize_grn")
  miss <- need[!vapply(need, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(miss)) .cater_eqtl_input_stop("Source CATER_MR.R before CATER_EQTL_INPUT.R; missing: %s", paste(miss, collapse=", "))
  invisible(TRUE)
}

.cater_eqtl_pick_gene_col <- function(x, gene_col = NULL) {
  if (!is.data.frame(x)) .cater_eqtl_input_stop("Combined eQTL input must be a data.frame after reading")
  nm <- names(x)
  if (!is.null(gene_col)) {
    if (length(gene_col) != 1L || !nzchar(gene_col) || !gene_col %in% nm)
      .cater_eqtl_input_stop("gene_col must name one existing column; requested: %s", paste(gene_col, collapse=","))
    return(gene_col)
  }
  candidates <- c("GENE", "gene", "gene_name", "gene_symbol", "symbol", "SYMBOL")
  hit <- match(tolower(candidates), tolower(nm), nomatch = 0L); hit <- hit[hit > 0L]
  if (!length(hit)) .cater_eqtl_input_stop("Combined eQTL table needs a gene column. Tried: %s", paste(candidates, collapse=", "))
  nm[hit[1L]]
}

.cater_read_combined_eqtl <- function(eqtl_table) {
  if (is.data.frame(eqtl_table)) return(eqtl_table)
  if (!is.character(eqtl_table) || length(eqtl_table) != 1L || !nzchar(eqtl_table))
    .cater_eqtl_input_stop("eqtl_table must be a data.frame or one tab-delimited .txt/.txt.gz path")
  if (!file.exists(eqtl_table)) .cater_eqtl_input_stop("Combined eQTL file does not exist: %s", eqtl_table)
  .cater_require_core(); .cater_read_table(eqtl_table)
}

.cater_validate_gene_names <- function(gene) {
  gene <- trimws(as.character(gene)); bad <- is.na(gene) | !nzchar(gene)
  if (any(bad)) .cater_eqtl_input_stop("Combined eQTL table has %d rows with missing/empty gene names", sum(bad))
  bad_path <- grepl("[/\\\\]", gene)
  if (any(bad_path)) .cater_eqtl_input_stop("Gene names cannot contain path separators (e.g. %s)", paste(utils::head(unique(gene[bad_path]),5L),collapse=", "))
  gene
}

.cater_needed_genes_from_grn <- function(grn, targets = NULL) {
  g <- .cater_standardize_grn(grn)
  needed <- unique(c(g$TF, g$Target, as.character(targets)))
  needed[!is.na(needed) & nzchar(needed)]
}

cater_prepare_eqtl_table <- function(eqtl_table,
                                     gene_col = NULL,
                                     genes = NULL,
                                     outdir = NULL,
                                     overwrite = FALSE,
                                     validate_sumstats = TRUE,
                                     qtl_n = NULL,
                                     verbose = TRUE) {
  .cater_require_core()
  x <- .cater_read_combined_eqtl(eqtl_table)
  if (!nrow(x)) .cater_eqtl_input_stop("Combined eQTL table has zero rows")
  gcol <- .cater_eqtl_pick_gene_col(x, gene_col)
  gene <- .cater_validate_gene_names(x[[gcol]])

  requested_genes <- NULL
  if (!is.null(genes)) {
    requested_genes <- unique(trimws(as.character(genes)))
    requested_genes <- requested_genes[!is.na(requested_genes) & nzchar(requested_genes)]
    keep <- gene %in% requested_genes
    x <- x[keep,,drop=FALSE]; gene <- gene[keep]
    if (!nrow(x)) .cater_eqtl_input_stop("None of the requested genes are present in the combined eQTL table")
  }

  temporary <- is.null(outdir)
  if (temporary) outdir <- tempfile("cater_eqtl_gene_cache_")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(outdir)) .cater_eqtl_input_stop("Could not create eQTL cache directory: %s", outdir)

  ord <- order(gene, method = "radix"); gene_ord <- gene[ord]; rr <- rle(gene_ord)
  ends <- cumsum(rr$lengths); starts <- c(1L, head(ends, -1L) + 1L)
  payload_cols <- setdiff(names(x), gcol)
  if (!length(payload_cols)) .cater_eqtl_input_stop("Combined eQTL table contains only the gene column")

  summary_rows <- vector("list", length(rr$values))
  for (i in seq_along(rr$values)) {
    g <- rr$values[[i]]; idx <- ord[starts[[i]]:ends[[i]]]; d <- x[idx, payload_cols, drop = FALSE]
    qc <- NULL
    if (isTRUE(validate_sumstats)) {
      std <- .cater_standardize_sumstats(d, n_default=qtl_n, label=paste0(g," eQTL"))
      qc <- attr(std,"qc")
      if (!nrow(std)) .cater_eqtl_input_stop("Gene '%s' has no valid full-summary rows after core validation", g)
    }
    f <- file.path(outdir, paste0(g, ".txt.gz"))
    if (file.exists(f) && !isTRUE(overwrite)) .cater_eqtl_input_stop("Cache file already exists for gene '%s': %s", g, f)
    con <- gzfile(f, "wt")
    tryCatch(utils::write.table(d, con, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE), finally=close(con))
    summary_rows[[i]] <- data.frame(gene=g,n_rows=nrow(d),n_valid_rows=if(is.null(qc)) NA_integer_ else qc$n_output,file=f,stringsAsFactors=FALSE)
    if (i %% 100L == 0L || i == length(rr$values)) .cater_eqtl_input_msg(verbose,"Materialized %d/%d genes",i,length(rr$values))
  }
  sm <- do.call(rbind, summary_rows); rownames(sm) <- NULL
  missing_requested <- if (is.null(requested_genes)) character() else setdiff(requested_genes, sm$gene)
  structure(list(eqtl_dir=normalizePath(outdir,winslash="/",mustWork=TRUE),gene_col=gcol,n_rows=nrow(x),n_genes=nrow(sm),genes=sm,
                 requested_genes=requested_genes,missing_requested_genes=missing_requested,temporary=temporary,
                 source=if(is.data.frame(eqtl_table))"data.frame" else normalizePath(eqtl_table,winslash="/",mustWork=TRUE)),class="cater_eqtl_cache")
}

cater_mr_from_table <- function(grn,
                                eqtl_table,
                                outcome,
                                ...,
                                gene_col = NULL,
                                genes = NULL,
                                eqtl_cache_dir = NULL,
                                eqtl_cache_overwrite = FALSE,
                                keep_eqtl_cache = FALSE,
                                validate_sumstats = TRUE,
                                materialize_verbose = TRUE) {
  .cater_require_core(); dots <- list(...)
  qtl_n <- if ("qtl_n" %in% names(dots)) dots$qtl_n else NULL
  targets <- if ("targets" %in% names(dots)) dots$targets else NULL
  if (is.null(genes)) genes <- .cater_needed_genes_from_grn(grn, targets)
  prep <- cater_prepare_eqtl_table(eqtl_table=eqtl_table,gene_col=gene_col,genes=genes,outdir=eqtl_cache_dir,
                                   overwrite=eqtl_cache_overwrite,validate_sumstats=validate_sumstats,qtl_n=qtl_n,verbose=materialize_verbose)
  if (isTRUE(prep$temporary) && !isTRUE(keep_eqtl_cache)) on.exit(unlink(prep$eqtl_dir, recursive=TRUE, force=TRUE), add=TRUE)
  do.call(cater_mr, c(list(grn=grn,eqtl_dir=prep$eqtl_dir,outcome=outcome),dots))
}

print.cater_eqtl_cache <- function(x, ...) {
  cat("CATER-MR combined eQTL cache\n")
  cat("  source:",x$source,"\n  gene column:",x$gene_col,"\n  rows materialized:",x$n_rows,"\n  genes:",x$n_genes,"\n  eqtl_dir:",x$eqtl_dir,"\n")
  if(length(x$missing_requested_genes)) cat("  requested genes absent:",length(x$missing_requested_genes),"\n")
  invisible(x)
}
