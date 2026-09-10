# CATER-MR scMORE GRN adapter
#
# CATER-MR obtains outcome-independent, cell-type-specific GRNs by subsetting a
# processed multiome Seurat object by cell type and then calling the upstream
# scMORE::createRegulon() function unchanged for each subset.
#
# Important distinction: upstream scMORE::createRegulon() itself constructs a
# global GRN for the cells supplied to it. The per-cell-type refit is a CATER-MR
# adapter around that exact engine; it is not claimed to be scMORE's native
# top-level scMore() workflow.

.CATER_SCMORE_REPO <- "mayunlong89/scMORE"
.CATER_SCMORE_AUDITED_SHA <- "f614736b9f49631471b0f18dfa415ee35e4d7b66"
.CATER_SCMORE_AUDITED_VERSION <- "2.0.0"

.cater_scmore_stop <- function(...) stop(sprintf(...), call. = FALSE)
.cater_scmore_msg <- function(verbose, ...) if (isTRUE(verbose)) message(sprintf(...))

.cater_scmore_provenance <- function(strict_upstream = TRUE) {
  if (!requireNamespace("scMORE", quietly = TRUE)) {
    .cater_scmore_stop(
      paste0(
        "scMORE is required. For the audited implementation install: ",
        "remotes::install_github('mayunlong89/scMORE@",
        .CATER_SCMORE_AUDITED_SHA, "')"
      )
    )
  }

  d <- utils::packageDescription("scMORE")
  remote_user <- d[["RemoteUsername"]]
  remote_repo <- d[["RemoteRepo"]]
  remote_sha <- d[["RemoteSha"]]
  version <- as.character(utils::packageVersion("scMORE"))

  verified <- !is.null(remote_user) && !is.null(remote_repo) && !is.null(remote_sha) &&
    identical(remote_user, "mayunlong89") && identical(remote_repo, "scMORE") &&
    identical(tolower(remote_sha), tolower(.CATER_SCMORE_AUDITED_SHA))

  if (isTRUE(strict_upstream) && !verified) {
    observed <- if (is.null(remote_sha)) "unavailable" else remote_sha
    .cater_scmore_stop(
      paste0(
        "strict_upstream=TRUE requires the audited scMORE GitHub commit ",
        .CATER_SCMORE_AUDITED_SHA, ". Installed scMORE version=", version,
        ", RemoteSha=", observed, ". Reinstall with remotes::install_github('",
        .CATER_SCMORE_REPO, "@", .CATER_SCMORE_AUDITED_SHA, "') or explicitly set ",
        "strict_upstream=FALSE after independently validating another upstream revision."
      )
    )
  }

  pkgver <- function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }

  list(
    scMORE_repository = .CATER_SCMORE_REPO,
    scMORE_audited_sha = .CATER_SCMORE_AUDITED_SHA,
    scMORE_remote_sha = if (is.null(remote_sha)) NA_character_ else remote_sha,
    scMORE_version = version,
    scMORE_sha_verified = verified,
    Pando_version = pkgver("Pando"),
    Seurat_version = pkgver("Seurat"),
    Signac_version = pkgver("Signac")
  )
}

.cater_scmore_validate_seurat <- function(single_cell) {
  if (!inherits(single_cell, "Seurat")) {
    .cater_scmore_stop("single_cell must be a processed Seurat object")
  }
  assays <- names(single_cell@assays)
  missing_assays <- setdiff(c("RNA", "peaks"), assays)
  if (length(missing_assays)) {
    .cater_scmore_stop(
      "scMORE::createRegulon() hard-codes assays 'RNA' and 'peaks'; missing: %s",
      paste(missing_assays, collapse = ", ")
    )
  }
  if (!requireNamespace("Signac", quietly = TRUE)) {
    .cater_scmore_stop("Signac is required by scMORE")
  }
  anno <- tryCatch(Signac::Annotation(single_cell[["peaks"]]), error = function(e) NULL)
  if (is.null(anno) || length(anno) == 0L) {
    .cater_scmore_stop(
      paste0(
        "The 'peaks' ChromatinAssay has no gene annotation. scMORE's documented ",
        "multiome workflow requires Annotation(single_cell[['peaks']]) to be populated ",
        "before createRegulon()."
      )
    )
  }
  invisible(TRUE)
}

.cater_scmore_cell_labels <- function(single_cell, celltype_col = NULL) {
  cells <- colnames(single_cell)
  if (!length(cells)) .cater_scmore_stop("single_cell has no cells")

  if (is.null(celltype_col)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      .cater_scmore_stop("SeuratObject is required to read Idents(single_cell)")
    }
    labels <- as.character(SeuratObject::Idents(single_cell))
    names(labels) <- cells
    source <- "Idents(single_cell)"
  } else {
    md <- single_cell[[]]
    if (!celltype_col %in% names(md)) {
      .cater_scmore_stop("celltype_col '%s' is not present in Seurat metadata", celltype_col)
    }
    labels <- as.character(md[[celltype_col]])
    names(labels) <- rownames(md)
    labels <- labels[match(cells, names(labels))]
    names(labels) <- cells
    source <- paste0("metadata$", celltype_col)
  }

  bad <- is.na(labels) | !nzchar(labels)
  if (any(bad)) {
    .cater_scmore_stop(
      "Cell-type labels contain %d missing/empty values; assign every cell before GRN inference",
      sum(bad)
    )
  }
  attr(labels, "source") <- source
  labels
}

.cater_scmore_create_args <- function(n_targets = 5,
                                      peak2gene_method = "Signac",
                                      infer_method = "glm",
                                      tss_upstream = 100000,
                                      tss_downstream = 0,
                                      exclude_exon_regions = TRUE,
                                      conserved_regions = NULL) {
  # These are the exact public createRegulon() arguments/defaults audited from
  # mayunlong89/scMORE. If conserved_regions is NULL it is intentionally omitted
  # so scMORE itself evaluates its own phastConsElements20Mammals.UCSC.hg38 default.
  args <- list(
    n_targets = n_targets,
    peak2gene_method = peak2gene_method,
    infer_method = infer_method,
    tss_upstream = tss_upstream,
    tss_downstream = tss_downstream,
    exclude_exon_regions = exclude_exon_regions
  )
  if (!is.null(conserved_regions)) args$conserved_regions <- conserved_regions
  args
}

.cater_scmore_call_create_regulon <- function(single_cell, args, create_fun = NULL) {
  if (is.null(create_fun)) create_fun <- getExportedValue("scMORE", "createRegulon")
  do.call(create_fun, c(list(single_cell = single_cell), args))
}

.cater_scmore_validate_output <- function(x, cell_type) {
  if (!is.list(x) || is.null(x$grn) || is.null(x$tf_names)) {
    .cater_scmore_stop(
      "scMORE::createRegulon() returned an unexpected object for cell type '%s'",
      cell_type
    )
  }
  if (!is.data.frame(x$grn) || !all(c("TF", "Target") %in% names(x$grn))) {
    .cater_scmore_stop(
      "scMORE::createRegulon() output for '%s' lacks data.frame grn with TF/Target columns",
      cell_type
    )
  }
  invisible(TRUE)
}

.cater_scmore_safe_name <- function(x) {
  y <- gsub("[^A-Za-z0-9._-]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  ifelse(nzchar(y), y, "celltype")
}

#' Build cell-type-specific GRNs with the exact scMORE createRegulon engine
#'
#' The processed multiome Seurat object is split by a cell-type annotation, and
#' each subset is passed unchanged to scMORE::createRegulon(). CATER-MR does not
#' reproduce or modify the internal scMORE/Pando GRN mathematics.
#'
#' Upstream scMORE createRegulon() performs, in order:
#' FindVariableFeatures(RNA), Pando::initiate_grn(), Pando::find_motifs(),
#' Pando::infer_grn(), Pando::find_modules(), NetworkModules(), extract_grn(),
#' then filters TF regulons by n_targets.
#'
#' @param single_cell Processed Seurat multiome object. Assays must be named
#'   exactly "RNA" and "peaks". The peaks ChromatinAssay must already contain
#'   gene annotation. scMORE assumes GRCh38/hg38 resources.
#' @param celltype_col Metadata column containing cell-type labels. If NULL,
#'   Seurat Idents(single_cell) are used.
#' @param cell_types Optional character vector selecting cell types to fit.
#' @param n_targets,peak2gene_method,infer_method,tss_upstream,tss_downstream,
#'   exclude_exon_regions,conserved_regions Arguments forwarded unchanged to
#'   scMORE::createRegulon(). Defaults mirror audited upstream code. When
#'   conserved_regions=NULL, the argument is omitted so scMORE uses its own
#'   phastConsElements20Mammals.UCSC.hg38 default.
#' @param strict_upstream If TRUE (default), require the installed scMORE package
#'   to report the audited GitHub RemoteSha.
#' @param on_error "stop" (default) or "record". The latter records a failed
#'   cell type without silently substituting a different GRN method.
#' @param outdir Optional directory for raw scMORE RDS outputs and edge tables.
#' @param gc_after_each Run gc() after each cell-type fit to limit peak memory.
#' @param verbose Print cell-type progress messages.
#'
#' @return A cater_scmore_grn object containing raw upstream outputs, per-cell
#'   type GRN tables, fit summary and software provenance.
cater_build_scmore_grn <- function(single_cell,
                                    celltype_col = NULL,
                                    cell_types = NULL,
                                    n_targets = 5,
                                    peak2gene_method = "Signac",
                                    infer_method = "glm",
                                    tss_upstream = 100000,
                                    tss_downstream = 0,
                                    exclude_exon_regions = TRUE,
                                    conserved_regions = NULL,
                                    strict_upstream = TRUE,
                                    on_error = c("stop", "record"),
                                    outdir = NULL,
                                    gc_after_each = TRUE,
                                    verbose = TRUE) {
  on_error <- match.arg(on_error)
  provenance <- .cater_scmore_provenance(strict_upstream = strict_upstream)
  .cater_scmore_validate_seurat(single_cell)
  labels <- .cater_scmore_cell_labels(single_cell, celltype_col = celltype_col)
  label_source <- attr(labels, "source")

  available <- unique(unname(labels))
  if (is.null(cell_types)) {
    cell_types <- available
  } else {
    cell_types <- unique(as.character(cell_types))
    unknown <- setdiff(cell_types, available)
    if (length(unknown)) {
      .cater_scmore_stop("Requested cell_types are absent: %s", paste(unknown, collapse = ", "))
    }
  }
  if (!length(cell_types)) .cater_scmore_stop("No cell types selected")

  create_args <- .cater_scmore_create_args(
    n_targets = n_targets,
    peak2gene_method = peak2gene_method,
    infer_method = infer_method,
    tss_upstream = tss_upstream,
    tss_downstream = tss_downstream,
    exclude_exon_regions = exclude_exon_regions,
    conserved_regions = conserved_regions
  )

  if (!is.null(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  raw_outputs <- setNames(vector("list", length(cell_types)), cell_types)
  grns <- setNames(vector("list", length(cell_types)), cell_types)
  summary_rows <- vector("list", length(cell_types))

  for (i in seq_along(cell_types)) {
    ct <- cell_types[[i]]
    cells <- names(labels)[labels == ct]
    .cater_scmore_msg(verbose, "[scMORE GRN] %s: %d cells", ct, length(cells))

    # This is the only CATER-MR algorithmic adaptation: define the cell population.
    # All subsequent GRN inference is delegated to the unmodified upstream function.
    ct_object <- base::subset(single_cell, cells = cells)
    fit <- tryCatch(
      .cater_scmore_call_create_regulon(ct_object, create_args),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      if (on_error == "stop") {
        .cater_scmore_stop("scMORE GRN failed for cell type '%s': %s", ct, conditionMessage(fit))
      }
      summary_rows[[i]] <- data.frame(
        cell_type = ct, n_cells = length(cells), n_edges = NA_integer_,
        n_tfs = NA_integer_, n_targets = NA_integer_, status = "SCMORE_FAILED",
        error_message = conditionMessage(fit), stringsAsFactors = FALSE
      )
      rm(ct_object, fit)
      if (isTRUE(gc_after_each)) gc(verbose = FALSE)
      next
    }

    .cater_scmore_validate_output(fit, ct)
    raw_outputs[[ct]] <- fit
    g <- fit$grn
    g$cell_type <- ct
    grns[[ct]] <- g

    summary_rows[[i]] <- data.frame(
      cell_type = ct,
      n_cells = length(cells),
      n_edges = nrow(g),
      n_tfs = length(unique(g$TF)),
      n_targets = length(unique(g$Target)),
      status = "OK",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )

    if (!is.null(outdir)) {
      safe <- .cater_scmore_safe_name(ct)
      saveRDS(fit, file.path(outdir, paste0("scMORE_raw_", safe, ".rds")))
      utils::write.table(
        g,
        file.path(outdir, paste0("scMORE_grn_", safe, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }

    rm(ct_object, fit, g)
    if (isTRUE(gc_after_each)) gc(verbose = FALSE)
  }

  fit_summary <- do.call(rbind, summary_rows)
  rownames(fit_summary) <- NULL

  if (!is.null(outdir)) {
    utils::write.table(
      fit_summary,
      file.path(outdir, "scMORE_grn_summary.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    saveRDS(provenance, file.path(outdir, "scMORE_provenance.rds"))
  }

  structure(
    list(
      grns = grns,
      scmore_outputs = raw_outputs,
      summary = fit_summary,
      provenance = provenance,
      celltype_source = label_source,
      createRegulon_args = create_args,
      mode = "per_celltype_subset_then_unmodified_scMORE_createRegulon"
    ),
    class = "cater_scmore_grn"
  )
}

#' Extract one cell-type GRN from cater_build_scmore_grn()
cater_get_scmore_grn <- function(x, cell_type = NULL, raw = FALSE) {
  if (!inherits(x, "cater_scmore_grn")) {
    .cater_scmore_stop("x must be a cater_scmore_grn object")
  }
  available <- names(x$grns)[!vapply(x$grns, is.null, logical(1))]
  if (is.null(cell_type)) {
    if (length(available) != 1L) {
      .cater_scmore_stop("Specify cell_type; available successful fits: %s", paste(available, collapse = ", "))
    }
    cell_type <- available[[1L]]
  }
  if (!cell_type %in% names(x$grns) || is.null(x$grns[[cell_type]])) {
    .cater_scmore_stop("No successful scMORE GRN for cell type '%s'", cell_type)
  }
  if (isTRUE(raw)) x$scmore_outputs[[cell_type]] else x$grns[[cell_type]]
}

print.cater_scmore_grn <- function(x, ...) {
  ok <- sum(x$summary$status == "OK", na.rm = TRUE)
  cat("CATER-MR scMORE cell-type GRNs\n")
  cat("  mode:", x$mode, "\n")
  cat("  successful:", ok, "/", nrow(x$summary), "cell types\n")
  cat("  scMORE version:", x$provenance$scMORE_version, "\n")
  cat("  scMORE RemoteSha:", x$provenance$scMORE_remote_sha, "\n")
  invisible(x)
}
