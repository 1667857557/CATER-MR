# CATER-MR scMORE GRN adapter
#
# CATER-MR obtains outcome-independent, cell-type-specific GRNs by subsetting a
# processed multiome Seurat object by cell type and then calling the upstream
# scMORE::createRegulon() function unchanged for each subset.
#
# Upstream scMORE::createRegulon() returns TF-target-region/evidence rows. CATER-MR
# itself requires a direct one-hop TF -> Target graph and gene coordinates for the
# target-centric cis + parent-TF-locus construction. This file therefore keeps the
# raw scMORE output unchanged and creates a deterministic CATER-ready graph layer.

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

.cater_scmore_check_genome_values <- function(genome_values, source = "peak annotation") {
  z <- unique(as.character(genome_values))
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) return("UNKNOWN")
  ok <- grepl("^(hg38|GRCh38)([._-].*)?$", z, ignore.case = TRUE)
  if (!all(ok)) {
    .cater_scmore_stop(
      "%s explicitly reports a non-hg38 genome assembly: %s. scMORE motif inference is hard-coded to BSgenome.Hsapiens.UCSC.hg38.",
      source, paste(z[!ok], collapse = ", ")
    )
  }
  paste(sort(unique(z)), collapse = ";")
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
  if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    g <- tryCatch(GenomeInfoDb::genome(anno), error = function(e) character())
    .cater_scmore_check_genome_values(g, "Signac peak gene annotation")
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

.cater_scmore_pick_col <- function(x, candidates, required = TRUE, label = "column") {
  nm <- names(x)
  hit <- match(tolower(candidates), tolower(nm), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(nm[hit[1L]])
  if (required) {
    .cater_scmore_stop("Missing %s. Tried: %s", label, paste(candidates, collapse = ", "))
  }
  NULL
}

.cater_scmore_standardize_gene_annotation <- function(gene_annotation) {
  if (!is.data.frame(gene_annotation)) .cater_scmore_stop("gene_annotation must be a data.frame")
  sym <- .cater_scmore_pick_col(gene_annotation, c("symbol","gene","gene_symbol","SYMBOL"), label="gene annotation symbol")
  chr <- .cater_scmore_pick_col(gene_annotation, c("chr","CHR","chrom","chromosome"), label="gene annotation chromosome")
  tss <- .cater_scmore_pick_col(gene_annotation, c("tss","TSS","txStart"), label="gene annotation TSS")
  genome_col <- .cater_scmore_pick_col(gene_annotation, c("genome","assembly","genome_build"), required=FALSE)
  if (!is.null(genome_col)) .cater_scmore_check_genome_values(gene_annotation[[genome_col]], "gene_annotation")

  out <- data.frame(
    symbol=as.character(gene_annotation[[sym]]),
    chr=sub("^chr","",as.character(gene_annotation[[chr]]),ignore.case=TRUE),
    tss=suppressWarnings(as.numeric(gene_annotation[[tss]])), stringsAsFactors=FALSE
  )
  out <- out[!is.na(out$symbol)&nzchar(out$symbol)&!is.na(out$chr)&nzchar(out$chr)&is.finite(out$tss),,drop=FALSE]
  key <- paste(out$chr, format(out$tss, scientific=FALSE, trim=TRUE), sep=":")
  amb <- vapply(split(key,out$symbol),function(z) length(unique(z))>1L,logical(1))
  if (any(amb)) {
    bad <- names(amb)[amb]
    .cater_scmore_stop(
      "gene_annotation has multiple distinct chr/TSS values for %d symbols (e.g. %s). Provide one gene-level TSS per symbol.",
      length(bad), paste(utils::head(bad,10L),collapse=", ")
    )
  }
  out <- out[!duplicated(out$symbol),,drop=FALSE]
  rownames(out) <- NULL
  out
}

.cater_scmore_annotation_from_object <- function(single_cell) {
  if (!requireNamespace("Signac",quietly=TRUE)) .cater_scmore_stop("Signac is required")
  anno <- Signac::Annotation(single_cell[["peaks"]])
  ad <- as.data.frame(anno)
  if (!nrow(ad)) .cater_scmore_stop("The peaks annotation is empty")
  sym <- .cater_scmore_pick_col(ad,c("gene_name","gene","symbol","gene_symbol","SYMBOL"),label="gene symbol in Signac Annotation(peaks)")
  seqcol <- .cater_scmore_pick_col(ad,c("seqnames","chr","chromosome"),label="chromosome in peak annotation")
  startcol <- .cater_scmore_pick_col(ad,c("start"),label="start in peak annotation")
  endcol <- .cater_scmore_pick_col(ad,c("end"),label="end in peak annotation")
  strandcol <- .cater_scmore_pick_col(ad,c("strand"),label="strand in peak annotation")
  strand <- as.character(ad[[strandcol]])
  keep <- strand %in% c("+","-") & !is.na(ad[[sym]]) & nzchar(as.character(ad[[sym]]))
  if (!any(keep)) .cater_scmore_stop("Cannot derive gene TSS from Annotation(peaks): no rows have both a gene symbol and '+'/'-' strand")
  tss <- ifelse(strand[keep]=="-", suppressWarnings(as.numeric(ad[[endcol]][keep])), suppressWarnings(as.numeric(ad[[startcol]][keep])))
  out <- data.frame(symbol=as.character(ad[[sym]][keep]), chr=sub("^chr","",as.character(ad[[seqcol]][keep]),ignore.case=TRUE), tss=tss, stringsAsFactors=FALSE)
  .cater_scmore_standardize_gene_annotation(out)
}

.cater_scmore_make_cater_grn <- function(raw_grn, gene_annotation, cell_type,
                                         drop_self_loops=TRUE,
                                         missing_coordinate=c("error","drop")) {
  missing_coordinate <- match.arg(missing_coordinate)
  if (!is.data.frame(raw_grn)||!all(c("TF","Target")%in%names(raw_grn))) .cater_scmore_stop("raw_grn must contain TF and Target")
  edges <- unique(data.frame(TF=as.character(raw_grn$TF),Target=as.character(raw_grn$Target),stringsAsFactors=FALSE))
  edges <- edges[!is.na(edges$TF)&nzchar(edges$TF)&!is.na(edges$Target)&nzchar(edges$Target),,drop=FALSE]
  n_raw_direct <- nrow(edges)
  n_self <- sum(edges$TF==edges$Target)
  if (isTRUE(drop_self_loops)&&n_self) edges <- edges[edges$TF!=edges$Target,,drop=FALSE]

  template <- data.frame(TF=character(),Target=character(),TF_chr=character(),TF_tss=numeric(),Target_chr=character(),Target_tss=numeric(),cell_type=character(),stringsAsFactors=FALSE)
  if (!nrow(edges)) {
    attr(template,"n_unique_raw_edges") <- n_raw_direct
    attr(template,"n_self_loops_dropped") <- if(isTRUE(drop_self_loops)) n_self else 0L
    attr(template,"n_missing_coordinate_edges") <- 0L
    return(template)
  }

  ann <- .cater_scmore_standardize_gene_annotation(gene_annotation)
  itf <- match(edges$TF,ann$symbol); itg <- match(edges$Target,ann$symbol)
  missing <- is.na(itf)|is.na(itg)
  if (any(missing)&&identical(missing_coordinate,"error")) {
    nodes <- unique(c(edges$TF[is.na(itf)],edges$Target[is.na(itg)]))
    .cater_scmore_stop("CATER-ready GRN for '%s' lacks chr/TSS coordinates for %d nodes (e.g. %s)",cell_type,length(nodes),paste(utils::head(nodes,10L),collapse=", "))
  }
  n_missing_edges <- sum(missing)
  if (any(missing)) { edges<-edges[!missing,,drop=FALSE]; itf<-itf[!missing]; itg<-itg[!missing] }
  out <- data.frame(TF=edges$TF,Target=edges$Target,TF_chr=ann$chr[itf],TF_tss=ann$tss[itf],Target_chr=ann$chr[itg],Target_tss=ann$tss[itg],cell_type=rep(cell_type,nrow(edges)),stringsAsFactors=FALSE)
  attr(out,"n_unique_raw_edges") <- n_raw_direct
  attr(out,"n_self_loops_dropped") <- if(isTRUE(drop_self_loops)) n_self else 0L
  attr(out,"n_missing_coordinate_edges") <- n_missing_edges
  out
}

.cater_scmore_safe_name <- function(x) {
  y <- gsub("[^A-Za-z0-9._-]+","_",x); y <- gsub("^_+|_+$","",y)
  ifelse(nzchar(y),y,"celltype")
}

.cater_scmore_output_names <- function(cell_types) {
  base <- vapply(cell_types,.cater_scmore_safe_name,character(1))
  rank_id <- match(cell_types,sort(unique(cell_types)))
  paste0(base,"__",sprintf("%03d",rank_id))
}

.cater_scmore_clean_message <- function(x) {
  if (!length(x)||is.na(x)) return(NA_character_)
  trimws(gsub("[\r\n\t]+"," ",as.character(x)))
}

#' Build CATER-ready cell-type GRNs with the exact scMORE createRegulon engine
#'
#' Raw scMORE TF-target-region/evidence rows are retained unchanged. CATER-MR
#' additionally constructs one unique direct TF -> Target edge per cell type and
#' appends gene-level hg38 chr/TSS coordinates needed for target cis and parent-TF loci.
cater_build_scmore_grn <- function(single_cell,
                                    celltype_col=NULL,
                                    cell_types=NULL,
                                    gene_annotation=NULL,
                                    n_targets=5,
                                    peak2gene_method="Signac",
                                    infer_method="glm",
                                    tss_upstream=100000,
                                    tss_downstream=0,
                                    exclude_exon_regions=TRUE,
                                    conserved_regions=NULL,
                                    drop_self_loops=TRUE,
                                    missing_coordinate=c("error","drop"),
                                    strict_upstream=TRUE,
                                    on_error=c("stop","record"),
                                    outdir=NULL,
                                    gc_after_each=TRUE,
                                    verbose=TRUE) {
  on_error <- match.arg(on_error); missing_coordinate <- match.arg(missing_coordinate)
  provenance <- .cater_scmore_provenance(strict_upstream=strict_upstream)
  .cater_scmore_validate_seurat(single_cell)

  if (is.null(gene_annotation)) {
    gene_annotation <- .cater_scmore_annotation_from_object(single_cell)
    annotation_source <- "Signac::Annotation(single_cell[['peaks']])"
  } else {
    gene_annotation <- .cater_scmore_standardize_gene_annotation(gene_annotation)
    annotation_source <- "user gene_annotation"
  }

  labels <- .cater_scmore_cell_labels(single_cell,celltype_col=celltype_col)
  label_source <- attr(labels,"source")
  available <- unique(unname(labels))
  if (is.null(cell_types)) cell_types <- available else {
    cell_types <- unique(as.character(cell_types)); unknown <- setdiff(cell_types,available)
    if (length(unknown)) .cater_scmore_stop("Requested cell_types are absent: %s",paste(unknown,collapse=", "))
  }
  if (!length(cell_types)) .cater_scmore_stop("No cell types selected")
  output_names <- setNames(.cater_scmore_output_names(cell_types),cell_types)

  create_args <- .cater_scmore_create_args(n_targets=n_targets,peak2gene_method=peak2gene_method,infer_method=infer_method,tss_upstream=tss_upstream,tss_downstream=tss_downstream,exclude_exon_regions=exclude_exon_regions,conserved_regions=conserved_regions)
  if (!is.null(outdir)) dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

  raw_outputs <- setNames(vector("list",length(cell_types)),cell_types)
  edge_evidence <- setNames(vector("list",length(cell_types)),cell_types)
  grns <- setNames(vector("list",length(cell_types)),cell_types)
  summary_rows <- vector("list",length(cell_types))

  for (i in seq_along(cell_types)) {
    ct <- cell_types[[i]]; cells <- names(labels)[labels==ct]
    .cater_scmore_msg(verbose,"[scMORE GRN] %s: %d cells",ct,length(cells))
    ct_object <- base::subset(single_cell,cells=cells)
    fit <- tryCatch(.cater_scmore_call_create_regulon(ct_object,create_args),error=function(e)e)
    if (inherits(fit,"error")) {
      if (on_error=="stop") .cater_scmore_stop("scMORE GRN failed for cell type '%s': %s",ct,conditionMessage(fit))
      summary_rows[[i]] <- data.frame(cell_type=ct,n_cells=length(cells),n_raw_rows=NA_integer_,n_unique_scMORE_edges=NA_integer_,n_cater_edges=NA_integer_,n_tfs=NA_integer_,n_targets=NA_integer_,n_self_loops_dropped=NA_integer_,n_missing_coordinate_edges=NA_integer_,status="SCMORE_FAILED",error_message=.cater_scmore_clean_message(conditionMessage(fit)),stringsAsFactors=FALSE)
      rm(ct_object,fit); if(isTRUE(gc_after_each)) gc(verbose=FALSE); next
    }

    processed <- tryCatch({
      .cater_scmore_validate_output(fit,ct)
      raw <- fit$grn; raw$cell_type <- rep(ct,nrow(raw))
      g <- .cater_scmore_make_cater_grn(fit$grn,gene_annotation,ct,drop_self_loops=drop_self_loops,missing_coordinate=missing_coordinate)
      list(raw=raw,g=g)
    },error=function(e)e)

    if (inherits(processed,"error")) {
      if (on_error=="stop") .cater_scmore_stop("CATER GRN contract failed for cell type '%s': %s",ct,conditionMessage(processed))
      raw_outputs[[ct]] <- fit
      summary_rows[[i]] <- data.frame(cell_type=ct,n_cells=length(cells),n_raw_rows=nrow(fit$grn),n_unique_scMORE_edges=length(unique(paste(fit$grn$TF,fit$grn$Target,sep="\r"))),n_cater_edges=NA_integer_,n_tfs=NA_integer_,n_targets=NA_integer_,n_self_loops_dropped=NA_integer_,n_missing_coordinate_edges=NA_integer_,status="CATER_CONTRACT_FAILED",error_message=.cater_scmore_clean_message(conditionMessage(processed)),stringsAsFactors=FALSE)
      rm(ct_object,fit,processed); if(isTRUE(gc_after_each)) gc(verbose=FALSE); next
    }

    raw_outputs[[ct]] <- fit; edge_evidence[[ct]] <- processed$raw; grns[[ct]] <- processed$g
    g <- processed$g
    status <- if(nrow(g)) "OK" else "EMPTY_GRN"
    summary_rows[[i]] <- data.frame(cell_type=ct,n_cells=length(cells),n_raw_rows=nrow(fit$grn),n_unique_scMORE_edges=attr(g,"n_unique_raw_edges"),n_cater_edges=nrow(g),n_tfs=length(unique(g$TF)),n_targets=length(unique(g$Target)),n_self_loops_dropped=attr(g,"n_self_loops_dropped"),n_missing_coordinate_edges=attr(g,"n_missing_coordinate_edges"),status=status,error_message=NA_character_,stringsAsFactors=FALSE)

    if (!is.null(outdir)) {
      safe <- output_names[[ct]]
      saveRDS(fit,file.path(outdir,paste0("scMORE_raw_",safe,".rds")))
      utils::write.table(processed$raw,file.path(outdir,paste0("scMORE_edge_evidence_",safe,".tsv")),sep="\t",quote=TRUE,row.names=FALSE)
      utils::write.table(g,file.path(outdir,paste0("CATER_grn_",safe,".tsv")),sep="\t",quote=TRUE,row.names=FALSE)
    }
    rm(ct_object,fit,processed,g); if(isTRUE(gc_after_each)) gc(verbose=FALSE)
  }

  fit_summary <- do.call(rbind,summary_rows); rownames(fit_summary)<-NULL
  if (!is.null(outdir)) {
    utils::write.table(gene_annotation,file.path(outdir,"CATER_gene_annotation.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
    utils::write.table(fit_summary,file.path(outdir,"scMORE_grn_summary.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
    saveRDS(provenance,file.path(outdir,"scMORE_provenance.rds"))
  }

  structure(list(grns=grns,edge_evidence=edge_evidence,scmore_outputs=raw_outputs,gene_annotation=gene_annotation,summary=fit_summary,provenance=provenance,celltype_source=label_source,gene_annotation_source=annotation_source,createRegulon_args=create_args,mode="per_celltype_scMORE_refit_then_CATER_direct_onehop_contract"),class="cater_scmore_grn")
}

#' Extract one CATER-ready cell-type GRN
cater_get_scmore_grn <- function(x,cell_type=NULL,raw=FALSE,evidence=FALSE) {
  if (!inherits(x,"cater_scmore_grn")) .cater_scmore_stop("x must be a cater_scmore_grn object")
  if (isTRUE(raw)&&isTRUE(evidence)) .cater_scmore_stop("Choose at most one of raw=TRUE or evidence=TRUE")
  available <- names(x$grns)[!vapply(x$grns,is.null,logical(1))]
  if (is.null(cell_type)) {
    if (length(available)!=1L) .cater_scmore_stop("Specify cell_type; available fits: %s",paste(available,collapse=", "))
    cell_type <- available[[1L]]
  }
  if (!cell_type%in%names(x$grns)||is.null(x$grns[[cell_type]])) .cater_scmore_stop("No successful CATER-ready scMORE GRN for cell type '%s'",cell_type)
  if(isTRUE(raw)) return(x$scmore_outputs[[cell_type]])
  if(isTRUE(evidence)) return(x$edge_evidence[[cell_type]])
  x$grns[[cell_type]]
}

print.cater_scmore_grn <- function(x,...) {
  ok <- sum(x$summary$status%in%c("OK","EMPTY_GRN"),na.rm=TRUE)
  cat("CATER-MR scMORE cell-type GRNs\n")
  cat("  mode:",x$mode,"\n")
  cat("  completed:",ok,"/",nrow(x$summary),"cell types\n")
  cat("  CATER-ready edges:",sum(x$summary$n_cater_edges,na.rm=TRUE),"\n")
  cat("  scMORE version:",x$provenance$scMORE_version,"\n")
  cat("  scMORE RemoteSha:",x$provenance$scMORE_remote_sha,"\n")
  invisible(x)
}
