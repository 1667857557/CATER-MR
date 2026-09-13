# CATER-MR scMORE GRN adapter v0.6
#
# scMORE performs the regulatory-network inference. CATER-MR only defines the
# cell population first, calls audited scMORE::createRegulon(), preserves all raw
# edge evidence, and converts it to the direct one-hop graph contract consumed by MR.

.CATER_SCMORE_REPO <- "mayunlong89/scMORE"
.CATER_SCMORE_AUDITED_SHA <- "f614736b9f49631471b0f18dfa415ee35e4d7b66"
.CATER_SCMORE_AUDITED_VERSION <- "2.0.0"

.cater_scmore_stop <- function(...) stop(sprintf(...), call.=FALSE)
.cater_scmore_msg <- function(verbose, ...) if(isTRUE(verbose)) message(sprintf(...))

.cater_scmore_provenance <- function(strict_upstream=TRUE) {
  if(!requireNamespace("scMORE",quietly=TRUE)) .cater_scmore_stop("scMORE is required; install the audited commit %s",.CATER_SCMORE_AUDITED_SHA)
  d<-utils::packageDescription("scMORE");sha<-d[["RemoteSha"]];user<-d[["RemoteUsername"]];repo<-d[["RemoteRepo"]]
  verified<-!is.null(sha)&&!is.null(user)&&!is.null(repo)&&identical(user,"mayunlong89")&&identical(repo,"scMORE")&&identical(tolower(sha),tolower(.CATER_SCMORE_AUDITED_SHA))
  if(isTRUE(strict_upstream)&&!verified) .cater_scmore_stop("strict_upstream=TRUE requires scMORE commit %s; installed RemoteSha=%s",.CATER_SCMORE_AUDITED_SHA,if(is.null(sha))"unavailable" else sha)
  pv<-function(p) if(requireNamespace(p,quietly=TRUE))as.character(utils::packageVersion(p)) else NA_character_
  list(scMORE_repository=.CATER_SCMORE_REPO,scMORE_audited_sha=.CATER_SCMORE_AUDITED_SHA,scMORE_remote_sha=if(is.null(sha))NA_character_ else sha,
       scMORE_version=as.character(utils::packageVersion("scMORE")),scMORE_sha_verified=verified,Pando_version=pv("Pando"),Seurat_version=pv("Seurat"),Signac_version=pv("Signac"),
       celltype_refit_status="CATER_ADAPTATION_NOT_SCMORE_TOPLEVEL_DEFAULT")
}

.cater_scmore_check_genome_values <- function(genome_values,source="annotation") {
  z<-unique(as.character(genome_values));z<-z[!is.na(z)&nzchar(z)]
  if(!length(z))return("UNKNOWN")
  ok<-grepl("^(hg38|GRCh38)([._-].*)?$",z,ignore.case=TRUE)
  if(!all(ok)) .cater_scmore_stop("%s explicitly reports a non-hg38 genome: %s",source,paste(z[!ok],collapse=", "))
  paste(sort(unique(z)),collapse=";")
}

.cater_scmore_validate_seurat <- function(single_cell) {
  if(!inherits(single_cell,"Seurat")) .cater_scmore_stop("single_cell must be a processed Seurat object")
  miss<-setdiff(c("RNA","peaks"),names(single_cell@assays));if(length(miss)) .cater_scmore_stop("scMORE::createRegulon() requires assays RNA and peaks; missing: %s",paste(miss,collapse=", "))
  if(!requireNamespace("Signac",quietly=TRUE)) .cater_scmore_stop("Signac is required")
  anno<-tryCatch(Signac::Annotation(single_cell[["peaks"]]),error=function(e)NULL)
  if(is.null(anno)||length(anno)==0L) .cater_scmore_stop("The peaks ChromatinAssay has no gene annotation")
  if(requireNamespace("GenomeInfoDb",quietly=TRUE)) .cater_scmore_check_genome_values(tryCatch(GenomeInfoDb::genome(anno),error=function(e)character()),"Signac peak annotation")
  invisible(TRUE)
}

.cater_scmore_cell_labels <- function(single_cell,celltype_col=NULL) {
  cells<-colnames(single_cell);if(!length(cells)) .cater_scmore_stop("single_cell has no cells")
  if(is.null(celltype_col)){
    if(!requireNamespace("SeuratObject",quietly=TRUE)) .cater_scmore_stop("SeuratObject is required")
    x<-as.character(SeuratObject::Idents(single_cell));names(x)<-cells;src<-"Idents(single_cell)"
  }else{
    md<-single_cell[[]];if(!celltype_col%in%names(md)) .cater_scmore_stop("celltype_col '%s' is absent",celltype_col)
    x<-as.character(md[[celltype_col]]);names(x)<-rownames(md);x<-x[match(cells,names(x))];names(x)<-cells;src<-paste0("metadata$",celltype_col)
  }
  if(any(is.na(x)|!nzchar(x))) .cater_scmore_stop("Every cell must have a nonempty cell-type label")
  attr(x,"source")<-src;x
}

.cater_scmore_default_regions <- function() {
  if(!requireNamespace("Pando",quietly=TRUE)) .cater_scmore_stop("Pando is required to load the default hg38 regulatory regions")
  if(!requireNamespace("GenomicRanges",quietly=TRUE)) .cater_scmore_stop("GenomicRanges is required to combine the default hg38 regulatory regions")
  e<-new.env(parent=baseenv())
  suppressWarnings(utils::data(list=c("phastConsElements20Mammals.UCSC.hg38","SCREEN.ccRE.UCSC.hg38"),package="Pando",envir=e))
  miss<-setdiff(c("phastConsElements20Mammals.UCSC.hg38","SCREEN.ccRE.UCSC.hg38"),ls(e,all.names=TRUE))
  if(length(miss)) .cater_scmore_stop("Pando is missing required hg38 region dataset(s): %s",paste(miss,collapse=", "))
  GenomicRanges::union(e$phastConsElements20Mammals.UCSC.hg38,e$SCREEN.ccRE.UCSC.hg38)
}

.cater_scmore_create_args <- function(n_targets=5,peak2gene_method="Signac",infer_method="glm",tss_upstream=100000,tss_downstream=0,exclude_exon_regions=TRUE,conserved_regions=NULL) {
  x<-list(n_targets=n_targets,peak2gene_method=peak2gene_method,infer_method=infer_method,tss_upstream=tss_upstream,tss_downstream=tss_downstream,exclude_exon_regions=exclude_exon_regions)
  if(!is.null(conserved_regions))x$conserved_regions<-conserved_regions;x
}

.cater_scmore_call_create_regulon <- function(single_cell,args,create_fun=NULL) {
  if(is.null(create_fun))create_fun<-getExportedValue("scMORE","createRegulon");do.call(create_fun,c(list(single_cell=single_cell),args))
}

.cater_scmore_validate_output <- function(x,cell_type) {
  if(!is.list(x)||is.null(x$grn)||is.null(x$tf_names)) .cater_scmore_stop("Unexpected scMORE output for '%s'",cell_type)
  if(!is.data.frame(x$grn)||!all(c("TF","Target")%in%names(x$grn))) .cater_scmore_stop("scMORE GRN for '%s' lacks TF/Target columns",cell_type)
  invisible(TRUE)
}

.cater_scmore_pick_col <- function(x,candidates,required=TRUE,label="column") {
  nm<-names(x);hit<-match(tolower(candidates),tolower(nm),nomatch=0L);hit<-hit[hit>0L]
  if(length(hit))return(nm[hit[1L]]);if(required).cater_scmore_stop("Missing %s; tried: %s",label,paste(candidates,collapse=", "));NULL
}

.cater_scmore_standardize_gene_annotation <- function(gene_annotation) {
  if(!is.data.frame(gene_annotation)) .cater_scmore_stop("gene_annotation must be a data.frame")
  sym<-.cater_scmore_pick_col(gene_annotation,c("symbol","gene","gene_symbol","SYMBOL"),label="gene symbol")
  chr<-.cater_scmore_pick_col(gene_annotation,c("chr","CHR","chrom","chromosome"),label="chromosome")
  tss<-.cater_scmore_pick_col(gene_annotation,c("tss","TSS"),label="gene-level TSS")
  genome_col<-.cater_scmore_pick_col(gene_annotation,c("genome","assembly","genome_build"),required=FALSE)
  if(!is.null(genome_col)).cater_scmore_check_genome_values(gene_annotation[[genome_col]],"gene_annotation")
  out<-data.frame(symbol=as.character(gene_annotation[[sym]]),chr=sub("^chr","",as.character(gene_annotation[[chr]]),ignore.case=TRUE),tss=suppressWarnings(as.numeric(gene_annotation[[tss]])),stringsAsFactors=FALSE)
  out<-out[!is.na(out$symbol)&nzchar(out$symbol)&!is.na(out$chr)&nzchar(out$chr)&is.finite(out$tss),,drop=FALSE]
  key<-paste(out$chr,format(out$tss,scientific=FALSE,trim=TRUE),sep=":");amb<-vapply(split(key,out$symbol),function(z)length(unique(z))>1L,logical(1))
  if(any(amb)){bad<-names(amb)[amb];.cater_scmore_stop("gene_annotation has multiple chr/TSS values for %d symbols (e.g. %s)",length(bad),paste(utils::head(bad,10L),collapse=", "))}
  out<-out[!duplicated(out$symbol),,drop=FALSE];rownames(out)<-NULL;out
}

.cater_scmore_annotation_from_object <- function(single_cell) {
  anno<-Signac::Annotation(single_cell[["peaks"]])
  tssgr<-Signac::GetTSSPositions(anno,biotypes=NULL);ad<-as.data.frame(tssgr)
  if(!nrow(ad)) .cater_scmore_stop("Signac::GetTSSPositions returned no genes")
  sym<-.cater_scmore_pick_col(ad,c("gene_name","gene","symbol","gene_symbol","SYMBOL"),label="gene symbol in TSS annotation")
  seqcol<-.cater_scmore_pick_col(ad,c("seqnames","chr","chromosome"),label="TSS chromosome")
  startcol<-.cater_scmore_pick_col(ad,c("start"),label="TSS coordinate")
  out<-data.frame(symbol=as.character(ad[[sym]]),chr=sub("^chr","",as.character(ad[[seqcol]]),ignore.case=TRUE),tss=suppressWarnings(as.numeric(ad[[startcol]])),stringsAsFactors=FALSE)
  .cater_scmore_standardize_gene_annotation(out)
}

.cater_scmore_make_cater_grn <- function(raw_grn,gene_annotation,cell_type,drop_self_loops=TRUE,missing_coordinate=c("error","drop")) {
  missing_coordinate<-match.arg(missing_coordinate)
  if(!is.data.frame(raw_grn)||!all(c("TF","Target")%in%names(raw_grn))) .cater_scmore_stop("raw_grn must contain TF and Target")
  edges<-unique(data.frame(TF=as.character(raw_grn$TF),Target=as.character(raw_grn$Target),stringsAsFactors=FALSE));edges<-edges[!is.na(edges$TF)&nzchar(edges$TF)&!is.na(edges$Target)&nzchar(edges$Target),,drop=FALSE]
  nraw<-nrow(edges);nself<-sum(edges$TF==edges$Target);if(isTRUE(drop_self_loops)&&nself)edges<-edges[edges$TF!=edges$Target,,drop=FALSE]
  template<-data.frame(TF=character(),Target=character(),TF_chr=character(),TF_tss=numeric(),Target_chr=character(),Target_tss=numeric(),cell_type=character(),stringsAsFactors=FALSE)
  if(!nrow(edges)){attr(template,"n_unique_raw_edges")<-nraw;attr(template,"n_self_loops_dropped")<-if(isTRUE(drop_self_loops))nself else 0L;attr(template,"n_missing_coordinate_edges")<-0L;return(template)}
  ann<-.cater_scmore_standardize_gene_annotation(gene_annotation);itf<-match(edges$TF,ann$symbol);itg<-match(edges$Target,ann$symbol);miss<-is.na(itf)|is.na(itg)
  if(any(miss)&&missing_coordinate=="error"){nodes<-unique(c(edges$TF[is.na(itf)],edges$Target[is.na(itg)]));.cater_scmore_stop("CATER-ready GRN '%s' lacks chr/TSS for %d nodes (e.g. %s)",cell_type,length(nodes),paste(utils::head(nodes,10L),collapse=", "))}
  nmiss<-sum(miss);if(any(miss)){edges<-edges[!miss,,drop=FALSE];itf<-itf[!miss];itg<-itg[!miss]}
  out<-data.frame(TF=edges$TF,Target=edges$Target,TF_chr=ann$chr[itf],TF_tss=ann$tss[itf],Target_chr=ann$chr[itg],Target_tss=ann$tss[itg],cell_type=rep(cell_type,nrow(edges)),stringsAsFactors=FALSE)
  attr(out,"n_unique_raw_edges")<-nraw;attr(out,"n_self_loops_dropped")<-if(isTRUE(drop_self_loops))nself else 0L;attr(out,"n_missing_coordinate_edges")<-nmiss;out
}

.cater_scmore_safe_name <- function(x){y<-gsub("[^A-Za-z0-9._-]+","_",x);y<-gsub("^_+|_+$","",y);ifelse(nzchar(y),y,"celltype")}
.cater_scmore_output_names <- function(cell_types){base<-vapply(cell_types,.cater_scmore_safe_name,character(1));rank<-match(cell_types,sort(unique(cell_types)));paste0(base,"__",sprintf("%03d",rank))}
.cater_scmore_clean_message <- function(x){if(!length(x)||is.na(x))return(NA_character_);trimws(gsub("[\r\n\t]+"," ",as.character(x)))}

cater_build_scmore_grn <- function(single_cell,celltype_col=NULL,cell_types=NULL,gene_annotation=NULL,n_targets=5,peak2gene_method="Signac",infer_method="glm",tss_upstream=100000,tss_downstream=0,exclude_exon_regions=TRUE,conserved_regions=.cater_scmore_default_regions(),drop_self_loops=TRUE,missing_coordinate=c("error","drop"),strict_upstream=TRUE,on_error=c("stop","record"),outdir=NULL,gc_after_each=TRUE,verbose=TRUE) {
  on_error<-match.arg(on_error);missing_coordinate<-match.arg(missing_coordinate);prov<-.cater_scmore_provenance(strict_upstream);.cater_scmore_validate_seurat(single_cell)
  if(is.null(gene_annotation)){gene_annotation<-.cater_scmore_annotation_from_object(single_cell);annsrc<-"Signac::GetTSSPositions(longest transcript)"}else{gene_annotation<-.cater_scmore_standardize_gene_annotation(gene_annotation);annsrc<-"user gene-level TSS annotation"}
  labels<-.cater_scmore_cell_labels(single_cell,celltype_col);src<-attr(labels,"source");available<-unique(unname(labels));if(is.null(cell_types))cell_types<-available else{cell_types<-unique(as.character(cell_types));bad<-setdiff(cell_types,available);if(length(bad)).cater_scmore_stop("Requested cell types absent: %s",paste(bad,collapse=", "))}
  args<-.cater_scmore_create_args(n_targets,peak2gene_method,infer_method,tss_upstream,tss_downstream,exclude_exon_regions,conserved_regions);if(!is.null(outdir))dir.create(outdir,recursive=TRUE,showWarnings=FALSE);namesafe<-setNames(.cater_scmore_output_names(cell_types),cell_types)
  raws<-edges<-grns<-setNames(vector("list",length(cell_types)),cell_types);rows<-vector("list",length(cell_types))
  for(i in seq_along(cell_types)){
    ct<-cell_types[[i]];cells<-names(labels)[labels==ct];.cater_scmore_msg(verbose,"[scMORE GRN] %s: %d cells",ct,length(cells));obj<-base::subset(single_cell,cells=cells)
    fit<-tryCatch(.cater_scmore_call_create_regulon(obj,args),error=function(e)e)
    if(inherits(fit,"error")){if(on_error=="stop").cater_scmore_stop("scMORE GRN failed for '%s': %s",ct,conditionMessage(fit));rows[[i]]<-data.frame(cell_type=ct,n_cells=length(cells),n_raw_rows=NA,n_unique_scMORE_edges=NA,n_cater_edges=NA,n_tfs=NA,n_targets=NA,n_self_loops_dropped=NA,n_missing_coordinate_edges=NA,status="SCMORE_FAILED",error_message=.cater_scmore_clean_message(conditionMessage(fit)));next}
    .cater_scmore_validate_output(fit,ct);raw<-fit$grn;raw$cell_type<-rep(ct,nrow(raw));g<-tryCatch(.cater_scmore_make_cater_grn(fit$grn,gene_annotation,ct,drop_self_loops,missing_coordinate),error=function(e)e)
    if(inherits(g,"error")){if(on_error=="stop").cater_scmore_stop("CATER GRN contract failed for '%s': %s",ct,conditionMessage(g));raws[[ct]]<-fit;edges[[ct]]<-raw;rows[[i]]<-data.frame(cell_type=ct,n_cells=length(cells),n_raw_rows=nrow(fit$grn),n_unique_scMORE_edges=length(unique(paste(fit$grn$TF,fit$grn$Target))),n_cater_edges=NA,n_tfs=NA,n_targets=NA,n_self_loops_dropped=NA,n_missing_coordinate_edges=NA,status="CATER_CONTRACT_FAILED",error_message=.cater_scmore_clean_message(conditionMessage(g)));next}
    raws[[ct]]<-fit;edges[[ct]]<-raw;grns[[ct]]<-g;rows[[i]]<-data.frame(cell_type=ct,n_cells=length(cells),n_raw_rows=nrow(fit$grn),n_unique_scMORE_edges=attr(g,"n_unique_raw_edges"),n_cater_edges=nrow(g),n_tfs=length(unique(g$TF)),n_targets=length(unique(g$Target)),n_self_loops_dropped=attr(g,"n_self_loops_dropped"),n_missing_coordinate_edges=attr(g,"n_missing_coordinate_edges"),status=if(nrow(g))"OK" else "EMPTY_GRN",error_message=NA_character_)
    if(!is.null(outdir)){safe<-namesafe[[ct]];saveRDS(fit,file.path(outdir,paste0("scMORE_raw_",safe,".rds")));utils::write.table(raw,file.path(outdir,paste0("scMORE_edge_evidence_",safe,".tsv")),sep="\t",quote=TRUE,row.names=FALSE);utils::write.table(g,file.path(outdir,paste0("CATER_grn_",safe,".tsv")),sep="\t",quote=TRUE,row.names=FALSE)}
    rm(obj,fit,raw,g);if(isTRUE(gc_after_each))gc(verbose=FALSE)
  }
  sm<-do.call(rbind,rows);rownames(sm)<-NULL
  if(!is.null(outdir)){utils::write.table(gene_annotation,file.path(outdir,"CATER_gene_annotation.tsv"),sep="\t",quote=TRUE,row.names=FALSE);utils::write.table(sm,file.path(outdir,"scMORE_grn_summary.tsv"),sep="\t",quote=TRUE,row.names=FALSE);saveRDS(prov,file.path(outdir,"scMORE_provenance.rds"))}
  structure(list(grns=grns,edge_evidence=edges,scmore_outputs=raws,gene_annotation=gene_annotation,summary=sm,provenance=prov,celltype_source=src,gene_annotation_source=annsrc,createRegulon_args=args,mode="per_celltype_scMORE_refit_CATER_ADAPTATION"),class="cater_scmore_grn")
}

cater_get_scmore_grn <- function(x,cell_type=NULL,raw=FALSE,evidence=FALSE){if(!inherits(x,"cater_scmore_grn")).cater_scmore_stop("x must be a cater_scmore_grn object");if(raw&&evidence).cater_scmore_stop("Choose only one of raw/evidence");ok<-names(x$grns)[!vapply(x$grns,is.null,logical(1))];if(is.null(cell_type)){if(length(ok)!=1L).cater_scmore_stop("Specify cell_type; available: %s",paste(ok,collapse=", "));cell_type<-ok[[1L]]};if(!cell_type%in%names(x$grns)||is.null(x$grns[[cell_type]])).cater_scmore_stop("No successful GRN for '%s'",cell_type);if(raw)return(x$scmore_outputs[[cell_type]]);if(evidence)return(x$edge_evidence[[cell_type]]);x$grns[[cell_type]]}

print.cater_scmore_grn <- function(x,...){cat("CATER-MR scMORE cell-type GRNs\n  mode:",x$mode,"\n  completed:",sum(x$summary$status%in%c("OK","EMPTY_GRN"),na.rm=TRUE),"/",nrow(x$summary),"\n  edges:",sum(x$summary$n_cater_edges,na.rm=TRUE),"\n");invisible(x)}
