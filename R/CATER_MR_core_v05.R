# CATER-MR v0.7
# Cis And Trans eQTLs guided by Regulatory networks for drug-target MR
# Direct implementation: GRN-gated cis/trans selection + standard LD-aware MR.

.CATER_VERSION <- "0.7.2"

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


.cater_standardize_grn <- function(grn) {
  if (is.list(grn) && !is.data.frame(grn) && !is.null(grn$grn)) grn <- grn$grn
  if (!is.data.frame(grn)) .cater_stop("grn must be a data.frame or a list containing $grn")
  tf <- .cater_pick_col(grn, c("TF","tf","regulator"), label="GRN TF")
  tg <- .cater_pick_col(grn, c("Target","target","gene"), label="GRN target")
  out <- data.frame(TF=as.character(grn[[tf]]), Target=as.character(grn[[tg]]),
                    stringsAsFactors=FALSE)
  out <- out[!is.na(out$TF)&!is.na(out$Target)&nzchar(out$TF)&nzchar(out$Target),,drop=FALSE]
  unique(out)
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


.cater_make_regions <- function(target, parents, annotation, cis_window, tf_window) {
  ta <- annotation[annotation$symbol==target,,drop=FALSE]
  if (!nrow(ta)) return(NULL)
  reg <- data.frame(type="cis",gene=target,chr=ta$chr[1],
                    start=max(1,ta$tss[1]-cis_window),end=ta$tss[1]+cis_window,
                    locus_id=NA_character_,stringsAsFactors=FALSE)
  if (length(parents)) {
    pa <- annotation[match(parents,annotation$symbol),,drop=FALSE]
    pa$gene <- parents
    pa <- pa[!is.na(pa$chr)&is.finite(pa$tss),,drop=FALSE]
    if (nrow(pa)) reg <- rbind(reg,data.frame(type="trans",gene=pa$gene,chr=pa$chr,
      start=pmax(1,pa$tss-tf_window),end=pa$tss+tf_window,
      locus_id=NA_character_,stringsAsFactors=FALSE))
  }

  # Define physical loci from all target-cis and parent-TF windows together.
  # Any connected interval component containing the target cis window is the
  # cis locus in full. This prevents a TF window that directly overlaps the
  # target cis region (and any TF window connected through it) from generating
  # an artificial trans instrument set inside the same local LD neighbourhood.
  component <- integer(nrow(reg)); next_component <- 0L
  for (ch in unique(reg$chr)) {
    oi <- which(reg$chr==ch)
    oi <- oi[order(reg$start[oi],reg$end[oi],reg$type[oi],reg$gene[oi])]
    current <- 0L; current_end <- -Inf
    for (idx in oi) {
      if (current==0L || reg$start[idx] > current_end) {
        next_component <- next_component + 1L
        current <- next_component
        current_end <- reg$end[idx]
      } else {
        current_end <- max(current_end,reg$end[idx])
      }
      component[idx] <- current
    }
  }
  for (cc in unique(component)) {
    jj <- which(component==cc)
    if (any(reg$type[jj]=="cis")) {
      reg$locus_id[jj] <- "cis"
    } else {
      reg$locus_id[jj] <- paste0("TF:",paste(sort(unique(reg$gene[jj])),collapse="+"))
    }
  }
  rownames(reg) <- NULL
  reg
}

.cater_candidate_map <- function(qtl, regions) {
  if (is.null(regions)||!nrow(qtl)) return(data.frame())
  if (!"locus_id" %in% names(regions) || any(is.na(regions$locus_id)|!nzchar(regions$locus_id)))
    .cater_stop("Candidate regions require a physical locus_id")

  # The cis locus is the full connected component containing the target cis
  # window, including any overlapping parent-TF windows.
  cis_regions <- regions[regions$locus_id=="cis",,drop=FALSE]
  is_cis <- rep(FALSE,nrow(qtl))
  if (nrow(cis_regions)) for (k in seq_len(nrow(cis_regions))) {
    is_cis <- is_cis | (qtl$chr==cis_regions$chr[k] &
                        qtl$pos>=cis_regions$start[k] & qtl$pos<=cis_regions$end[k])
  }

  tr <- regions[regions$type=="trans" & regions$locus_id!="cis",,drop=FALSE]
  hits <- vector("list",nrow(qtl)); loci <- vector("list",nrow(qtl))
  if (nrow(tr)) for (k in seq_len(nrow(tr))) {
    z <- qtl$chr==tr$chr[k] & qtl$pos>=tr$start[k] & qtl$pos<=tr$end[k] & !is_cis
    if (any(z)) {
      hits[z] <- lapply(hits[z], function(x) c(x,tr$gene[k]))
      loci[z] <- lapply(loci[z], function(x) c(x,tr$locus_id[k]))
    }
  }
  is_trans <- lengths(hits)>0L
  keep <- is_cis|is_trans
  if (!any(keep)) return(data.frame())
  parents <- vapply(hits[keep],function(x) if(length(x)) paste(sort(unique(x)),collapse=";") else "",character(1))
  physical_locus <- vapply(loci[keep],function(x) {
    u <- sort(unique(x[nzchar(x)]))
    if (!length(u)) return("")
    if (length(u)!=1L) .cater_stop("A trans SNP was assigned to multiple physical TF loci")
    u
  },character(1))
  data.frame(snp=qtl$snp[keep],source=ifelse(is_cis[keep],"cis","trans"),
             parent_tf=parents,
             locus_id=ifelse(is_cis[keep],"cis",physical_locus),
             stringsAsFactors=FALSE)
}

.cater_write_cojo_ma <- function(qtl,path) {
  ok <- is.finite(qtl$eaf)&qtl$eaf>0&qtl$eaf<1&is.finite(qtl$n)&qtl$n>0
  if (!all(ok)) .cater_stop("Manc-COJO requires valid EAF/freq and N; %d rows invalid",sum(!ok))
  ma <- data.frame(SNP=qtl$snp,A1=qtl$a1,A2=qtl$a2,freq=qtl$eaf,
                   b=qtl$beta,se=qtl$se,p=qtl$p,N=qtl$n,stringsAsFactors=FALSE)
  utils::write.table(ma,path,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
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


.cater_run_cojo <- function(qtl,candidate_map,ld_bfile,manc_cojo_bin,cojo_p,
                            cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,
                            plink_bin="plink",enable_ld_diagnosis=TRUE,
                            ld_diag_loglr=2,ld_diag_abs_z=2) {
  if (length(ld_bfile)!=1L || !nzchar(ld_bfile)) .cater_stop("Current CATER-MR expects one COJO LD cohort")
  if (!all(file.exists(paste0(ld_bfile,c(".bed",".bim",".fam"))))) .cater_stop("ld_bfile must be a PLINK prefix: %s",ld_bfile)
  exe <- Sys.which(manc_cojo_bin); if (!nzchar(exe)) .cater_stop("Cannot find Manc-COJO executable '%s'",manc_cojo_bin)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE); .cater_remove_prefix_outputs(prefix)
  diag_res <- list(candidate_map=candidate_map,diagnostics=data.frame(),removed=character())
  if (isTRUE(enable_ld_diagnosis)) {
    diag_res <- .cater_ld_diagnosis_filter(qtl,candidate_map,ld_bfile,plink_bin,cojo_threads,prefix,
      loglr_cutoff=ld_diag_loglr,abs_z_cutoff=ld_diag_abs_z,verbose=verbose)
    candidate_map <- diag_res$candidate_map
    if (nrow(diag_res$diagnostics)) utils::write.table(diag_res$diagnostics,paste0(prefix,".ld_diagnosis.tsv"),
      sep="\t",quote=FALSE,row.names=FALSE)
  }
  finish <- function(selected=data.frame(),selected_step=selected,joint_removed=character(),ld=matrix(numeric(),0,0))
    list(selected=selected,selected_step=selected_step,joint_removed=joint_removed,ld=ld,
         ld_diagnosis=diag_res$diagnostics,ld_diagnosis_removed=diag_res$removed)
  if (!nrow(candidate_map)) return(finish())
  ma <- paste0(prefix,".sumstat"); extract <- paste0(prefix,".candidate.snplist")
  .cater_write_cojo_ma(qtl,ma); writeLines(unique(candidate_map$snp),extract)
  sp <- paste0(prefix,".select")
  args <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",extract,"--cojo-slct",
            "--cojo-p",format(cojo_p,scientific=TRUE),"--cojo-wind",as.character(as.integer(cojo_wind_kb)),
            "--cojo-collinear",as.character(cojo_collinear),"--thread-num",as.character(as.integer(cojo_threads)),"--out",sp)
  st <- system2(exe,args=.cater_quote_args(args),stdout=paste0(sp,".stdout"),stderr=paste0(sp,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO selection failed for %s",prefix)
  selected_step <- .cater_read_jma_ref(paste0(sp,".jma.cojo"))
  if (!nrow(selected_step)) return(finish(selected=selected_step,selected_step=selected_step))
  if (nrow(selected_step)==1L) return(finish(selected=selected_step,selected_step=selected_step,
    ld=matrix(1,1,1,dimnames=list(selected_step$snp,selected_step$snp))))
  jp <- paste0(prefix,".joint"); .cater_remove_prefix_outputs(jp)
  args2 <- c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",paste0(sp,".jma.cojo"),"2","header",
             "--cojo-joint","--thread-num",as.character(as.integer(cojo_threads)),"--output-all","--out",jp)
  st <- system2(exe,args=.cater_quote_args(args2),stdout=paste0(jp,".stdout"),stderr=paste0(jp,".stderr"))
  if (!identical(st,0L)) .cater_stop("Manc-COJO joint/LD failed for %s",prefix)
  joint_ref <- .cater_read_jma_ref(paste0(jp,".jma.cojo"))
  if (!nrow(joint_ref)) .cater_stop("Manc-COJO joint retained no SNPs for %s",prefix)
  if (length(setdiff(joint_ref$snp,selected_step$snp))) .cater_stop("Manc-COJO joint output contains unselected SNPs for %s",prefix)
  ld <- .cater_read_manc_ldr(paste0(jp,".ldr.cojo"),joint_ref$snp)
  finish(selected=joint_ref,selected_step=selected_step,joint_removed=setdiff(selected_step$snp,joint_ref$snp),ld=ld)
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


.cater_align_to_ld <- function(dat,ref) {
  rr<-ref[match(dat$snp,ref$snp),,drop=FALSE]
  ok<-!is.na(rr$snp); dat<-dat[ok,,drop=FALSE]; rr<-rr[ok,,drop=FALSE]
  if(!nrow(dat)) return(dat)
  same<-dat$a1==rr$ld_a1&dat$a2==rr$ld_a2
  swap<-dat$a1==rr$ld_a2&dat$a2==rr$ld_a1
  keep<-same|swap; dat<-dat[keep,,drop=FALSE]; rr<-rr[keep,,drop=FALSE]; swap<-swap[keep]
  if(!nrow(dat)) return(dat)
  dat$beta[swap]<--dat$beta[swap]
  dat$eaf[swap&is.finite(dat$eaf)]<-1-dat$eaf[swap&is.finite(dat$eaf)]
  dat$a1<-rr$ld_a1;dat$a2<-rr$ld_a2
  dat
}

.cater_is_palindromic <- function(a1,a2) paste0(a1,a2)%in%c("AT","TA","CG","GC")

.cater_harmonize <- function(exp,outcome,drop_palindromic=TRUE) {
  y<-outcome[match(exp$snp,outcome$snp),,drop=FALSE]
  keep<-!is.na(y$snp); exp<-exp[keep,,drop=FALSE];y<-y[keep,,drop=FALSE]
  if(!nrow(exp)) return(data.frame())
  same<-exp$a1==y$a1&exp$a2==y$a2; swap<-exp$a1==y$a2&exp$a2==y$a1
  keep<-same|swap
  if(drop_palindromic) keep<-keep&!.cater_is_palindromic(exp$a1,exp$a2)
  exp<-exp[keep,,drop=FALSE];y<-y[keep,,drop=FALSE];swap<-swap[keep]
  if(!nrow(exp)) return(data.frame())
  data.frame(snp=exp$snp,source=exp$source,parent_tf=exp$parent_tf,locus_id=exp$locus_id,
             bx=exp$beta,bx_se=exp$se,by=ifelse(swap,-y$beta,y$beta),by_se=y$se,
             stringsAsFactors=FALSE)
}

.cater_subset_ld <- function(ld,snps) {
  if(!length(snps)) return(matrix(numeric(),0,0))
  if(length(snps)==1L) return(matrix(1,1,1,dimnames=list(snps,snps)))
  if(!all(snps%in%rownames(ld))||!all(snps%in%colnames(ld))) .cater_stop("LD matrix missing selected SNPs")
  ld[snps,snps,drop=FALSE]
}

.cater_inv <- function(x,tol=1e-10) {
  if(!length(x)) return(NULL)
  ev<-eigen((x+t(x))/2,symmetric=TRUE,only.values=TRUE)$values
  if(any(!is.finite(ev))||min(ev)<=tol*max(ev)) return(NULL)
  solve(x)
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


.cater_cis_trans_het <- function(h,ld) {
  ci<-which(h$source=="cis"); tr<-which(h$source=="trans")
  if(!length(ci)||!length(tr)) return(c(z=NA_real_,p=NA_real_))
  dc<-h[ci,,drop=FALSE];dt<-h[tr,,drop=FALSE]
  Rc<-.cater_subset_ld(ld,dc$snp);Rt<-.cater_subset_ld(ld,dt$snp)
  fc<-.cater_givw(dc,Rc);ft<-.cater_givw(dt,Rt)
  if(fc$status!="OK"||ft$status!="OK") return(c(z=NA_real_,p=NA_real_))
  Oab<-diag(dc$by_se,nrow=length(ci))%*%ld[dc$snp,dt$snp,drop=FALSE]%*%
    diag(dt$by_se,nrow=length(tr))
  Oc<-diag(dc$by_se,nrow=length(ci))%*%Rc%*%diag(dc$by_se,nrow=length(ci))
  Ot<-diag(dt$by_se,nrow=length(tr))%*%Rt%*%diag(dt$by_se,nrow=length(tr))
  ic<-.cater_inv(Oc);it<-.cater_inv(Ot)
  if(is.null(ic)||is.null(it)) return(c(z=NA_real_,p=NA_real_))
  wc<-as.numeric(ic%*%dc$bx)/as.numeric(crossprod(dc$bx,ic%*%dc$bx))
  wt<-as.numeric(it%*%dt$bx)/as.numeric(crossprod(dt$bx,it%*%dt$bx))
  cv<-as.numeric(t(wc)%*%Oab%*%wt)
  vd<-fc$se^2+ft$se^2-2*cv
  if(!is.finite(vd)||vd<=0) return(c(z=NA_real_,p=NA_real_))
  z<-(fc$beta-ft$beta)/sqrt(vd)
  c(z=z,p=2*stats::pnorm(-abs(z)))
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


.cater_read_gene <- function(gene,eqtl_dir,qtl_n) {
  f<-file.path(eqtl_dir,paste0(gene,".txt.gz"))
  if(!file.exists(f)) return(NULL)
  .cater_standardize_sumstats(.cater_read_table(f),n_default=qtl_n,label=paste0(gene," eQTL"))
}

.cater_exact_effects <- function(gene,snps,eqtl_dir,qtl_n,ld_ref) {
  q<-.cater_read_gene(gene,eqtl_dir,qtl_n)
  if(is.null(q)) return(NULL)
  d<-q[match(snps,q$snp),,drop=FALSE];d<-d[!is.na(d$snp),,drop=FALSE]
  if(!nrow(d)) return(data.frame())
  .cater_align_to_ld(d,ld_ref)
}

.cater_omnibus <- function(beta,se,ld) {
  ok<-is.finite(beta)&is.finite(se)&se>0
  beta<-beta[ok];se<-se[ok]
  if(!length(beta)) return(c(Q=NA_real_,df=NA_real_,p=NA_real_))
  R<-ld[ok,ok,drop=FALSE]
  if(length(beta)==1L) {q<-(beta/se)^2;return(c(Q=q,df=1,p=stats::pchisq(q,1,lower.tail=FALSE)))}
  S<-diag(se,nrow=length(se))%*%R%*%diag(se,nrow=length(se));inv<-.cater_inv(S)
  if(is.null(inv)) return(c(Q=NA_real_,df=NA_real_,p=NA_real_))
  q<-as.numeric(crossprod(beta,inv%*%beta));df<-qr(R)$rank
  c(Q=q,df=df,p=stats::pchisq(q,df,lower.tail=FALSE))
}

.cater_tf_anchor <- function(target,h,cojo_ref,ld,eqtl_dir,qtl_n) {
  tr<-h[h$source=="trans",,drop=FALSE]
  if(!nrow(tr)) return(data.frame())
  out<-list()
  for(i in seq_len(nrow(tr))){
    tfs<-strsplit(tr$parent_tf[i],";",fixed=TRUE)[[1L]]
    for(tf in tfs){
      d<-.cater_exact_effects(tf,tr$snp[i],eqtl_dir,qtl_n,cojo_ref)
      if(is.null(d)) {
        out[[length(out)+1L]]<-data.frame(target=target,snp=tr$snp[i],tf=tf,beta_tf=NA,se_tf=NA,p_tf=NA,status="QTL_NOT_AVAILABLE")
      } else if(!nrow(d)) {
        out[[length(out)+1L]]<-data.frame(target=target,snp=tr$snp[i],tf=tf,beta_tf=NA,se_tf=NA,p_tf=NA,status="SNP_NOT_AVAILABLE")
      } else {
        out[[length(out)+1L]]<-data.frame(target=target,snp=tr$snp[i],tf=tf,beta_tf=d$beta[1],se_tf=d$se[1],
          p_tf=2*stats::pnorm(-abs(d$beta[1]/d$se[1])),status="OK")
      }
    }
  }
  do.call(rbind,out)
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
  tab$co_perturbation_detected<-is.finite(tab$q)&tab$q<sibling_fdr
  tab$interpretation<-ifelse(tab$co_perturbation_detected,"DETECTED_SIBLING_COPERTURBATION",ifelse(tab$testable,"NO_DETECTED_SIBLING_COPERTURBATION","NOT_TESTABLE"))
  tab$active<-tab$co_perturbation_detected
  active<-tab$sibling[tab$active]
  complete<-all(tab$complete)
  list(table=tab,active=active,co_perturbed=active,n_candidate=nrow(tab),n_incomplete=sum(!tab$complete),n_qtl_missing=qtl_missing,n_iv_missing=sum(tab$n_iv_missing),complete=complete,testable=any(tab$testable),status=if(!complete)"INCOMPLETE" else if(length(active))"COMPLETE_DETECTED_COPERTURBATION" else "COMPLETE_NO_DETECTED_COPERTURBATION")
}


.cater_validate_exposure_corr <- function(exposure_corr,exposure_names) {
  p<-length(exposure_names)
  if(is.null(exposure_corr)) {
    C<-diag(p);dimnames(C)<-list(exposure_names,exposure_names);return(C)
  }
  C<-as.matrix(exposure_corr)
  if(!is.numeric(C)||nrow(C)!=ncol(C)||is.null(rownames(C))||is.null(colnames(C))||
     !all(exposure_names%in%rownames(C))||!all(exposure_names%in%colnames(C)))
    .cater_stop("exposure_corr must be a named square numeric matrix covering all MVMR exposures")
  C<-C[exposure_names,exposure_names,drop=FALSE]
  if(any(!is.finite(C))) .cater_stop("exposure_corr contains non-finite values")
  if(max(abs(C-t(C)))>1e-8) .cater_stop("exposure_corr must be symmetric")
  if(any(abs(diag(C)-1)>1e-6)) .cater_stop("exposure_corr must have unit diagonal")
  if(any(abs(C)>1+1e-8)) .cater_stop("exposure_corr entries must lie in [-1,1]")
  ev<-eigen((C+t(C))/2,symmetric=TRUE,only.values=TRUE)$values
  if(min(ev) < -1e-8 * max(1,max(abs(ev)))) .cater_stop("exposure_corr must be positive semidefinite")
  C
}

.cater_residual_exposure_cov <- function(SE,ld,q,C) {
  M<-sweep(SE,2,q,"*")
  V<-(M%*%C%*%t(M))*ld
  (V+t(V))/2
}

.cater_conditional_f <- function(B,SE,ld,exposure_names,exposure_corr=NULL,
                                  max_iter=200L,tol=1e-9) {
  B<-as.matrix(B);SE<-as.matrix(SE);ld<-as.matrix(ld)
  p<-ncol(B);m<-nrow(B)
  if(length(exposure_names)!=p) .cater_stop("exposure_names must match the MVMR exposure columns")
  if(!all(dim(SE)==dim(B))) .cater_stop("SE must have the same dimensions as B")
  if(!all(dim(ld)==c(m,m))) .cater_stop("LD must be an m x m signed-correlation matrix")
  if(any(!is.finite(B))||any(!is.finite(SE))||any(SE<=0))
    .cater_stop("B/SE supplied to conditional F must be finite with SE > 0")
  if(any(!is.finite(ld))||max(abs(ld-t(ld)))>1e-8)
    .cater_stop("LD supplied to conditional F must be finite and symmetric")
  ans<-setNames(rep(NA_real_,p),exposure_names)
  df<-m-p+1L
  if(p<2L||df<=0L) return(ans)
  C<-.cater_validate_exposure_corr(exposure_corr,exposure_names)
  qstat<-setNames(rep(NA_real_,p),exposure_names)
  converged<-setNames(rep(FALSE,p),exposure_names)
  deltas<-setNames(vector("list",p),exposure_names)
  iterations<-setNames(rep(NA_integer_,p),exposure_names)
  for(i in seq_len(p)){
    other<-setdiff(seq_len(p),i);X<-B[,other,drop=FALSE];y<-B[,i]
    delta<-tryCatch(stats::lm.fit(x=X,y=y)$coefficients,
                    error=function(e) rep(0,length(other)))
    if(length(delta)!=length(other)||any(!is.finite(delta))) delta<-rep(0,length(other))
    ok<-FALSE
    for(iter in seq_len(max_iter)){
      q<-numeric(p);q[i]<-1;q[other]<--delta
      V<-.cater_residual_exposure_cov(SE,ld,q,C);W<-.cater_inv(V)
      if(is.null(W)) break
      A<-crossprod(X,W%*%X);Ai<-.cater_inv(A)
      if(is.null(Ai)) break
      delta_new<-as.numeric(Ai%*%crossprod(X,W%*%y))
      if(any(!is.finite(delta_new))) break
      scale<-1+max(abs(delta),abs(delta_new))
      if(max(abs(delta_new-delta))<=tol*scale){delta<-delta_new;ok<-TRUE;iterations[i]<-iter;break}
      delta<-0.5*delta+0.5*delta_new
    }
    if(!ok) next
    q<-numeric(p);q[i]<-1;q[other]<--delta
    V<-.cater_residual_exposure_cov(SE,ld,q,C);W<-.cater_inv(V)
    if(is.null(W)) next
    r<-y-as.numeric(X%*%delta);Q<-as.numeric(crossprod(r,W%*%r))
    if(!is.finite(Q)||Q<0) next
    qstat[i]<-Q;ans[i]<-Q/df;converged[i]<-TRUE;deltas[[i]]<-delta
  }
  attr(ans,"Q")<-qstat
  attr(ans,"df")<-setNames(rep(df,p),exposure_names)
  attr(ans,"converged")<-converged
  attr(ans,"delta")<-deltas
  attr(ans,"iterations")<-iterations
  attr(ans,"method") <- "CORRELATED_IV_CONDITIONAL_F"
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


.cater_run_sibling_cis <- function(gene,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
                                   cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,
                                   plink_bin="plink",enable_ld_diagnosis=TRUE,
                                   ld_diag_loglr=2,ld_diag_abs_z=2) {
  q<-.cater_read_gene(gene,eqtl_dir,qtl_n);if(is.null(q)) return(NULL)
  reg<-.cater_make_regions(gene,character(),annotation,cis_window,cis_window)
  if(is.null(reg)) return(NULL)
  cmap<-.cater_candidate_map(q,reg);if(!nrow(cmap)) return(list(qtl=q,selected=data.frame()))
  co<-.cater_run_cojo(q,cmap,ld_bfile,manc_cojo_bin,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose,
    plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z)
  list(qtl=q,selected=co$selected)
}

.cater_build_mvmr <- function(target,target_qtl,target_sel,active_siblings,annotation,eqtl_dir,qtl_n,
                              ld_bfile,manc_cojo_bin,cis_window,cojo_p,cojo_wind_kb,cojo_collinear,
                              cojo_threads,outdir,outcome,drop_palindromic,exposure_corr,min_cond_F,
                              mvmr_max_r2,mvmr_max_condition,verbose,plink_bin="plink",
                              enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2) {
  sibco<-list();union<-target_sel$snp
  for(z in active_siblings){
    pr<-file.path(outdir,"cojo_mvmr",target,z)
    x<-tryCatch(.cater_run_sibling_cis(z,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
      cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,pr,verbose,
      plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,
      ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z),error=function(e)NULL)
    sibco[[z]]<-x
    if(!is.null(x)&&nrow(x$selected)) union<-unique(c(union,x$selected$snp))
  }
  if(length(union)<2L) return(list(status="MVMR_NO_INSTRUMENT_UNION"))
  jl<-tryCatch(.cater_joint_ld(target_qtl,union,ld_bfile,manc_cojo_bin,cojo_threads,
    file.path(outdir,"cojo_mvmr",target,"union_ld"),verbose),error=function(e)NULL)
  if(is.null(jl)||!nrow(jl$selected)) return(list(status="MVMR_LD_FAILED"))
  snps<-jl$selected$snp; exposures<-c(target,active_siblings)
  B<-matrix(NA_real_,length(snps),length(exposures),dimnames=list(snps,exposures))
  SE<-B
  for(e in exposures){
    q<-if(e==target) target_qtl else .cater_read_gene(e,eqtl_dir,qtl_n)
    if(is.null(q)) next
    d<-q[match(snps,q$snp),,drop=FALSE];d<-d[!is.na(d$snp),,drop=FALSE]
    d<-.cater_align_to_ld(d,jl$selected)
    B[d$snp,e]<-d$beta;SE[d$snp,e]<-d$se
  }
  y<-outcome[match(snps,outcome$snp),,drop=FALSE]
  same<-!is.na(y$snp)&jl$selected$ld_a1==y$a1&jl$selected$ld_a2==y$a2
  swap<-!is.na(y$snp)&jl$selected$ld_a1==y$a2&jl$selected$ld_a2==y$a1
  ok<-(same|swap)&complete.cases(B)&complete.cases(SE)&is.finite(y$se)
  if(drop_palindromic) ok<-ok&!.cater_is_palindromic(jl$selected$ld_a1,jl$selected$ld_a2)
  if(sum(ok)<=length(exposures)) return(list(status="MVMR_UNDERIDENTIFIED"))
  by<-ifelse(swap[ok],-y$beta[ok],y$beta[ok]);seY<-y$se[ok]
  R<-jl$ld[snps[ok],snps[ok],drop=FALSE]
  fit<-.cater_mvmr_fit(B[ok,,drop=FALSE],SE[ok,,drop=FALSE],by,seY,R,exposures,exposure_corr)
  if(fit$status!="OK") return(fit)
  fit$n_iv<-sum(ok)
  fit$n_union_requested<-length(union)
  fit$n_union_ld_available<-length(snps)
  fit$n_union_complete_case<-sum(ok)
  off<-R;diag(off)<-0
  fit$max_r2<-if(length(off)) max(off^2,na.rm=TRUE) else 0
  ld_gate<-is.null(mvmr_max_r2)||(length(mvmr_max_r2)==1L&&is.finite(mvmr_max_r2)&&
    is.finite(fit$max_r2)&&fit$max_r2<=mvmr_max_r2)
  fit$ld_gate_pass<-ld_gate
  fit$primary_eligible<-is.finite(fit$conditional_F[target])&&fit$conditional_F[target]>=min_cond_F&&
    isTRUE(fit$conditional_F_converged[target])&&!is.null(exposure_corr)&&
    is.finite(fit$condition)&&fit$condition<=mvmr_max_condition&&ld_gate
  fit
}

.cater_trans_information_fraction <- function(cis_fit,combined_fit) {
  Ia<-combined_fit$information
  if(!is.finite(Ia)||Ia<=0) return(NA_real_)
  Ic<-cis_fit$information
  if(isTRUE(cis_fit$n_iv==0L)) Ic<-0
  if(!is.finite(Ic)) return(NA_real_)
  max(0,min(1,(Ia-Ic)/Ia))
}

.cater_primary_decision <- function(fits,net=NULL,target=NULL,has_trans=FALSE,
                                    sibling_screen_performed=FALSE,sibling_screen_complete=TRUE,
                                    n_active_siblings=0L,primary_policy=c("cis_anchor","screened_cater"),
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
      use_fit("cis",fits$cis,"TRANS_PLEIOTROPY_UNRESOLVED","DETECTED_SIBLING_COPERTURBATION")
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
    use_fit("combined",fits$combined,"OK_CATER","NO_DETECTED_SIBLING_COPERTURBATION")
  } else {
    # screened_cater is an opt-in promotion policy, not permission to discard a
    # valid conservative cis estimate when the augmented cis+trans fit fails.
    use_fit("cis",fits$cis,"CATER_COMBINED_FAILED_CIS_FALLBACK","COMBINED_FIT_FAILED")
  }
  out
}


#' Run CATER-MR
#'
#' Primary inputs: cell-type-specific one-hop GRN, directory of SYMBOL.txt.gz
#' full-summary eQTL files, and an in-memory outcome summary object.
.cater_atomic_write_table <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(.Machine$integer.max, 1L))
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)
  utils::write.table(x, tmp, ...)
  if (!file.rename(tmp, path)) .cater_stop("Could not atomically replace output: %s", path)
  invisible(path)
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

.cater_remove_prefix_outputs <- function(prefix) {
  z <- Sys.glob(paste0(prefix, "*")); if (length(z)) unlink(z, recursive=TRUE, force=TRUE)
  invisible(NULL)
}

.cater_quote_args <- function(args) {
  vapply(as.character(args), function(z) if (grepl("[[:space:]]",z)) shQuote(z) else z,
         character(1), USE.NAMES=FALSE)
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



.cater_descendants <- function(grn,target) {
  seen <- character(); frontier <- target
  repeat {
    nxt <- unique(grn$Target[grn$TF %in% frontier]); nxt <- setdiff(nxt,c(target,seen))
    if(!length(nxt)) break
    seen <- unique(c(seen,nxt)); frontier <- nxt
  }
  seen
}

.cater_validate_input_manifest <- function(input_manifest, strict=TRUE) {
  keys <- c("grn_build","eqtl_build","ld_build","outcome_build","eqtl_n_unit","eqtl_full_summary","eqtl_ancestry","outcome_ancestry","ld_ancestry","trans_qtl_qc","outcome_ld_ancestry")
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

`%||%` <- function(x,y) if(is.null(x)) y else x

cater_mr <- function(grn,eqtl_dir,outcome,gene_annotation=NULL,ld_bfile,
                     manc_cojo_bin="manc_cojo",targets=NULL,cell_type=NA_character_,trait=NA_character_,
                     cis_window=1e6,tf_window=cis_window,cojo_p=5e-8,cojo_wind_kb=10000L,
                     cojo_collinear=0.9,cojo_threads=1L,qtl_n=NULL,
                     sibling_fdr=0.05,enable_sibling_screen=TRUE,enable_mvmr=TRUE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2) {
  primary_policy<-match.arg(primary_policy)
  if(length(plink_bin)!=1L||is.na(plink_bin)||!nzchar(as.character(plink_bin))) .cater_stop("plink_bin must be a non-empty scalar")
  if(!is.logical(enable_ld_diagnosis)||length(enable_ld_diagnosis)!=1L||is.na(enable_ld_diagnosis)) .cater_stop("enable_ld_diagnosis must be TRUE or FALSE")
  if(length(ld_diag_loglr)!=1L||!is.finite(ld_diag_loglr)||ld_diag_loglr<0) .cater_stop("ld_diag_loglr must be a finite non-negative scalar")
  if(length(ld_diag_abs_z)!=1L||!is.finite(ld_diag_abs_z)||ld_diag_abs_z<0) .cater_stop("ld_diag_abs_z must be a finite non-negative scalar")
  if(!dir.exists(eqtl_dir)) .cater_stop("eqtl_dir does not exist: %s",eqtl_dir)
  if(!is.null(mvmr_max_r2) && (length(mvmr_max_r2)!=1L || !is.finite(mvmr_max_r2) ||
                               mvmr_max_r2<0 || mvmr_max_r2>1))
    .cater_stop("mvmr_max_r2 must be NULL or a single value in [0,1]")
  if(length(min_cond_F)!=1L||!is.finite(min_cond_F)||min_cond_F<=0)
    .cater_stop("min_cond_F must be a positive finite scalar")
  if(length(mvmr_max_condition)!=1L||!is.finite(mvmr_max_condition)||mvmr_max_condition<=1)
    .cater_stop("mvmr_max_condition must be a finite scalar > 1")
  grn0<-grn;grn<-.cater_standardize_grn(grn)
  ann<-.cater_standardize_annotation(gene_annotation);if(is.null(ann)) ann<-.cater_annotation_from_grn(grn0)
  if(is.null(ann)) .cater_stop("Gene coordinates required via gene_annotation or GRN coordinate columns")
  outcome<-.cater_standardize_sumstats(outcome,label="outcome",require_position=FALSE)
  if(is.null(targets)) targets<-sort(unique(c(grn$TF,grn$Target)))
  targets<-unique(as.character(targets))
  for(d in c(outdir,file.path(outdir,"cojo"),file.path(outdir,"instruments"),file.path(outdir,"diagnostics"),
             file.path(outdir,"mechanism"),file.path(outdir,"pleiotropy"),file.path(outdir,"cojo_mvmr")))
    dir.create(d,recursive=TRUE,showWarnings=FALSE)

  long<-list();summ<-list()
  empty_summary<-function(x,status) data.frame(cell_type=cell_type,target=x,trait=trait,status=status,
    n_parent_tf=0,n_candidate_snp=0,n_ld_diag_removed=0,n_ld_diag_susie=0,
    n_ld_reference_missing=0,n_ld_allele_mismatch=0,n_ld_nonfinite=0,ld_diag_lambda_max=NA_real_,
    n_cojo_signal=0,n_cis_signal=0,n_trans_signal=0,n_ld_aligned_signal=0,
    n_mr_iv=0,n_mr_cis_iv=0,n_mr_trans_iv=0,
    effective_F=NA,trans_information_fraction=NA,n_tf_anchor_tested=0,n_tf_anchor_fdr=0,
    n_sibling_candidate=0,n_sibling_active=0,n_sibling_incomplete=0,n_sibling_qtl_missing=0,
    n_sibling_iv_missing=0,sibling_screen_performed=FALSE,sibling_screen_complete=NA,
    cis_trans_heterogeneity_p=NA,max_tf_locus_weight=NA,leave_one_tf_max_delta=NA,
    network_status=NA_character_,beta_network=NA,se_network=NA,p_network=NA,
    conditional_F_X=NA,conditional_F_converged=NA,mvmr_rank=NA,mvmr_condition=NA,
    mvmr_n_iv=NA_integer_,mvmr_union_requested=NA_integer_,mvmr_union_ld_available=NA_integer_,
    mvmr_union_complete_case=NA_integer_,mvmr_max_r2=NA,mvmr_ld_gate_pass=NA,
    primary_model=NA_character_,primary_beta=NA,primary_se=NA,primary_p=NA,stringsAsFactors=FALSE)

  for(target in targets){
    .cater_msg(verbose,"[%s] starting",target)
    srow<-empty_summary(target,"START")
    parents<-unique(grn$TF[grn$Target==target]);srow$n_parent_tf<-length(parents)
    qtl<-.cater_read_gene(target,eqtl_dir,qtl_n)
    if(is.null(qtl)){srow$status<-"NO_EQTL_FILE";summ[[target]]<-srow;next}
    regions<-.cater_make_regions(target,parents,ann,cis_window,tf_window)
    if(is.null(regions)){srow$status<-"NO_TARGET_ANNOTATION";summ[[target]]<-srow;next}
    cmap<-.cater_candidate_map(qtl,regions);srow$n_candidate_snp<-nrow(cmap)
    if(!nrow(cmap)){srow$status<-"NO_CANDIDATE_SNP";summ[[target]]<-srow;next}
    co<-tryCatch(.cater_run_cojo(qtl,cmap,ld_bfile,manc_cojo_bin,cojo_p,cojo_wind_kb,
      cojo_collinear,cojo_threads,file.path(outdir,"cojo",target),verbose,
      plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,
      ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z),error=function(e)e)
    if(inherits(co,"error")){
      warning(sprintf("[%s] %s",target,conditionMessage(co)))
      srow$status<-if(grepl("SuSiE-RSS|LD diagnosis|PLINK",conditionMessage(co))) "LD_DIAGNOSIS_FAILED" else "COJO_FAILED"
      summ[[target]]<-srow;next
    }
    if(nrow(co$ld_diagnosis)){
      dd<-co$ld_diagnosis
      srow$n_ld_diag_removed<-sum(dd$remove,na.rm=TRUE)
      srow$n_ld_diag_susie<-sum(dd$status=="SUSIE_LD_INCONSISTENT",na.rm=TRUE)
      srow$n_ld_reference_missing<-sum(dd$status=="LD_REFERENCE_MISSING",na.rm=TRUE)
      srow$n_ld_allele_mismatch<-sum(dd$status=="LD_ALLELE_MISMATCH",na.rm=TRUE)
      srow$n_ld_nonfinite<-sum(dd$status=="LD_NONFINITE",na.rm=TRUE)
      if(any(is.finite(dd$lambda))) srow$ld_diag_lambda_max<-max(dd$lambda[is.finite(dd$lambda)])
    }
    if(!nrow(co$selected)){
      srow$status<-if(nrow(co$ld_diagnosis)&&all(co$ld_diagnosis$remove)) "NO_CANDIDATE_AFTER_LD_DIAGNOSIS" else "NO_COJO_SIGNAL"
      summ[[target]]<-srow;next
    }
    co_meta<-cmap[match(co$selected$snp,cmap$snp),,drop=FALSE]
    srow$n_cojo_signal<-nrow(co$selected)
    srow$n_cis_signal<-sum(co_meta$source=="cis",na.rm=TRUE)
    srow$n_trans_signal<-sum(co_meta$source=="trans",na.rm=TRUE)
    sel<-qtl[match(co$selected$snp,qtl$snp),,drop=FALSE];sel<-.cater_align_to_ld(sel,co$selected)
    srow$n_ld_aligned_signal<-nrow(sel)
    if(!nrow(sel)){srow$status<-"NO_LD_ALLELE_MATCH";summ[[target]]<-srow;next}
    meta<-cmap[match(sel$snp,cmap$snp),,drop=FALSE]
    sel$source<-meta$source;sel$parent_tf<-meta$parent_tf;sel$locus_id<-meta$locus_id
    h<-.cater_harmonize(sel,outcome,drop_palindromic)
    srow$n_mr_iv<-nrow(h)
    srow$n_mr_cis_iv<-if(nrow(h)) sum(h$source=="cis") else 0L
    srow$n_mr_trans_iv<-if(nrow(h)) sum(h$source=="trans") else 0L
    if(!nrow(h)){srow$status<-"NO_HARMONIZED_IV";summ[[target]]<-srow;next}
    co$ld<-.cater_subset_ld(co$ld,h$snp)
    utils::write.table(h,file.path(outdir,"instruments",paste0(target,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)

    fits<-list()
    for(m in c("cis","trans","combined")){
      d<-if(m=="combined") h else h[h$source==m,,drop=FALSE]
      R<-if(nrow(d)) .cater_subset_ld(co$ld,d$snp) else matrix(numeric(),0,0)
      f<-.cater_givw(d,R); fits[[m]]<-f
      long[[paste(target,m,sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model=m,f,stringsAsFactors=FALSE)
    }
    srow$effective_F<-fits$combined$effective_F
    srow$trans_information_fraction<-.cater_trans_information_fraction(fits$cis,fits$combined)
    ht<-.cater_cis_trans_het(h,co$ld);srow$cis_trans_heterogeneity_p<-ht["p"]
    ldgn<-.cater_locus_diagnostics(h,co$ld);srow$max_tf_locus_weight<-ldgn$max_weight;srow$leave_one_tf_max_delta<-ldgn$leave_one_max_delta
    if(nrow(ldgn$table)) utils::write.table(ldgn$table,file.path(outdir,"diagnostics",paste0(target,"_loci.tsv")),sep="\t",quote=FALSE,row.names=FALSE)

    anchor<-.cater_tf_anchor(target,h,co$selected,co$ld,eqtl_dir,qtl_n)
    if(nrow(anchor)){
      anchor$q_tf<-NA_real_;aa<-which(anchor$status=="OK"&is.finite(anchor$p_tf))
      if(length(aa)) anchor$q_tf[aa]<-p.adjust(anchor$p_tf[aa],method="BH")
      anchor$evidence<-ifelse(anchor$status!="OK","TF_ANCHOR_NOT_TESTABLE",ifelse(is.finite(anchor$q_tf)&anchor$q_tf<0.05,"TF_ANCHOR_SUPPORTED","TF_ANCHOR_NOT_SUPPORTED"))
      utils::write.table(anchor,file.path(outdir,"mechanism",paste0(target,"_tf_anchor.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
      srow$n_tf_anchor_tested<-sum(anchor$status=="OK")
      srow$n_tf_anchor_fdr<-sum(anchor$status=="OK"&is.finite(anchor$q_tf)&anchor$q_tf<0.05)
    }

    has_trans<-any(h$source=="trans")
    sib<-list(table=data.frame(),active=character(),n_candidate=0L,n_incomplete=0L,
              n_qtl_missing=0L,n_iv_missing=0L,complete=TRUE)
    if(enable_sibling_screen&&has_trans){
      srow$sibling_screen_performed<-TRUE
      sib<-.cater_sibling_screen(target,h,grn,co$selected,co$ld,eqtl_dir,qtl_n,sibling_fdr)
      srow$sibling_screen_complete<-sib$complete
      if(nrow(sib$table)) utils::write.table(sib$table,file.path(outdir,"pleiotropy",paste0(target,"_siblings.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    } else if(has_trans) {
      srow$sibling_screen_complete<-FALSE
    }
    srow$n_sibling_candidate<-sib$n_candidate;srow$n_sibling_active<-length(sib$active)
    srow$n_sibling_incomplete<-sib$n_incomplete;srow$n_sibling_qtl_missing<-sib$n_qtl_missing
    srow$n_sibling_iv_missing<-sib$n_iv_missing

    net<-NULL
    if(enable_mvmr&&length(sib$active)){
      net<-.cater_build_mvmr(target,qtl,co$selected,sib$active,ann,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
        cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,outdir,outcome,drop_palindromic,
        exposure_corr,min_cond_F,mvmr_max_r2,mvmr_max_condition,verbose,
        plink_bin=plink_bin,enable_ld_diagnosis=enable_ld_diagnosis,
        ld_diag_loglr=ld_diag_loglr,ld_diag_abs_z=ld_diag_abs_z)
      srow$network_status<-net$status
      if(identical(net$status,"OK")){
        srow$beta_network<-net$beta[target];srow$se_network<-net$se[target];srow$p_network<-net$p[target]
        srow$conditional_F_X<-net$conditional_F[target]
        srow$conditional_F_converged<-isTRUE(net$conditional_F_converged[target])
        srow$mvmr_rank<-net$rank;srow$mvmr_condition<-net$condition;srow$mvmr_n_iv<-net$n_iv
        srow$mvmr_union_requested<-net$n_union_requested
        srow$mvmr_union_ld_available<-net$n_union_ld_available
        srow$mvmr_union_complete_case<-net$n_union_complete_case
        srow$mvmr_max_r2<-net$max_r2;srow$mvmr_ld_gate_pass<-net$ld_gate_pass
        long[[paste(target,"network",sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model="network",
          n_iv=net$n_iv,beta=net$beta[target],se=net$se[target],p=net$p[target],Q=net$Q,Q_p=net$Q_p,
          information=NA,joint_wald_per_df=NA,effective_F=net$conditional_F[target],precision_information=NA,mean_F=NA,min_F=NA,status="SENSITIVITY_ONLY",
          stringsAsFactors=FALSE)
      }
    }

    if(!is.null(net)&&identical(net$status,"OK"))
      net$primary_eligible<-isTRUE(net$primary_eligible)&&isTRUE(sib$complete)
    dec<-.cater_primary_decision(fits,net=net,target=target,has_trans=has_trans,
      sibling_screen_performed=srow$sibling_screen_performed,
      sibling_screen_complete=if(has_trans) isTRUE(srow$sibling_screen_complete) else TRUE,
      n_active_siblings=length(sib$active),primary_policy=primary_policy,
      sibling_screen_testable=if(has_trans) isTRUE(sib$testable %||% FALSE) else TRUE)
    srow$primary_model<-dec$model;srow$primary_beta<-dec$beta;srow$primary_se<-dec$se;srow$primary_p<-dec$p
    srow$status<-dec$status
    summ[[target]]<-srow
  }

  longdf<-if(length(long)) do.call(rbind,long) else data.frame()
  sumdf<-if(length(summ)) do.call(rbind,summ) else data.frame()
  if(nrow(longdf)){
    longdf$q<-NA_real_
    for(m in unique(longdf$model)){
      ii<-which(longdf$model==m&is.finite(longdf$p))
      if(length(ii)) longdf$q[ii]<-p.adjust(longdf$p[ii],method="BH")
    }
    utils::write.table(longdf,file.path(outdir,"cater_mr_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  }
  if(nrow(sumdf)){
    sumdf$primary_q<-NA_real_;ii<-which(is.finite(sumdf$primary_p))
    if(length(ii)) sumdf$primary_q[ii]<-p.adjust(sumdf$primary_p[ii],method="BH")
    utils::write.table(sumdf,file.path(outdir,"cater_mr_target_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  }
  structure(list(results=longdf,targets=sumdf),class="cater_mr_result")
}
