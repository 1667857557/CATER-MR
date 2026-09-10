# CATER-MR v0.4
# Cis And Trans eQTLs guided by Regulatory networks for drug-target MR
# Target-centric Manc-COJO + GIVW + same-IV pleiotropy screen + triggered local MVMR.

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
    chr = if (is.null(chr)) NA_character_ else sub("^chr","",as.character(x[[chr]]),ignore.case=TRUE),
    pos = if (is.null(pos)) NA_real_ else suppressWarnings(as.numeric(x[[pos]])),
    a1 = toupper(as.character(x[[a1]])),
    a2 = toupper(as.character(x[[a2]])),
    beta = suppressWarnings(as.numeric(x[[b]])),
    se = suppressWarnings(as.numeric(x[[se]])),
    stringsAsFactors = FALSE
  )
  out$p <- if (is.null(p)) 2 * stats::pnorm(-abs(out$beta/out$se)) else suppressWarnings(as.numeric(x[[p]]))
  out$eaf <- if (is.null(eaf)) NA_real_ else suppressWarnings(as.numeric(x[[eaf]]))
  out$n <- if (is.null(n)) {
    if (is.null(n_default)) NA_real_ else rep(as.numeric(n_default), nrow(x))
  } else suppressWarnings(as.numeric(x[[n]]))

  keep <- !is.na(out$snp) & nzchar(out$snp) &
    (!require_position | is.finite(out$pos)) &
    is.finite(out$beta) & is.finite(out$se) & out$se > 0 &
    is.finite(out$p) & out$p >= 0 & out$p <= 1 &
    out$a1 %in% c("A","C","G","T") & out$a2 %in% c("A","C","G","T")
  out <- out[keep,,drop=FALSE]
  dup <- duplicated(out$snp) | duplicated(out$snp, fromLast=TRUE)
  out <- out[!dup,,drop=FALSE]
  rownames(out) <- NULL
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
    data.frame(symbol=as.character(grn_raw[[tf]]),chr=as.character(grn_raw[[tf_chr]]),tss=as.numeric(grn_raw[[tf_tss]])),
    data.frame(symbol=as.character(grn_raw[[tg]]),chr=as.character(grn_raw[[tg_chr]]),tss=as.numeric(grn_raw[[tg_tss]]))
  )
  a$chr <- sub("^chr","",a$chr,ignore.case=TRUE)
  a <- a[!is.na(a$symbol)&nzchar(a$symbol)&is.finite(a$tss),,drop=FALSE]
  a[!duplicated(a$symbol),,drop=FALSE]
}

.cater_standardize_annotation <- function(annotation) {
  if (is.null(annotation)) return(NULL)
  if (!is.data.frame(annotation)) .cater_stop("gene_annotation must be a data.frame")
  sym <- .cater_pick_col(annotation,c("symbol","gene","gene_symbol","SYMBOL"),label="gene annotation symbol")
  chr <- .cater_pick_col(annotation,c("chr","CHR","chrom","chromosome"),label="gene annotation chromosome")
  tss <- .cater_pick_col(annotation,c("tss","TSS","txStart"),label="gene annotation TSS")
  out <- data.frame(symbol=as.character(annotation[[sym]]),
                    chr=sub("^chr","",as.character(annotation[[chr]]),ignore.case=TRUE),
                    tss=as.numeric(annotation[[tss]]),stringsAsFactors=FALSE)
  out <- out[!is.na(out$symbol)&nzchar(out$symbol)&is.finite(out$tss),,drop=FALSE]
  out[!duplicated(out$symbol),,drop=FALSE]
}

.cater_make_regions <- function(target, parents, annotation, cis_window, tf_window) {
  ta <- annotation[annotation$symbol==target,,drop=FALSE]
  if (!nrow(ta)) return(NULL)
  reg <- data.frame(type="cis",gene=target,chr=ta$chr[1],
                    start=max(1,ta$tss[1]-cis_window),end=ta$tss[1]+cis_window,
                    stringsAsFactors=FALSE)
  if (length(parents)) {
    pa <- annotation[match(parents,annotation$symbol),,drop=FALSE]
    pa$gene <- parents
    pa <- pa[!is.na(pa$chr)&is.finite(pa$tss),,drop=FALSE]
    if (nrow(pa)) reg <- rbind(reg,data.frame(type="trans",gene=pa$gene,chr=pa$chr,
      start=pmax(1,pa$tss-tf_window),end=pa$tss+tf_window,stringsAsFactors=FALSE))
  }
  rownames(reg) <- NULL
  reg
}

.cater_candidate_map <- function(qtl, regions) {
  if (is.null(regions)||!nrow(qtl)) return(data.frame())
  cis <- regions[regions$type=="cis",,drop=FALSE]
  is_cis <- qtl$chr==cis$chr[1] & qtl$pos>=cis$start[1] & qtl$pos<=cis$end[1]
  tr <- regions[regions$type=="trans",,drop=FALSE]
  hits <- vector("list",nrow(qtl))
  if (nrow(tr)) for (k in seq_len(nrow(tr))) {
    z <- qtl$chr==tr$chr[k] & qtl$pos>=tr$start[k] & qtl$pos<=tr$end[k] & !is_cis
    if (any(z)) hits[z] <- lapply(hits[z], function(x) c(x,tr$gene[k]))
  }
  is_trans <- lengths(hits)>0L
  keep <- is_cis|is_trans
  if (!any(keep)) return(data.frame())
  parents <- vapply(hits[keep],function(x) if(length(x)) paste(sort(unique(x)),collapse=";") else "",character(1))
  data.frame(snp=qtl$snp[keep],source=ifelse(is_cis[keep],"cis","trans"),
             parent_tf=parents,
             locus_id=ifelse(is_cis[keep],"cis",paste0("TF:",parents)),
             stringsAsFactors=FALSE)
}

.cater_write_cojo_ma <- function(qtl,path) {
  ok <- is.finite(qtl$eaf)&qtl$eaf>0&qtl$eaf<1&is.finite(qtl$n)&qtl$n>0
  if (!all(ok)) .cater_stop("Manc-COJO requires valid EAF/freq and N; %d rows invalid",sum(!ok))
  ma <- data.frame(SNP=qtl$snp,A1=qtl$a1,A2=qtl$a2,freq=qtl$eaf,
                   b=qtl$beta,se=qtl$se,p=qtl$p,N=qtl$n,stringsAsFactors=FALSE)
  utils::write.table(ma,path,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}

.cater_read_manc_ldr <- function(path,selected) {
  selected <- as.character(selected)
  if (!length(selected)) return(matrix(numeric(),0,0))
  if (length(selected)==1L) return(matrix(1,1,1,dimnames=list(selected,selected)))
  if (!file.exists(path)) .cater_stop("Manc-COJO LD output not found: %s",path)
  lines <- trimws(readLines(path,warn=FALSE)); lines <- lines[nzchar(lines)]
  ld <- diag(length(selected)); dimnames(ld) <- list(selected,selected)
  i <- 1L
  while (i<=length(lines)) {
    if (startsWith(lines[i],"#")) {i<-i+1L;next}
    if (!grepl("^SNP(\\s|$)",lines[i])) {i<-i+1L;next}
    header <- strsplit(lines[i],"\\s+")[[1L]]; block <- header[-1L]; i<-i+1L
    rows <- list()
    while(i<=length(lines)&&!startsWith(lines[i],"#")&&!grepl("^SNP(\\s|$)",lines[i])){
      tok <- strsplit(lines[i],"\\s+")[[1L]]
      if(length(tok)>=2L) rows[[length(rows)+1L]]<-tok
      i<-i+1L
    }
    if(!length(rows)||!length(block)) next
    rn<-vapply(rows,`[`,character(1),1L)
    vals<-do.call(rbind,lapply(rows,function(z) as.numeric(z[-1L])))
    if(ncol(vals)!=length(block)) .cater_stop("Malformed Manc-COJO .ldr.cojo block in %s",path)
    rownames(vals)<-rn;colnames(vals)<-block
    common<-intersect(intersect(rn,block),selected)
    if(length(common)) ld[common,common]<-vals[common,common,drop=FALSE]
  }
  ld
}

.cater_read_jma_ref <- function(path) {
  if (!file.exists(path)) return(data.frame())
  x <- utils::read.table(path,header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  snp <- .cater_pick_col(x,"SNP",label="Manc-COJO SNP")
  a1 <- .cater_pick_col(x,"A1",label="Manc-COJO A1")
  a2 <- .cater_pick_col(x,"A2",label="Manc-COJO A2")
  data.frame(snp=as.character(x[[snp]]),ld_a1=toupper(as.character(x[[a1]])),
             ld_a2=toupper(as.character(x[[a2]])),stringsAsFactors=FALSE)
}

.cater_run_cojo <- function(qtl,candidate_map,ld_bfile,manc_cojo_bin,cojo_p,
                            cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose) {
  if(length(ld_bfile)!=1L||!nzchar(ld_bfile)) .cater_stop("Current CATER-MR expects one LD cohort")
  if(!all(file.exists(paste0(ld_bfile,c(".bed",".bim",".fam")))))
    .cater_stop("ld_bfile must be a PLINK bed/bim/fam prefix: %s",ld_bfile)
  exe<-Sys.which(manc_cojo_bin); if(!nzchar(exe)) .cater_stop("Cannot find Manc-COJO executable '%s'",manc_cojo_bin)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
  ma<-paste0(prefix,".sumstat"); extract<-paste0(prefix,".candidate.snplist")
  .cater_write_cojo_ma(qtl,ma); writeLines(unique(candidate_map$snp),extract)
  sp<-paste0(prefix,".select")
  args<-c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",extract,"--cojo-slct",
          "--cojo-p",format(cojo_p,scientific=TRUE),"--cojo-wind",as.character(as.integer(cojo_wind_kb)),
          "--cojo-collinear",as.character(cojo_collinear),"--thread-num",as.character(as.integer(cojo_threads)),
          "--out",sp)
  .cater_msg(verbose,"Manc-COJO selection: %s %s",exe,paste(args,collapse=" "))
  st<-system2(exe,args=args,stdout=paste0(sp,".stdout"),stderr=paste0(sp,".stderr"))
  if(!identical(st,0L)) .cater_stop("Manc-COJO selection failed for %s; see stderr",prefix)
  jma<-paste0(sp,".jma.cojo"); selected<-.cater_read_jma_ref(jma)
  if(!nrow(selected)) return(list(selected=selected,ld=matrix(numeric(),0,0)))
  if(nrow(selected)==1L) return(list(selected=selected,ld=matrix(1,1,1,dimnames=list(selected$snp,selected$snp))))
  jp<-paste0(prefix,".joint")
  args2<-c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",jma,"2","header","--cojo-joint",
           "--thread-num",as.character(as.integer(cojo_threads)),"--output-all","--out",jp)
  st<-system2(exe,args=args2,stdout=paste0(jp,".stdout"),stderr=paste0(jp,".stderr"))
  if(!identical(st,0L)) .cater_stop("Manc-COJO joint/LD failed for %s",prefix)
  ld<-.cater_read_manc_ldr(paste0(jp,".ldr.cojo"),selected$snp)
  list(selected=selected,ld=ld)
}

.cater_joint_ld <- function(qtl,snps,ld_bfile,manc_cojo_bin,cojo_threads,prefix,verbose) {
  snps<-unique(as.character(snps)); snps<-snps[snps%in%qtl$snp]
  if(!length(snps)) return(list(selected=data.frame(),ld=matrix(numeric(),0,0)))
  if(length(snps)==1L) {
    q<-qtl[match(snps,qtl$snp),,drop=FALSE]
    ref<-data.frame(snp=q$snp,ld_a1=q$a1,ld_a2=q$a2,stringsAsFactors=FALSE)
    return(list(selected=ref,ld=matrix(1,1,1,dimnames=list(snps,snps))))
  }
  exe<-Sys.which(manc_cojo_bin); if(!nzchar(exe)) .cater_stop("Cannot find Manc-COJO executable '%s'",manc_cojo_bin)
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
  ma<-paste0(prefix,".sumstat"); ex<-paste0(prefix,".snplist")
  .cater_write_cojo_ma(qtl,ma); writeLines(snps,ex)
  args<-c("--bfile",ld_bfile,"--cojo-file",ma,"--extract",ex,"--cojo-joint",
          "--thread-num",as.character(as.integer(cojo_threads)),"--output-all","--out",prefix)
  .cater_msg(verbose,"Manc-COJO MVMR LD: %s %s",exe,paste(args,collapse=" "))
  st<-system2(exe,args=args,stdout=paste0(prefix,".stdout"),stderr=paste0(prefix,".stderr"))
  if(!identical(st,0L)) .cater_stop("Manc-COJO joint LD failed for %s",prefix)
  ref<-.cater_read_jma_ref(paste0(prefix,".jma.cojo"))
  if(!nrow(ref)) return(list(selected=ref,ld=matrix(numeric(),0,0)))
  ld<-.cater_read_manc_ldr(paste0(prefix,".ldr.cojo"),ref$snp)
  list(selected=ref,ld=ld)
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
  if(!nrow(dat)) return(list(information=NA_real_,effective_F=NA_real_,mean_F=NA_real_,min_F=NA_real_))
  f<-(dat$bx/dat$bx_se)^2
  Dx<-diag(dat$bx_se,nrow=nrow(dat)); Sx<-Dx%*%ld%*%Dx; inv<-.cater_inv(Sx)
  info<-if(is.null(inv)) NA_real_ else as.numeric(crossprod(dat$bx,inv%*%dat$bx))
  list(information=info,effective_F=if(is.finite(info)) info/nrow(dat) else NA_real_,
       mean_F=mean(f),min_F=min(f))
}

.cater_givw <- function(dat,ld) {
  st<-.cater_strength(dat,ld); n<-nrow(dat)
  if(!n) return(data.frame(n_iv=0L,beta=NA,se=NA,p=NA,Q=NA,Q_p=NA,
    information=NA,effective_F=NA,mean_F=NA,min_F=NA,status="NO_IV"))
  if(n==1L) {
    b<-dat$by/dat$bx;s<-abs(dat$by_se/dat$bx)
    return(data.frame(n_iv=1L,beta=b,se=s,p=2*stats::pnorm(-abs(b/s)),Q=NA,Q_p=NA,
      information=st$information,effective_F=st$effective_F,mean_F=st$mean_F,min_F=st$min_F,status="OK"))
  }
  Dy<-diag(dat$by_se,nrow=n); Oy<-Dy%*%ld%*%Dy; inv<-.cater_inv(Oy)
  if(is.null(inv)) return(data.frame(n_iv=n,beta=NA,se=NA,p=NA,Q=NA,Q_p=NA,
    information=st$information,effective_F=st$effective_F,mean_F=st$mean_F,min_F=st$min_F,status="LD_SINGULAR"))
  den<-as.numeric(crossprod(dat$bx,inv%*%dat$bx))
  if(!is.finite(den)||den<=0) return(data.frame(n_iv=n,beta=NA,se=NA,p=NA,Q=NA,Q_p=NA,
    information=st$information,effective_F=st$effective_F,mean_F=st$mean_F,min_F=st$min_F,status="INVALID_INFORMATION"))
  b<-as.numeric(crossprod(dat$bx,inv%*%dat$by))/den;s<-sqrt(1/den)
  r<-dat$by-b*dat$bx;q<-as.numeric(crossprod(r,inv%*%r))
  data.frame(n_iv=n,beta=b,se=s,p=2*stats::pnorm(-abs(b/s)),Q=q,
    Q_p=stats::pchisq(q,n-1,lower.tail=FALSE),information=st$information,
    effective_F=st$effective_F,mean_F=st$mean_F,min_F=st$min_F,status="OK")
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

.cater_locus_diagnostics <- function(h,ld) {
  tr<-h[h$source=="trans",,drop=FALSE]
  if(!nrow(tr)) return(list(table=data.frame(),max_weight=NA_real_,leave_one_max_delta=NA_real_))
  full<-.cater_givw(h,ld); total_info<-.cater_strength(h,ld)$information
  keys<-unique(tr$locus_id); tabs<-list(); deltas<-numeric()
  for(k in keys){
    d<-h[h$locus_id==k,,drop=FALSE];R<-.cater_subset_ld(ld,d$snp)
    keep<-h$locus_id!=k
    fminus<-if(any(keep)) .cater_givw(h[keep,,drop=FALSE],.cater_subset_ld(ld,h$snp[keep])) else NULL
    info_minus<-if(any(keep)) .cater_strength(h[keep,,drop=FALSE],.cater_subset_ld(ld,h$snp[keep]))$information else 0
    delta<-if(!is.null(fminus)&&full$status=="OK"&&fminus$status=="OK") abs(full$beta-fminus$beta) else NA_real_
    deltas<-c(deltas,delta)
    iw<-if(is.finite(total_info)&&total_info>0&&is.finite(info_minus)) (total_info-info_minus)/total_info else NA_real_
    if(is.finite(iw)) iw<-max(0,min(1,iw))
    tabs[[k]]<-data.frame(locus_id=k,n_iv=nrow(d),incremental_information=if(is.finite(total_info)&&is.finite(info_minus)) total_info-info_minus else NA_real_,
      information_weight=iw,leave_one_delta=delta,stringsAsFactors=FALSE)
  }
  z<-do.call(rbind,tabs)
  list(table=z,max_weight=if(nrow(z)) max(z$information_weight,na.rm=TRUE) else NA_real_,
       leave_one_max_delta=if(any(is.finite(deltas))) max(deltas,na.rm=TRUE) else NA_real_)
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

.cater_sibling_screen <- function(target,h,grn,cojo_ref,ld,eqtl_dir,qtl_n,sibling_fdr) {
  tr<-h[h$source=="trans",,drop=FALSE]
  if(!nrow(tr)) return(list(table=data.frame(),active=character(),n_candidate=0L,n_missing=0L))
  children_x<-unique(grn$Target[grn$TF==target])
  tf_to_snps<-list()
  for(i in seq_len(nrow(tr))) for(tf in strsplit(tr$parent_tf[i],";",fixed=TRUE)[[1L]])
    tf_to_snps[[tf]]<-unique(c(tf_to_snps[[tf]],tr$snp[i]))
  sib_parents<-list()
  for(tf in names(tf_to_snps)){
    sibs<-setdiff(unique(grn$Target[grn$TF==tf]),c(target,children_x))
    for(z in sibs) sib_parents[[z]]<-unique(c(sib_parents[[z]],tf))
  }
  if(!length(sib_parents)) return(list(table=data.frame(),active=character(),n_candidate=0L,n_missing=0L))
  out<-list();missing<-0L
  for(z in names(sib_parents)){
    tfs<-sib_parents[[z]]
    snps<-unique(unlist(tf_to_snps[tfs],use.names=FALSE))
    d<-.cater_exact_effects(z,snps,eqtl_dir,qtl_n,cojo_ref)
    if(is.null(d)) {
      missing<-missing+1L
      out[[z]]<-data.frame(sibling=z,parent_tf=paste(sort(tfs),collapse=";"),n_iv=0,Q=NA,df=NA,p=NA,status="QTL_NOT_AVAILABLE")
      next
    }
    d<-d[match(intersect(snps,d$snp),d$snp),,drop=FALSE]
    if(!nrow(d)){
      out[[z]]<-data.frame(sibling=z,parent_tf=paste(sort(tfs),collapse=";"),n_iv=0,Q=NA,df=NA,p=NA,status="SNP_NOT_AVAILABLE")
      next
    }
    R<-.cater_subset_ld(ld,d$snp); om<-.cater_omnibus(d$beta,d$se,R)
    out[[z]]<-data.frame(sibling=z,parent_tf=paste(sort(tfs),collapse=";"),n_iv=nrow(d),
                         Q=om["Q"],df=om["df"],p=om["p"],status=if(is.finite(om["p"]))"OK" else "TEST_FAILED")
  }
  tab<-do.call(rbind,out);rownames(tab)<-NULL
  tab$q<-NA_real_;ii<-which(is.finite(tab$p));if(length(ii)) tab$q[ii]<-p.adjust(tab$p[ii],method="BH")
  tab$active<-is.finite(tab$q)&tab$q<sibling_fdr
  list(table=tab,active=tab$sibling[tab$active],n_candidate=nrow(tab),n_missing=missing)
}

.cater_conditional_f <- function(B,SE,exposure_names,exposure_corr=NULL) {
  p<-ncol(B);m<-nrow(B)
  if(p<2L||m<1L) return(setNames(rep(NA_real_,p),exposure_names))
  if(!is.null(exposure_corr)){
    if(is.null(rownames(exposure_corr))||!all(exposure_names%in%rownames(exposure_corr)))
      .cater_stop("exposure_corr must have named rows/columns covering MVMR exposures")
    C<-as.matrix(exposure_corr[exposure_names,exposure_names,drop=FALSE])
  } else C<-diag(p)
  ans<-rep(NA_real_,p)
  for(i in seq_len(p)){
    other<-setdiff(seq_len(p),i)
    fit<-tryCatch(stats::lm.fit(x=B[,other,drop=FALSE],y=B[,i]),error=function(e)NULL)
    if(is.null(fit)||any(!is.finite(fit$coefficients))) next
    delta<-fit$coefficients
    resid<-B[,i]-as.numeric(B[,other,drop=FALSE]%*%delta)
    v<-numeric(p);v[i]<--1;v[other]<-delta
    sig2<-numeric(m)
    for(j in seq_len(m)){
      covj<-C*outer(SE[j,],SE[j,])
      sig2[j]<-as.numeric(t(v)%*%covj%*%v)
    }
    if(all(is.finite(sig2)&sig2>0)) ans[i]<-mean(resid^2/sig2)
  }
  setNames(ans,exposure_names)
}

.cater_mvmr_fit <- function(B,seB,by,seY,ld,exposure_names,exposure_corr=NULL) {
  m<-nrow(B);p<-ncol(B)
  if(m<=p) return(list(status="MVMR_UNDERIDENTIFIED"))
  Dy<-diag(seY,nrow=m);Oy<-Dy%*%ld%*%Dy;invY<-.cater_inv(Oy)
  if(is.null(invY)) return(list(status="MVMR_LD_SINGULAR"))
  info<-t(B)%*%invY%*%B
  if(qr(info)$rank<p) return(list(status="MVMR_RANK_DEFICIENT",rank=qr(info)$rank,condition=Inf))
  iv<-.cater_inv(info);if(is.null(iv)) return(list(status="MVMR_RANK_DEFICIENT",rank=qr(info)$rank,condition=Inf))
  th<-as.numeric(iv%*%t(B)%*%invY%*%by); names(th)<-exposure_names
  ses<-sqrt(diag(iv)); names(ses)<-exposure_names
  ps<-2*stats::pnorm(-abs(th/ses));names(ps)<-exposure_names
  r<-by-as.numeric(B%*%th);Q<-as.numeric(t(r)%*%invY%*%r)
  cf<-.cater_conditional_f(B,seB,exposure_names,exposure_corr)
  list(status="OK",beta=th,se=ses,p=ps,Q=Q,Q_p=stats::pchisq(Q,m-p,lower.tail=FALSE),
       rank=qr(info)$rank,condition=kappa(info),conditional_F=cf,
       covariance_assumption=if(is.null(exposure_corr))"zero_within_SNP_exposure_covariance" else "user_exposure_corr")
}

.cater_run_sibling_cis <- function(gene,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
                                   cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose) {
  q<-.cater_read_gene(gene,eqtl_dir,qtl_n);if(is.null(q)) return(NULL)
  reg<-.cater_make_regions(gene,character(),annotation,cis_window,cis_window)
  if(is.null(reg)) return(NULL)
  cmap<-.cater_candidate_map(q,reg);if(!nrow(cmap)) return(list(qtl=q,selected=data.frame()))
  co<-.cater_run_cojo(q,cmap,ld_bfile,manc_cojo_bin,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,prefix,verbose)
  list(qtl=q,selected=co$selected)
}

.cater_build_mvmr <- function(target,target_qtl,target_sel,active_siblings,annotation,eqtl_dir,qtl_n,
                              ld_bfile,manc_cojo_bin,cis_window,cojo_p,cojo_wind_kb,cojo_collinear,
                              cojo_threads,outdir,outcome,drop_palindromic,exposure_corr,min_cond_F,
                              mvmr_max_r2,mvmr_max_condition,verbose) {
  sibco<-list();union<-target_sel$snp
  for(z in active_siblings){
    pr<-file.path(outdir,"cojo_mvmr",target,z)
    x<-tryCatch(.cater_run_sibling_cis(z,annotation,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
      cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,pr,verbose),error=function(e)NULL)
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
  off<-R;diag(off)<-0
  fit$max_r2<-if(length(off)) max(off^2,na.rm=TRUE) else 0
  fit$primary_eligible<-is.finite(fit$conditional_F[target])&&fit$conditional_F[target]>=min_cond_F&&
    !is.null(exposure_corr)&&is.finite(fit$condition)&&fit$condition<=mvmr_max_condition&&
    is.finite(fit$max_r2)&&fit$max_r2<=mvmr_max_r2
  fit
}

#' Run CATER-MR
#'
#' Primary inputs: cell-type-specific one-hop GRN, directory of SYMBOL.txt.gz
#' full-summary eQTL files, and an in-memory outcome summary object.
cater_mr <- function(grn,eqtl_dir,outcome,gene_annotation=NULL,ld_bfile,
                     manc_cojo_bin="manc_cojo",targets=NULL,cell_type=NA_character_,trait=NA_character_,
                     cis_window=1e6,tf_window=cis_window,cojo_p=5e-8,cojo_wind_kb=10000L,
                     cojo_collinear=0.9,cojo_threads=1L,qtl_n=NULL,
                     sibling_fdr=0.05,enable_sibling_screen=TRUE,enable_mvmr=TRUE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=0.01,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE) {
  if(!dir.exists(eqtl_dir)) .cater_stop("eqtl_dir does not exist: %s",eqtl_dir)
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
    n_parent_tf=0,n_candidate_snp=0,n_cojo_signal=0,n_cis_signal=0,n_trans_signal=0,
    effective_F=NA,trans_information_fraction=NA,n_tf_anchor_tested=0,n_tf_anchor_fdr=0,
    n_sibling_candidate=0,n_sibling_active=0,n_sibling_qtl_missing=0,
    cis_trans_heterogeneity_p=NA,max_tf_locus_weight=NA,leave_one_tf_max_delta=NA,
    beta_network=NA,se_network=NA,p_network=NA,conditional_F_X=NA,mvmr_rank=NA,mvmr_condition=NA,
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
      cojo_collinear,cojo_threads,file.path(outdir,"cojo",target),verbose),error=function(e)e)
    if(inherits(co,"error")){warning(sprintf("[%s] %s",target,conditionMessage(co)));srow$status<-"COJO_FAILED";summ[[target]]<-srow;next}
    if(!nrow(co$selected)){srow$status<-"NO_COJO_SIGNAL";summ[[target]]<-srow;next}
    sel<-qtl[match(co$selected$snp,qtl$snp),,drop=FALSE];sel<-.cater_align_to_ld(sel,co$selected)
    if(!nrow(sel)){srow$status<-"NO_LD_ALLELE_MATCH";summ[[target]]<-srow;next}
    meta<-cmap[match(sel$snp,cmap$snp),,drop=FALSE]
    sel$source<-meta$source;sel$parent_tf<-meta$parent_tf;sel$locus_id<-meta$locus_id
    h<-.cater_harmonize(sel,outcome,drop_palindromic)
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
    srow$n_cojo_signal<-nrow(h);srow$n_cis_signal<-sum(h$source=="cis");srow$n_trans_signal<-sum(h$source=="trans")
    srow$effective_F<-fits$combined$effective_F
    Ic<-fits$cis$information;Ia<-fits$combined$information
    srow$trans_information_fraction<-if(is.finite(Ia)&&Ia>0&&is.finite(Ic)) max(0,min(1,(Ia-Ic)/Ia)) else NA_real_
    ht<-.cater_cis_trans_het(h,co$ld);srow$cis_trans_heterogeneity_p<-ht["p"]
    ldgn<-.cater_locus_diagnostics(h,co$ld);srow$max_tf_locus_weight<-ldgn$max_weight;srow$leave_one_tf_max_delta<-ldgn$leave_one_max_delta
    if(nrow(ldgn$table)) utils::write.table(ldgn$table,file.path(outdir,"diagnostics",paste0(target,"_loci.tsv")),sep="\t",quote=FALSE,row.names=FALSE)

    anchor<-.cater_tf_anchor(target,h,co$selected,co$ld,eqtl_dir,qtl_n)
    if(nrow(anchor)){
      anchor$q_tf<-NA_real_;aa<-which(anchor$status=="OK"&is.finite(anchor$p_tf))
      if(length(aa)) anchor$q_tf[aa]<-p.adjust(anchor$p_tf[aa],method="BH")
      utils::write.table(anchor,file.path(outdir,"mechanism",paste0(target,"_tf_anchor.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
      srow$n_tf_anchor_tested<-sum(anchor$status=="OK")
      srow$n_tf_anchor_fdr<-sum(anchor$status=="OK"&is.finite(anchor$q_tf)&anchor$q_tf<0.05)
    }

    sib<-list(table=data.frame(),active=character(),n_candidate=0L,n_missing=0L)
    if(enable_sibling_screen&&any(h$source=="trans")){
      sib<-.cater_sibling_screen(target,h,grn,co$selected,co$ld,eqtl_dir,qtl_n,sibling_fdr)
      if(nrow(sib$table)) utils::write.table(sib$table,file.path(outdir,"pleiotropy",paste0(target,"_siblings.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    }
    srow$n_sibling_candidate<-sib$n_candidate;srow$n_sibling_active<-length(sib$active);srow$n_sibling_qtl_missing<-sib$n_missing

    net<-NULL
    if(enable_mvmr&&length(sib$active)){
      net<-.cater_build_mvmr(target,qtl,co$selected,sib$active,ann,eqtl_dir,qtl_n,ld_bfile,manc_cojo_bin,
        cis_window,cojo_p,cojo_wind_kb,cojo_collinear,cojo_threads,outdir,outcome,drop_palindromic,
        exposure_corr,min_cond_F,mvmr_max_r2,mvmr_max_condition,verbose)
      if(identical(net$status,"OK")){
        srow$beta_network<-net$beta[target];srow$se_network<-net$se[target];srow$p_network<-net$p[target]
        srow$conditional_F_X<-net$conditional_F[target];srow$mvmr_rank<-net$rank;srow$mvmr_condition<-net$condition
        long[[paste(target,"network",sep=":")]]<-data.frame(cell_type=cell_type,target=target,trait=trait,model="network",
          n_iv=net$n_iv,beta=net$beta[target],se=net$se[target],p=net$p[target],Q=net$Q,Q_p=net$Q_p,
          information=NA,effective_F=net$conditional_F[target],mean_F=NA,min_F=NA,status=if(net$primary_eligible)"OK" else "SENSITIVITY_ONLY",
          stringsAsFactors=FALSE)
      }
    }

    if(length(sib$active)){
      if(!is.null(net)&&identical(net$status,"OK")&&isTRUE(net$primary_eligible)){
        srow$primary_model<-"network";srow$primary_beta<-srow$beta_network;srow$primary_se<-srow$se_network;srow$primary_p<-srow$p_network
        srow$status<-"OK_NETWORK_ADJUSTED"
      } else if(fits$cis$status=="OK"){
        srow$primary_model<-"cis";srow$primary_beta<-fits$cis$beta;srow$primary_se<-fits$cis$se;srow$primary_p<-fits$cis$p
        srow$status<-"TRANS_PLEIOTROPY_UNRESOLVED"
      } else {
        srow$primary_model<-"combined_sensitivity";srow$primary_beta<-fits$combined$beta;srow$primary_se<-fits$combined$se;srow$primary_p<-fits$combined$p
        srow$status<-"TRANS_PLEIOTROPY_UNRESOLVED_NO_CIS"
      }
    } else if(fits$combined$status=="OK"){
      srow$primary_model<-if(any(h$source=="trans"))"combined" else "cis"
      srow$primary_beta<-fits$combined$beta;srow$primary_se<-fits$combined$se;srow$primary_p<-fits$combined$p
      srow$status<-if(any(h$source=="trans"))"OK_CATER" else "OK_CIS_ONLY"
      if(sib$n_missing>0L) srow$status<-"SIBLING_QTL_INCOMPLETE"
    } else {
      srow$status<-"MR_FAILED"
    }
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
