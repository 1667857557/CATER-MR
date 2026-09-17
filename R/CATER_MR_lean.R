# CATER-MR Lean v0.9 override layer
#
# This file intentionally keeps the audited v0.8 implementation as a compatibility
# substrate, but changes the default causal workflow to:
#   cis primary MR + outcome-independent GRN-mediated trans qualification +
#   signed-LD augmented MR.
# Advanced sibling/MVMR analyses remain available as sensitivity analyses only.

.CATER_VERSION <- "0.9.0"

.cater_target_candidates_v08 <- .cater_target_candidates
.cater_primary_decision_v08 <- .cater_primary_decision
.cater_mr_v08 <- cater_mr

.CATER_LEAN_STATE <- new.env(parent=emptyenv())
.CATER_LEAN_STATE$active <- FALSE
.CATER_LEAN_STATE$config <- NULL
.CATER_LEAN_STATE$diagnostics <- list()

.cater_preld_tf_anchor <- function(tf, snp, target_row, eqtl_dir, qtl_n) {
  d <- .cater_read_gene(tf, eqtl_dir, qtl_n)
  if (is.null(d)) {
    return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=NA_real_,p_tf=NA_real_,status="QTL_NOT_AVAILABLE",stringsAsFactors=FALSE))
  }
  d <- d[d$snp == snp,,drop=FALSE]
  if (!nrow(d)) {
    return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=NA_real_,p_tf=NA_real_,status="SNP_NOT_AVAILABLE",stringsAsFactors=FALSE))
  }
  same <- d$a1[1] == target_row$a1[1] && d$a2[1] == target_row$a2[1]
  swap <- d$a1[1] == target_row$a2[1] && d$a2[1] == target_row$a1[1]
  if (!(same || swap)) {
    return(data.frame(tf=tf,beta_tf=NA_real_,se_tf=d$se[1],p_tf=d$p[1],status="ALLELE_MISMATCH",stringsAsFactors=FALSE))
  }
  beta <- if (swap) -d$beta[1] else d$beta[1]
  p <- d$p[1]
  if (!is.finite(p)) p <- 2*stats::pnorm(-abs(beta/d$se[1]))
  data.frame(tf=tf,beta_tf=beta,se_tf=d$se[1],p_tf=p,status="OK",stringsAsFactors=FALSE)
}

.cater_qualify_trans_candidates <- function(target, qtl, map, eqtl_dir, qtl_n,
                                             trans_hits=NULL, tf_anchor_p=5e-8,
                                             max_reported_trans_targets=1L,
                                             trans_set=c("core","extended")) {
  trans_set <- match.arg(trans_set)
  if (length(tf_anchor_p)!=1L || !is.finite(tf_anchor_p) || tf_anchor_p<=0 || tf_anchor_p>=1)
    .cater_stop("tf_anchor_p must be in (0,1)")
  if (length(max_reported_trans_targets)!=1L || !is.finite(max_reported_trans_targets) ||
      max_reported_trans_targets < 1 || max_reported_trans_targets != as.integer(max_reported_trans_targets))
    .cater_stop("max_reported_trans_targets must be a positive integer")

  if (!nrow(qtl) || !nrow(map)) {
    return(list(qtl=qtl,map=map,table=data.frame(),anchors=data.frame(),
                counts=list(n_trans_raw=0L,n_trans_anchor_pass=0L,n_trans_unique_parent=0L,
                            n_trans_core=0L,n_trans_extended=0L,n_parent_tf_core=0L)))
  }
  ti <- which(map$source == "trans")
  if (!length(ti)) {
    return(list(qtl=qtl,map=map,table=data.frame(),anchors=data.frame(),
                counts=list(n_trans_raw=0L,n_trans_anchor_pass=0L,n_trans_unique_parent=0L,
                            n_trans_core=0L,n_trans_extended=0L,n_parent_tf_core=0L)))
  }

  rows <- vector("list",length(ti)); anchor_rows <- list()
  for (ii in seq_along(ti)) {
    mi <- ti[ii]; s <- map$snp[mi]
    qi <- match(s,qtl$snp)
    tr <- qtl[qi,,drop=FALSE]
    parents <- .cater_parent_tokens(map$parent_tf[mi])[[1L]]
    aa <- if (length(parents)) do.call(rbind,lapply(parents,function(tf)
      .cater_preld_tf_anchor(tf,s,tr,eqtl_dir,qtl_n))) else data.frame()
    if (nrow(aa)) {
      aa$target <- target; aa$snp <- s
      anchor_rows[[length(anchor_rows)+1L]] <- aa
      pass <- aa$status == "OK" & is.finite(aa$p_tf) & aa$p_tf < tf_anchor_p
      qtf <- unique(aa$tf[pass])
    } else qtf <- character()
    nqual <- length(qtf)
    nrep <- NA_integer_
    if (!is.null(trans_hits) && nrow(trans_hits))
      nrep <- length(unique(trans_hits$gene[trans_hits$snp == s]))
    specificity <- is.na(nrep) || nrep <= max_reported_trans_targets
    extended <- nqual == 1L
    core <- extended && specificity
    reason <- if (nqual == 0L) "TRANS_NO_TF_CIS_ANCHOR" else if (nqual > 1L) "TRANS_AMBIGUOUS_PARENT" else if (!specificity) "TRANS_REPORTED_HOTSPOT" else "TRANS_CORE"
    chosen <- if (nqual == 1L) qtf[[1L]] else NA_character_
    aone <- if (nqual == 1L) aa[match(chosen,aa$tf),,drop=FALSE] else NULL
    rows[[ii]] <- data.frame(target=target,snp=s,raw_parent_tf=map$parent_tf[mi],
      qualified_parent=chosen,n_qualifying_parent=nqual,
      beta_tf=if(is.null(aone)) NA_real_ else aone$beta_tf[1],
      se_tf=if(is.null(aone)) NA_real_ else aone$se_tf[1],
      p_tf=if(is.null(aone)) NA_real_ else aone$p_tf[1],
      n_reported_trans_targets=nrep,anchor_pass=nqual>=1L,
      unique_parent=extended,specificity_pass=specificity,
      eligible_extended=extended,eligible_core=core,reason=reason,
      stringsAsFactors=FALSE)
  }
  tab <- do.call(rbind,rows); rownames(tab) <- NULL
  anchors <- if(length(anchor_rows)) do.call(rbind,anchor_rows) else data.frame()

  use <- if(trans_set=="core") tab$eligible_core else tab$eligible_extended
  keep_snps <- tab$snp[use]
  cis_idx <- which(map$source == "cis")
  trans_idx <- ti[match(keep_snps,map$snp[ti],nomatch=0L)]
  trans_idx <- trans_idx[trans_idx>0L]
  keep_idx <- c(cis_idx,trans_idx)
  qout <- qtl[match(map$snp[keep_idx],qtl$snp),,drop=FALSE]
  mout <- map[keep_idx,,drop=FALSE]
  if (length(trans_idx)) {
    dd <- tab[match(mout$snp[mout$source=="trans"],tab$snp),,drop=FALSE]
    mout$parent_tf[mout$source=="trans"] <- dd$qualified_parent
    mout$provenance[mout$source=="trans"] <- if(trans_set=="core") "REPORTED_TRANS_CORE" else "REPORTED_TRANS_EXTENDED"
    qout$provenance[qout$snp %in% mout$snp[mout$source=="trans"]] <- if(trans_set=="core") "REPORTED_TRANS_CORE" else "REPORTED_TRANS_EXTENDED"
  }

  counts <- list(
    n_trans_raw=nrow(tab),
    n_trans_anchor_pass=sum(tab$anchor_pass),
    n_trans_unique_parent=sum(tab$unique_parent),
    n_trans_core=sum(tab$eligible_core),
    n_trans_extended=sum(tab$eligible_extended),
    n_parent_tf_core=length(unique(tab$qualified_parent[tab$eligible_core & !is.na(tab$qualified_parent)]))
  )
  list(qtl=qout,map=mout,table=tab,anchors=anchors,counts=counts)
}

# During an ordinary direct internal call, preserve v0.8 behavior for regression
# compatibility. The public cater_mr() wrapper activates Lean qualification before
# LD selection and before outcome harmonization.
.cater_target_candidates <- function(target,parents,annotation,eqtl_dir,qtl_n,trans_hits,
                                     instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc) {
  raw <- .cater_target_candidates_v08(target,parents,annotation,eqtl_dir,qtl_n,trans_hits,
                                      instrument_p,cis_window,tf_window,drop_palindromic,exclude_mhc)
  if (!isTRUE(.CATER_LEAN_STATE$active) || !nrow(raw$qtl)) return(raw)
  cfg <- .CATER_LEAN_STATE$config
  z <- .cater_qualify_trans_candidates(target,raw$qtl,raw$map,eqtl_dir,qtl_n,trans_hits,
                                        cfg$tf_anchor_p,cfg$max_reported_trans_targets,cfg$trans_set)
  .CATER_LEAN_STATE$diagnostics[[target]] <- z
  list(qtl=z$qtl,map=z$map,regions=raw$regions,qualification=z$table)
}

# Primary inference is intentionally conservative in Lean v0.9: cis is primary when
# available. Trans and augmented estimates are evidence/precision augmentations.
.cater_primary_decision <- function(fits, ...) {
  if (!isTRUE(.CATER_LEAN_STATE$active)) return(.cater_primary_decision_v08(fits,...))
  if (!is.null(fits$cis) && identical(fits$cis$status,"OK")) {
    return(list(model="cis",beta=fits$cis$beta,se=fits$cis$se,p=fits$cis$p,
                status="OK_CIS_PRIMARY",evidence_status="CIS_ANCHOR_PRIMARY"))
  }
  list(model=NA_character_,beta=NA_real_,se=NA_real_,p=NA_real_,
       status="NO_CIS_PRIMARY",evidence_status="TRANS_SENSITIVITY_ONLY")
}

.cater_lean_fit_row <- function(results,target,model) {
  d <- results[results$target==target & results$model==model,,drop=FALSE]
  if(!nrow(d)) return(NULL)
  d[1,,drop=FALSE]
}

.cater_lean_postprocess <- function(res,outdir,tf_anchor_p,max_reported_trans_targets,trans_set) {
  if (!is.list(res)) return(res)
  if (is.data.frame(res$results) && nrow(res$results)) {
    if (!any(res$results$model=="augmented") && any(res$results$model=="combined")) {
      aug <- res$results[res$results$model=="combined",,drop=FALSE]
      aug$model <- "augmented"
      res$results <- rbind(res$results,aug)
    }
    res$results$q <- NA_real_
    for (m in unique(res$results$model)) {
      ii <- which(res$results$model==m & is.finite(res$results$p))
      if(length(ii)) res$results$q[ii] <- stats::p.adjust(res$results$p[ii],method="BH")
    }
  }
  if (is.data.frame(res$targets) && nrow(res$targets)) {
    add_num <- c("n_trans_raw","n_trans_anchor_pass","n_trans_unique_parent","n_trans_core","n_trans_extended","n_parent_tf_core")
    for(nm in add_num) res$targets[[nm]] <- 0L
    res$targets$trans_set_used <- trans_set
    res$targets$tf_anchor_p <- tf_anchor_p
    res$targets$max_reported_trans_targets <- as.integer(max_reported_trans_targets)
    res$targets$cis_effective_F <- NA_real_
    res$targets$trans_effective_F <- NA_real_
    res$targets$augmented_effective_F <- NA_real_
    res$targets$primary_effective_F <- NA_real_
    res$targets$augmented_se_reduction <- NA_real_
    res$targets$augmented_precision_fraction <- NA_real_
    for(i in seq_len(nrow(res$targets))) {
      target <- res$targets$target[i]
      z <- .CATER_LEAN_STATE$diagnostics[[target]]
      if(!is.null(z)) {
        for(nm in add_num) res$targets[[nm]][i] <- as.integer(z$counts[[nm]])
        if(nrow(z$table)) {
          fn <- file.path(outdir,"mechanism",paste0(target,"_trans_qualification.tsv"))
          utils::write.table(z$table,fn,sep="\t",quote=FALSE,row.names=FALSE)
        }
        if(nrow(z$anchors)) {
          fn <- file.path(outdir,"mechanism",paste0(target,"_trans_anchor_candidates.tsv"))
          utils::write.table(z$anchors,fn,sep="\t",quote=FALSE,row.names=FALSE)
        }
      }
      if(is.data.frame(res$results) && nrow(res$results)) {
        fc <- .cater_lean_fit_row(res$results,target,"cis")
        ft <- .cater_lean_fit_row(res$results,target,"trans")
        fa <- .cater_lean_fit_row(res$results,target,"augmented")
        if(!is.null(fc)) res$targets$cis_effective_F[i] <- fc$effective_F
        if(!is.null(ft)) res$targets$trans_effective_F[i] <- ft$effective_F
        if(!is.null(fa)) res$targets$augmented_effective_F[i] <- fa$effective_F
        res$targets$primary_effective_F[i] <- res$targets$cis_effective_F[i]
        if(!is.null(fc) && !is.null(fa) && identical(fc$status,"OK") && identical(fa$status,"OK") &&
           is.finite(fc$se) && fc$se>0 && is.finite(fa$se))
          res$targets$augmented_se_reduction[i] <- 1-fa$se/fc$se
        if(!is.null(fc) && !is.null(fa) && is.finite(fa$precision_information) && fa$precision_information>0) {
          base <- if(is.finite(fc$precision_information)) fc$precision_information else 0
          d <- fa$precision_information-base
          if(is.finite(d) && d>=-sqrt(.Machine$double.eps)*max(1,abs(fa$precision_information),abs(base)))
            res$targets$augmented_precision_fraction[i] <- max(0,min(1,d/fa$precision_information))
        }
      }
    }
    # Backward-compatible generic field now reflects the actual primary model.
    res$targets$effective_F <- res$targets$primary_effective_F
    res$targets$primary_model <- ifelse(is.finite(res$targets$primary_p),"cis",NA_character_)
  }
  if(is.data.frame(res$results) && nrow(res$results))
    utils::write.table(res$results,file.path(outdir,"cater_mr_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  if(is.data.frame(res$targets) && nrow(res$targets))
    utils::write.table(res$targets,file.path(outdir,"cater_mr_target_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  res
}

# Public Lean v0.9 wrapper. Historical positional slots are preserved exactly;
# new Lean controls are appended at the end. Advanced network analyses are off by
# default and cannot replace the cis primary estimate.
cater_mr <- function(grn,eqtl_dir,outcome,gene_annotation=NULL,ld_bfile,
                     manc_cojo_bin="manc_cojo",targets=NULL,cell_type=NA_character_,trait=NA_character_,
                     cis_window=1e6,tf_window=cis_window,cojo_p=5e-8,cojo_wind_kb=10000L,
                     cojo_collinear=0.9,cojo_threads=1L,qtl_n=NULL,
                     sibling_fdr=0.05,enable_sibling_screen=FALSE,enable_mvmr=FALSE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2,
                     exclude_mhc=TRUE,trans_eqtl=NULL,trans_gene_col=NULL,instrument_p=5e-8,
                     trans_reporting_p=5e-8,ld_clump_r2=.01,ld_clump_kb=10000L,ld_threads=1L,
                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,
                     allow_network_primary=FALSE,sibling_screen_independent=FALSE,
                     input_manifest=NULL,strict_input_contract=FALSE,
                     tf_anchor_p=instrument_p,max_reported_trans_targets=1L,
                     trans_set=c("core","extended")) {
  trans_set <- match.arg(trans_set)
  if(length(tf_anchor_p)!=1L||!is.finite(tf_anchor_p)||tf_anchor_p<=0||tf_anchor_p>=1)
    .cater_stop("tf_anchor_p must be in (0,1)")
  if(length(max_reported_trans_targets)!=1L||!is.finite(max_reported_trans_targets)||
     max_reported_trans_targets<1||max_reported_trans_targets!=as.integer(max_reported_trans_targets))
    .cater_stop("max_reported_trans_targets must be a positive integer")
  if(!missing(primary_policy) && !identical(match.arg(primary_policy),"cis_anchor"))
    warning("primary_policy is compatibility-only in Lean v0.9; cis remains the primary model",call.=FALSE)
  if(isTRUE(allow_network_primary))
    warning("allow_network_primary is ignored in Lean v0.9; network fits are sensitivity-only",call.=FALSE)

  old_active <- .CATER_LEAN_STATE$active
  old_config <- .CATER_LEAN_STATE$config
  old_diag <- .CATER_LEAN_STATE$diagnostics
  on.exit({
    .CATER_LEAN_STATE$active <- old_active
    .CATER_LEAN_STATE$config <- old_config
    .CATER_LEAN_STATE$diagnostics <- old_diag
  },add=TRUE)
  .CATER_LEAN_STATE$active <- TRUE
  .CATER_LEAN_STATE$config <- list(tf_anchor_p=tf_anchor_p,
                                   max_reported_trans_targets=as.integer(max_reported_trans_targets),
                                   trans_set=trans_set)
  .CATER_LEAN_STATE$diagnostics <- list()

  mc <- as.list(match.call(expand.dots=FALSE))[-1L]
  mc[c("tf_anchor_p","max_reported_trans_targets","trans_set")] <- NULL
  if(!"enable_sibling_screen" %in% names(mc)) mc$enable_sibling_screen <- FALSE
  if(!"enable_mvmr" %in% names(mc)) mc$enable_mvmr <- FALSE
  if(!"primary_policy" %in% names(mc)) mc$primary_policy <- "cis_anchor"
  res <- do.call(.cater_mr_v08,mc,envir=parent.frame())
  .cater_lean_postprocess(res,outdir,tf_anchor_p,max_reported_trans_targets,trans_set)
}
