from pathlib import Path
import re

core = Path("R/CATER_MR_core_v05.R")
text = core.read_text()
text = text.replace('.CATER_VERSION <- "0.7.1"', '.CATER_VERSION <- "0.7.2"', 1)
pattern = re.compile(r'\.cater_make_regions <- function\(target, parents, annotation, cis_window, tf_window\) \{.*?\n\}\n\n\.cater_candidate_map <- function\(qtl, regions\) \{.*?\n\}\n\n(?=\.cater_write_cojo_ma)', re.S)
replacement = r'''.cater_make_regions <- function(target, parents, annotation, cis_window, tf_window) {
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

'''
text2, n = pattern.subn(replacement, text)
assert n == 1, f"core region/candidate replacement count={n}"
core.write_text(text2)

# Extend the existing hardening test with target-cis/TF overlap semantics.
test = Path("tests/v06_hardening_smoke.R")
t = test.read_text()
marker = 'cat("CATER-MR direct-core evidence-safety tests passed\\n")'
assert marker in t
block = r'''# 19. A TF window connected to the target cis window is absorbed into the full cis locus.
# A directly overlaps target cis; B overlaps A, so the whole connected component is cis.
ann_cis_tf <- data.frame(symbol=c("X","A","B","C"),chr=rep("1",4),
                         tss=c(10000,11500,13000,20000),stringsAsFactors=FALSE)
reg_cis_tf <- .cater_make_regions("X",c("A","B","C"),ann_cis_tf,cis_window=1000,tf_window=1000)
stopifnot(reg_cis_tf$locus_id[reg_cis_tf$gene=="X"]=="cis",
          reg_cis_tf$locus_id[reg_cis_tf$gene=="A"]=="cis",
          reg_cis_tf$locus_id[reg_cis_tf$gene=="B"]=="cis",
          reg_cis_tf$locus_id[reg_cis_tf$gene=="C"]=="TF:C")
q_cis_tf <- data.frame(snp=paste0("ct",1:4),chr=rep("1",4),
                       pos=c(9500,11500,13000,20000),stringsAsFactors=FALSE)
cm_cis_tf <- .cater_candidate_map(q_cis_tf,reg_cis_tf)
stopifnot(identical(cm_cis_tf$source,c("cis","cis","cis","trans")),
          identical(cm_cis_tf$locus_id,c("cis","cis","cis","TF:C")),
          identical(cm_cis_tf$parent_tf,c("","","","C")))
gg_cis_tf <- .cater_ld_diagnosis_groups(q_cis_tf,cm_cis_tf)
stopifnot(length(unique(gg_cis_tf[1:3]))==1L,gg_cis_tf[4]!=gg_cis_tf[1])

'''
t = t.replace(marker, block + marker, 1)
test.write_text(t)

# Update method documentation to make cis-component precedence explicit.
doc = Path("LD_DIAGNOSIS.md")
d = doc.read_text()
old = "Trans candidates are defined by local parent-TF windows. TF windows on the same chromosome that overlap are merged into one physical locus component **before candidate SNP assignment**. For example, overlapping A and B windows share `locus_id=TF:A+B`. The original TF windows remain available for biological attribution, so a SNP can still have `parent_tf=A`, `B`, or `A;B` according to its actual genomic position."
new = "Candidate loci are defined from the target cis window and all parent-TF windows together. Overlapping intervals on the same chromosome are merged into connected physical components before candidate SNP assignment. Any component containing the target cis window is treated in full as the cis locus (`locus_id=cis`), including TF-window extensions connected to it; those SNPs are not reused as trans instruments. Components containing only TF windows remain local trans loci such as `TF:A+B`."
assert old in d
d = d.replace(old,new,1)
d = d.replace("1. merges overlapping parent-TF windows by genomic interval into connected physical loci; the target cis window remains its own locus;",
              "1. merges the target cis window and parent-TF windows by genomic interval into connected physical loci; any component containing the target cis window is the full cis locus, while TF-only components remain trans loci;",1)
d = d.replace("`locus_id` is the physical cis/TF locus used for LD diagnosis and is also propagated with the retained instruments. `parent_tf` remains the biological attribution field. Thus overlapping TF windows are one locus for LD/COJO bookkeeping without erasing whether an individual SNP lies in A, B, or both windows.",
              "`locus_id` is the physical cis/TF locus used for LD diagnosis and is also propagated with the retained instruments. TF-only components retain `parent_tf` for biological attribution. If a TF component overlaps the target cis component, the entire connected component is classified as cis and its SNPs are deliberately not assigned a trans parent, avoiding cis/trans double interpretation within one local LD neighbourhood.",1)
doc.write_text(d)
