from pathlib import Path
import re

core = Path("R/CATER_MR_core_v05.R")
text = core.read_text()
pattern = re.compile(r'\.cater_make_regions <- function\(target, parents, annotation, cis_window, tf_window\) \{.*?\n\}\n\n\.cater_candidate_map <- function\(qtl, regions\) \{.*?\n\}\n\n(?=\.cater_write_cojo_ma)', re.S)
replacement = r'''.cater_make_regions <- function(target, parents, annotation, cis_window, tf_window) {
  ta <- annotation[annotation$symbol==target,,drop=FALSE]
  if (!nrow(ta)) return(NULL)
  reg <- data.frame(type="cis",gene=target,chr=ta$chr[1],
                    start=max(1,ta$tss[1]-cis_window),end=ta$tss[1]+cis_window,
                    locus_id="cis",stringsAsFactors=FALSE)
  if (length(parents)) {
    pa <- annotation[match(parents,annotation$symbol),,drop=FALSE]
    pa$gene <- parents
    pa <- pa[!is.na(pa$chr)&is.finite(pa$tss),,drop=FALSE]
    if (nrow(pa)) {
      tr <- data.frame(type="trans",gene=pa$gene,chr=pa$chr,
        start=pmax(1,pa$tss-tf_window),end=pa$tss+tf_window,
        locus_id=NA_character_,stringsAsFactors=FALSE)
      # Merge overlapping TF windows into physical locus components before SNP
      # assignment. Original TF windows remain separate rows so parent_tf still
      # records which TF window(s) each SNP actually occupies.
      component <- integer(nrow(tr)); next_component <- 0L
      for (ch in unique(tr$chr)) {
        oi <- which(tr$chr==ch)
        oi <- oi[order(tr$start[oi],tr$end[oi],tr$gene[oi])]
        current <- 0L; current_end <- -Inf
        for (idx in oi) {
          if (current==0L || tr$start[idx] > current_end) {
            next_component <- next_component + 1L
            current <- next_component
            current_end <- tr$end[idx]
          } else {
            current_end <- max(current_end,tr$end[idx])
          }
          component[idx] <- current
        }
      }
      for (cc in unique(component)) {
        jj <- which(component==cc)
        tr$locus_id[jj] <- paste0("TF:",paste(sort(unique(tr$gene[jj])),collapse="+"))
      }
      reg <- rbind(reg,tr)
    }
  }
  rownames(reg) <- NULL
  reg
}

.cater_candidate_map <- function(qtl, regions) {
  if (is.null(regions)||!nrow(qtl)) return(data.frame())
  cis <- regions[regions$type=="cis",,drop=FALSE]
  is_cis <- qtl$chr==cis$chr[1] & qtl$pos>=cis$start[1] & qtl$pos<=cis$end[1]
  tr <- regions[regions$type=="trans",,drop=FALSE]
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

helper = Path("R/LD_diagnosis_precojo.R")
h = helper.read_text()
gpat = re.compile(r'\.cater_ld_diagnosis_groups <- function\(q, candidate_map\) \{.*?\n\}\n\n(?=\.cater_validate_diag_ld)', re.S)
grep = r'''.cater_ld_diagnosis_groups <- function(q, candidate_map) {
  if (nrow(q)!=nrow(candidate_map)) .cater_stop("LD diagnosis grouping requires aligned qtl and candidate_map rows")
  if (!nrow(candidate_map)) return(character())
  if (any(is.na(candidate_map$locus_id) | !nzchar(candidate_map$locus_id)))
    .cater_stop("LD diagnosis requires one physical locus_id per candidate")
  paste(q$chr,candidate_map$locus_id,sep="|")
}

'''
h2, n = gpat.subn(grep, h)
assert n == 1, f"diagnosis grouping replacement count={n}"
h2 = h2.replace('rows <- data.frame(snp=ids,chr=qq$chr,locus_id=cm$locus_id,ld_locus_id=rep(ug[gidx],length(ids)),source=cm$source,',
                'rows <- data.frame(snp=ids,chr=qq$chr,locus_id=cm$locus_id,source=cm$source,')
helper.write_text(h2)

test = Path("tests/v06_hardening_smoke.R")
t = test.read_text()
tpat = re.compile(r'# 18\..*?\n(?=cat\("CATER-MR direct-core evidence-safety tests passed)', re.S)
trep = r'''# 18. Overlapping TF windows are merged by genomic interval before candidate assignment.
ann_locus <- data.frame(symbol=c("X","A","B","C"),chr=rep("1",4),
                        tss=c(1000,10000,11500,20000),stringsAsFactors=FALSE)
reg_locus <- .cater_make_regions("X",c("A","B","C"),ann_locus,cis_window=50,tf_window=1000)
tr_locus <- reg_locus[reg_locus$type=="trans",,drop=FALSE]
stopifnot(tr_locus$locus_id[tr_locus$gene=="A"]=="TF:A+B",
          tr_locus$locus_id[tr_locus$gene=="B"]=="TF:A+B",
          tr_locus$locus_id[tr_locus$gene=="C"]=="TF:C")
q_locus <- data.frame(snp=paste0("g",1:4),chr=rep("1",4),
                      pos=c(9500,10750,12000,20000),stringsAsFactors=FALSE)
cm_locus <- .cater_candidate_map(q_locus,reg_locus)
stopifnot(identical(cm_locus$parent_tf,c("A","A;B","B","C")),
          identical(cm_locus$locus_id,c("TF:A+B","TF:A+B","TF:A+B","TF:C")))
gg <- .cater_ld_diagnosis_groups(q_locus,cm_locus)
stopifnot(gg[1]==gg[2],gg[2]==gg[3],gg[4]!=gg[1])

'''
t2, n = tpat.subn(trep, t)
assert n == 1, f"test replacement count={n}"
test.write_text(t2)

doc = Path("LD_DIAGNOSIS.md")
d = doc.read_text()
needle = "CATER-MR performs an LD/summary-statistic consistency gate before every Manc-COJO selection step by default.\n"
insert = needle + "\nTrans candidates are defined by local parent-TF windows. TF windows on the same chromosome that overlap are merged into one physical locus component before candidate SNP assignment (for example, overlapping A and B windows share `locus_id=TF:A+B`). The original windows are retained for `parent_tf`, so a SNP can still be attributed to A, B, or A;B according to its actual position while LD/SuSiE treats the full overlapping component as one locus.\n"
assert needle in d
d = d.replace(needle, insert, 1)
d = d.replace("For each candidate cis/trans locus, CATER-MR:", "For each physical candidate cis/TF locus, CATER-MR:")
doc.write_text(d)
