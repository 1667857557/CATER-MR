source("CATER_MR.R")
source("SCMORE_GRN.R")

# CATER-ready graph must collapse repeated scMORE TF-target evidence rows,
# remove self loops, and append the exact coordinate columns consumed by cater_mr().
raw <- data.frame(
  TF = c("TF1", "TF1", "X", "TF2"),
  Target = c("X", "X", "X", "Y"),
  Regions = c("chr1-10-20", "chr1-30-40", "chr1-50-60", "chr2-10-20"),
  Pval = c(.01, .02, .03, .04),
  stringsAsFactors = FALSE
)
ann <- data.frame(
  symbol = c("TF1", "X", "TF2", "Y"),
  chr = c("1", "1", "2", "2"),
  tss = c(1000, 5000, 2000, 8000),
  stringsAsFactors = FALSE
)
g <- .cater_scmore_make_cater_grn(raw, ann, "Microglia")
stopifnot(nrow(g) == 2L)
stopifnot(all(c("TF","Target","TF_chr","TF_tss","Target_chr","Target_tss","cell_type") %in% names(g)))
stopifnot(sum(g$TF == "TF1" & g$Target == "X") == 1L)
stopifnot(!any(g$TF == g$Target))
stopifnot(attr(g, "n_unique_raw_edges") == 3L)
stopifnot(attr(g, "n_self_loops_dropped") == 1L)
stopifnot(identical(unique(g$cell_type), "Microglia"))

# The generated graph must satisfy the core CATER-MR GRN + coordinate contract
# without a separate gene_annotation argument.
g_core <- .cater_standardize_grn(g)
a_core <- .cater_annotation_from_grn(g)
stopifnot(nrow(g_core) == 2L)
stopifnot(all(c("TF1", "X", "TF2", "Y") %in% a_core$symbol))
stopifnot(a_core$tss[a_core$symbol == "TF1"] == 1000)

# Empty scMORE output is a valid empty biological result, not a scalar-assignment error.
empty_raw <- data.frame(TF=character(), Target=character(), Regions=character(), Pval=numeric())
empty_g <- .cater_scmore_make_cater_grn(empty_raw, ann, "EmptyCell")
stopifnot(nrow(empty_g) == 0L)
stopifnot(identical(names(empty_g), c("TF","Target","TF_chr","TF_tss","Target_chr","Target_tss","cell_type")))

# Missing coordinates are strict by default because parent TF TSS is required
# to define the one-hop trans locus. Explicit drop mode is available but counted.
bad_ann <- ann[ann$symbol != "TF2",]
err <- try(.cater_scmore_make_cater_grn(raw, bad_ann, "Microglia"), silent=TRUE)
stopifnot(inherits(err, "try-error"))
g_drop <- .cater_scmore_make_cater_grn(raw, bad_ann, "Microglia", missing_coordinate="drop")
stopifnot(attr(g_drop, "n_missing_coordinate_edges") == 1L)

# Ambiguous transcript/TSS mappings are never silently resolved.
amb <- rbind(ann, data.frame(symbol="X", chr="1", tss=7000))
err2 <- try(.cater_scmore_standardize_gene_annotation(amb), silent=TRUE)
stopifnot(inherits(err2, "try-error"))

# Explicit incompatible genomes are rejected; hg38 and unknown metadata are accepted.
stopifnot(.cater_scmore_check_genome_values("hg38") == "hg38")
stopifnot(.cater_scmore_check_genome_values(character()) == "UNKNOWN")
err3 <- try(.cater_scmore_check_genome_values("hg19"), silent=TRUE)
stopifnot(inherits(err3, "try-error"))

# Sanitized output names must remain unique even when labels collapse to the same base name.
nm <- .cater_scmore_output_names(c("CD8+ T", "CD8 T"))
stopifnot(length(unique(nm)) == 2L)

# Persisted error messages cannot contain literal tabs/newlines that would corrupt TSV rows.
msg <- .cater_scmore_clean_message("first\tsecond\nthird")
stopifnot(!grepl("[\t\r\n]", msg))

# CATER-MR must explicitly override scMORE's phastCons-only upstream default so that
# Pando::initiate_grn() receives the genomic union of conserved and SCREEN hg38 regions.
default_expr <- paste(deparse(formals(cater_build_scmore_grn)$conserved_regions), collapse="")
stopifnot(identical(default_expr, ".cater_scmore_default_regions()"))
region_loader <- paste(deparse(body(.cater_scmore_default_regions)), collapse=" ")
stopifnot(grepl("GenomicRanges::union", region_loader, fixed=TRUE))
stopifnot(!"conserved_regions" %in% names(.cater_scmore_create_args(conserved_regions=NULL)))

if(requireNamespace("Pando",quietly=TRUE) && requireNamespace("GenomicRanges",quietly=TRUE)) {
  default_regions <- .cater_scmore_default_regions()
  e <- new.env(parent=baseenv())
  suppressWarnings(utils::data(
    list=c("phastConsElements20Mammals.UCSC.hg38","SCREEN.ccRE.UCSC.hg38"),
    package="Pando", envir=e
  ))
  expected_regions <- GenomicRanges::union(
    e$phastConsElements20Mammals.UCSC.hg38,
    e$SCREEN.ccRE.UCSC.hg38
  )
  stopifnot(inherits(default_regions, "GRanges"))
  stopifnot(identical(default_regions, expected_regions))
} else {
  default_regions_err <- try(.cater_scmore_default_regions(), silent=TRUE)
  stopifnot(inherits(default_regions_err, "try-error"))
}

cat("scMORE -> CATER-MR GRN contract smoke tests passed\n")
