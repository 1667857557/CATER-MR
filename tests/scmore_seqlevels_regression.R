source("SCMORE_GRN.R")
needed <- c("Seurat", "Signac", "GenomicRanges", "GenomeInfoDb")
stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly=TRUE)))

# Sequence-name aliases must fail instead of being converted into empty levels.
alias_error <- tryCatch(
  .cater_scmore_reject_seqname_aliases(c("chr1", "GL000195.1"), c("1", "chrM")),
  error=function(e)e
)
stopifnot(inherits(alias_error,"error"))
stopifnot(grepl("conflicting aliases",conditionMessage(alias_error),fixed=TRUE))
stopifnot(grepl("chr1/1",conditionMessage(alias_error),fixed=TRUE))
stopifnot(isTRUE(.cater_scmore_reject_seqname_aliases(c("chr1", "GL000195.1"), c("chr1", "chrM"))))

# Both operands have an exclusive sequence level. No private data needed.
peaks <- GenomicRanges::GRanges(c("chr1", "GL000195.1"), IRanges::IRanges(c(10,30),c(20,40)))
anno <- GenomicRanges::GRanges(c("chr1", "chrM"), IRanges::IRanges(c(12,50),c(15,60)), strand="+")
anno$type <- "exon"; anno$gene_id <- c("g1","g2")
anno$tx_id <- c("tx1","tx2")
anno$gene_name <- c("G1","G2"); anno$gene_biotype <- "protein_coding"
GenomeInfoDb::genome(anno) <- "hg38"
GenomeInfoDb::seqlengths(anno) <- c(1000,100)
GenomeInfoDb::isCircular(anno) <- c(FALSE,TRUE)
counts <- Matrix::Matrix(matrix(c(1,2,3,4),2),sparse=TRUE)
dimnames(counts) <- list(c("G1","G2"),c("c1","c2"))
obj <- Seurat::CreateSeuratObject(counts=counts)
rownames(counts) <- c("chr1-10-20","GL000195.1-30-40")
obj[["peaks"]] <- Signac::CreateChromatinAssay(counts=counts,annotation=anno,min.cells=0,min.features=0)
original <- obj
fixed <- .cater_scmore_prepare_seqlevels(obj)
fa <- Signac::Annotation(fixed[["peaks"]])
stopifnot(identical(obj, original))
annotation_rows <- function(a) {
  d <- as.data.frame(a); d$seqnames <- as.character(d$seqnames); d
}
stopifnot(identical(annotation_rows(fa),annotation_rows(anno)))
stopifnot(identical(GenomeInfoDb::seqinfo(fa)[GenomeInfoDb::seqlevels(anno)],GenomeInfoDb::seqinfo(anno)))
stopifnot(identical(.cater_scmore_prepare_seqlevels(fixed),fixed))
stopifnot(identical(rownames(fixed[["peaks"]]),rownames(obj[["peaks"]])))
delivered <- .cater_scmore_call_create_regulon(obj,list(),create_fun=function(single_cell)single_cell)
stopifnot(identical(delivered,fixed))

old <- tryCatch(suppressWarnings(GenomicRanges::subtract(peaks,anno,ignore.strand=TRUE)),error=function(e)e)
if(as.character(utils::packageVersion("GenomicRanges"))=="1.64.0")stopifnot(inherits(old,"error"))
if(inherits(old,"error"))stopifnot(grepl("Level set",conditionMessage(old),fixed=TRUE))
out <- unlist(GenomicRanges::subtract(peaks,fa,ignore.strand=TRUE))
expected <- GenomicRanges::GRanges(c("chr1","chr1","GL000195.1"),IRanges::IRanges(c(10,16,30),c(11,20,40)))
coords <- function(x)paste(as.character(GenomicRanges::seqnames(x)),IRanges::start(x),IRanges::end(x),sep=":")
stopifnot(identical(sort(coords(out)),sort(coords(expected))))
stopifnot(length(GenomicRanges::findOverlaps(out,fa,ignore.strand=TRUE))==0L)

if(requireNamespace("Pando",quietly=TRUE)) {
  before <- tryCatch(suppressWarnings(Pando::initiate_grn(obj,regions=peaks)),error=function(e)e)
  if(as.character(utils::packageVersion("GenomicRanges"))=="1.64.0" && as.character(utils::packageVersion("Pando"))=="1.1.1")stopifnot(inherits(before,"error"))
  if(inherits(before,"error"))stopifnot(grepl("Level set",conditionMessage(before),fixed=TRUE))
  after <- Pando::initiate_grn(fixed,regions=peaks)
  stopifnot(identical(sort(coords(after@grn@regions@ranges)),sort(coords(expected))))
  cat("Real Pando::initiate_grn integration passed\n")
}
cat("scMORE seqlevels regression passed\n")
