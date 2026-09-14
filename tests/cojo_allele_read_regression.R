source('CATER_MR.R')
local({
  p <- tempfile(fileext='.jma.cojo')
  on.exit(unlink(p))
  for (alleles in list(c('T','C'), c('C','T'), c('T','G'))) {
    writeLines(c('Chr SNP bp A1 A2 freq b se p n bJ bJ_se pJ',
      paste('1 rs1 100',paste(alleles,collapse=' '),'.2 .5 .1 1e-9 342 .5 .1 1e-9')),p)
    ref <- .cater_read_jma_ref(p)
    stopifnot(identical(ref$ld_a1,alleles[1]),identical(ref$ld_a2,alleles[2]))
    d <- data.frame(snp='rs1',a1=alleles[1],a2=alleles[2],beta=.5,eaf=.2)
    aligned <- .cater_align_to_ld(d,ref)
    stopifnot(nrow(aligned)==1L,aligned$beta==.5,aligned$eaf==.2)
    d$a1 <- alleles[2];d$a2 <- alleles[1]
    aligned <- .cater_align_to_ld(d,ref)
    stopifnot(nrow(aligned)==1L,aligned$beta== -.5,aligned$eaf==.8)
  }
  writeLines(c('SNP A1 A2','rs1 T C','rs2 T G'),p)
  stopifnot(identical(.cater_read_jma_ref(p)$ld_a1,c('T','T')))
  writeLines('SNP A1 A2',p)
  stopifnot(nrow(.cater_read_jma_ref(p))==0L)
})
cat('COJO allele reader regression passed\n')
