source("CATER_MR.R")

fml <- names(formals(cater_mr))
expected <- c("grn","eqtl_dir","outcome","gene_annotation","ld_bfile",
              "manc_cojo_bin","targets","cell_type","trait","cis_window","tf_window",
              "cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads","qtl_n")
stopifnot(identical(fml[seq_along(expected)],expected))
pi <- match("plink_bin",fml)
stopifnot(identical(fml[(pi+1L):(pi+3L)],c("enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z")))

# Positional matching must assign the historical 6th/7th and COJO slots correctly.
mc <- match.call(cater_mr, quote(cater_mr(G,E,O,A,L,"manc_cojo",c("X","Y"),"ct","trait",
                                         1e6,1e6,5e-8,10000L,0.9,4L,500)), expand.dots=FALSE)
stopifnot(identical(as.character(mc$manc_cojo_bin),"manc_cojo"))
stopifnot(identical(mc$targets,quote(c("X","Y"))))
stopifnot(identical(as.character(mc$cell_type),"ct"))
stopifnot(identical(as.character(mc$trait),"trait"))
stopifnot(identical(as.numeric(mc$cojo_p),5e-8))
stopifnot(identical(as.integer(mc$cojo_wind_kb),10000L))
stopifnot(identical(as.numeric(mc$cojo_collinear),0.9))
stopifnot(identical(as.integer(mc$cojo_threads),4L))
stopifnot(identical(as.numeric(mc$qtl_n),500))

cat("CATER-MR positional compatibility tests passed\n")
