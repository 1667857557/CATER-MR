source("CATER_MR.R")

required <- c("README.md","EQTL_INPUT.md","SCMORE_GRN.md","MATHEMATICS.md",
              "R/CATER_MR_core.R","R/CATER_MR_lean.R","R/LD_helpers.R")
stopifnot(all(file.exists(required)))

obsolete_files <- c("METHODOLOGY_V08.md","LD_DIAGNOSIS.md",
                    "R/CATER_MR_core_v05.R","R/LD_diagnosis_precojo.R")
stopifnot(!any(file.exists(obsolete_files)))

for (f in c("README.md","EQTL_INPUT.md","SCMORE_GRN.md","CATER_MR.R","R/LD_helpers.R")) {
  x <- paste(readLines(f,warn=FALSE),collapse="\n")
  stopifnot(!grepl("pre-COJO|Manc-COJO|enable_ld_diagnosis|CATER_MR_core_v05|LD_diagnosis_precojo",x,ignore.case=TRUE))
}
core_txt <- paste(readLines("R/CATER_MR_core.R",warn=FALSE),collapse="\n")
stopifnot(!grepl("pre-COJO|Manc-COJO|CATER_MR_core_v05|LD_diagnosis_precojo",core_txt,ignore.case=TRUE))
stopifnot(!grepl(".cater_susie_ld_diagnosis",core_txt,fixed=TRUE))
stopifnot(!grepl(".cater_read_manc_ldr",core_txt,fixed=TRUE))

fml <- names(formals(cater_mr))
legacy_prefix <- c("grn","eqtl_dir","outcome","gene_annotation","ld_bfile",
                   "manc_cojo_bin","targets","cell_type","trait","cis_window","tf_window",
                   "cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads","qtl_n")
stopifnot(identical(fml[seq_along(legacy_prefix)],legacy_prefix))
plink_i <- match("plink_bin",fml)
stopifnot(identical(fml[(plink_i+1L):(plink_i+3L)],c("enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z")))
stopifnot("ld_threads" %in% fml)
stopifnot(all(c("tf_anchor_p","max_reported_trans_targets","trans_set") %in% fml))
stopifnot(identical(formals(cater_mr)$enable_sibling_screen,FALSE))
stopifnot(identical(formals(cater_mr)$enable_mvmr,FALSE))
stopifnot(exists(".cater_qualify_trans_candidates",mode="function"))
stopifnot(exists(".cater_plink_ld",mode="function"))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))

cat("CATER-MR documentation/API cleanup smoke tests passed\n")
