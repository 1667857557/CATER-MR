source("CATER_MR.R")

required <- c("README.md","EQTL_INPUT.md","SCMORE_GRN.md","MATHEMATICS.md",
              "R/CATER_MR_core.R","R/LD_helpers.R")
stopifnot(all(file.exists(required)))

obsolete_files <- c("METHODOLOGY_V08.md","LD_DIAGNOSIS.md",
                    "R/CATER_MR_core_v05.R","R/LD_diagnosis_precojo.R")
stopifnot(!any(file.exists(obsolete_files)))

for (f in c("README.md","EQTL_INPUT.md","SCMORE_GRN.md","CATER_MR.R","R/CATER_MR_core.R","R/LD_helpers.R")) {
  x <- paste(readLines(f,warn=FALSE),collapse="\n")
  stopifnot(!grepl("pre-COJO|Manc-COJO|enable_ld_diagnosis|CATER_MR_core_v05|LD_diagnosis_precojo",x,ignore.case=TRUE))
}

fml <- names(formals(cater_mr))
stopifnot("ld_threads" %in% fml)
stopifnot(!any(c("manc_cojo_bin","cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads",
                 "enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z") %in% fml))
stopifnot(exists(".cater_plink_ld",mode="function"))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))

cat("CATER-MR documentation/API cleanup smoke tests passed\n")
