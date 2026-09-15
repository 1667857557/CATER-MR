from pathlib import Path

core = Path("R/CATER_MR_core.R")
s = core.read_text()

s = s.replace('.CATER_VERSION <- "0.8.1"', '.CATER_VERSION <- "0.8.2"', 1)

old_sig = '''cater_mr <- function(grn,eqtl_dir,outcome,gene_annotation=NULL,ld_bfile,
                     targets=NULL,cell_type=NA_character_,trait=NA_character_,
                     cis_window=1e6,tf_window=cis_window,qtl_n=NULL,
                     sibling_fdr=0.05,enable_sibling_screen=TRUE,enable_mvmr=TRUE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     exclude_mhc=TRUE,trans_eqtl=NULL,trans_gene_col=NULL,instrument_p=5e-8,
                     trans_reporting_p=5e-8,ld_clump_r2=.01,ld_clump_kb=10000L,ld_threads=1L,
                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,
                     allow_network_primary=FALSE,sibling_screen_independent=FALSE,
                     input_manifest=NULL,strict_input_contract=FALSE) {'''

new_sig = '''cater_mr <- function(grn,eqtl_dir,outcome,gene_annotation=NULL,ld_bfile,
                     manc_cojo_bin="manc_cojo",targets=NULL,cell_type=NA_character_,trait=NA_character_,
                     cis_window=1e6,tf_window=cis_window,cojo_p=5e-8,cojo_wind_kb=10000L,
                     cojo_collinear=0.9,cojo_threads=1L,qtl_n=NULL,
                     sibling_fdr=0.05,enable_sibling_screen=TRUE,enable_mvmr=TRUE,
                     exposure_corr=NULL,min_cond_F=10,mvmr_max_r2=NULL,mvmr_max_condition=1e4,
                     outdir="CATER_MR_results",drop_palindromic=TRUE,verbose=TRUE,
                     primary_policy=c("cis_anchor","screened_cater"),plink_bin="plink",
                     enable_ld_diagnosis=TRUE,ld_diag_loglr=2,ld_diag_abs_z=2,
                     exclude_mhc=TRUE,trans_eqtl=NULL,trans_gene_col=NULL,instrument_p=5e-8,
                     trans_reporting_p=5e-8,ld_clump_r2=.01,ld_clump_kb=10000L,ld_threads=1L,
                     cross_effect_lookup=NULL,accept_experimental_conditional_f=FALSE,
                     allow_network_primary=FALSE,sibling_screen_independent=FALSE,
                     input_manifest=NULL,strict_input_contract=FALSE) {'''

if old_sig not in s:
    raise SystemExit("current cater_mr signature not found")
s = s.replace(old_sig, new_sig, 1)

marker = '''  primary_policy<-match.arg(primary_policy)\n'''
compat = '''  legacy_args<-c(
    if(!missing(manc_cojo_bin)) "manc_cojo_bin",
    if(!missing(cojo_p)) "cojo_p",
    if(!missing(cojo_wind_kb)) "cojo_wind_kb",
    if(!missing(cojo_collinear)) "cojo_collinear",
    if(!missing(cojo_threads)) "cojo_threads",
    if(!missing(enable_ld_diagnosis)) "enable_ld_diagnosis",
    if(!missing(ld_diag_loglr)) "ld_diag_loglr",
    if(!missing(ld_diag_abs_z)) "ld_diag_abs_z")
  if(length(legacy_args)) warning(sprintf("Deprecated compatibility argument(s) ignored by CATER-MR v0.8: %s",paste(legacy_args,collapse=", ")),call.=FALSE)
  primary_policy<-match.arg(primary_policy)
'''
if marker not in s:
    raise SystemExit("cater_mr body marker not found")
s = s.replace(marker, compat, 1)
core.write_text(s)

# Core hardening: preserve old positional slots while keeping old implementations absent.
p = Path("tests/core_hardening_smoke.R")
t = p.read_text()
old = '''# Current LD helper is loaded; obsolete COJO/SuSiE-RSS API is absent.
fml <- names(formals(cater_mr))
stopifnot("ld_threads" %in% fml)
stopifnot(exists(".cater_plink_ld",mode="function"),exists(".cater_read_plink_bim",mode="function"))
stopifnot(!any(c("manc_cojo_bin","cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads",
                 "enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z") %in% fml))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))
'''
new = '''# Current LD helper is loaded; obsolete implementations are absent, while
# deprecated formal slots preserve the pre-cleanup positional ABI.
fml <- names(formals(cater_mr))
legacy_prefix <- c("grn","eqtl_dir","outcome","gene_annotation","ld_bfile",
                   "manc_cojo_bin","targets","cell_type","trait","cis_window","tf_window",
                   "cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads","qtl_n")
stopifnot(identical(fml[seq_along(legacy_prefix)],legacy_prefix))
plink_i <- match("plink_bin",fml)
stopifnot(identical(fml[(plink_i+1L):(plink_i+3L)],c("enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z")))
stopifnot("ld_threads" %in% fml)
stopifnot(exists(".cater_plink_ld",mode="function"),exists(".cater_read_plink_bim",mode="function"))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))
'''
if old not in t:
    raise SystemExit("core hardening legacy block not found")
p.write_text(t.replace(old,new,1))

# Documentation/API cleanup: compatibility formals are allowed; obsolete behavior remains absent.
p = Path("tests/documentation_cleanup_smoke.R")
t = p.read_text()
old_scan = '''for (f in c("README.md","EQTL_INPUT.md","SCMORE_GRN.md","CATER_MR.R","R/CATER_MR_core.R","R/LD_helpers.R")) {
  x <- paste(readLines(f,warn=FALSE),collapse="\\n")
  stopifnot(!grepl("pre-COJO|Manc-COJO|enable_ld_diagnosis|CATER_MR_core_v05|LD_diagnosis_precojo",x,ignore.case=TRUE))
}
'''
new_scan = '''for (f in c("README.md","EQTL_INPUT.md","SCMORE_GRN.md","CATER_MR.R","R/LD_helpers.R")) {
  x <- paste(readLines(f,warn=FALSE),collapse="\\n")
  stopifnot(!grepl("pre-COJO|Manc-COJO|enable_ld_diagnosis|CATER_MR_core_v05|LD_diagnosis_precojo",x,ignore.case=TRUE))
}
core_txt <- paste(readLines("R/CATER_MR_core.R",warn=FALSE),collapse="\\n")
stopifnot(!grepl("pre-COJO|Manc-COJO|CATER_MR_core_v05|LD_diagnosis_precojo|\\.cater_susie_ld_diagnosis|\\.cater_read_manc_ldr",
                 core_txt,ignore.case=TRUE))
'''
if old_scan not in t:
    raise SystemExit("documentation cleanup scan block not found")
t = t.replace(old_scan,new_scan,1)
old = '''fml <- names(formals(cater_mr))
stopifnot("ld_threads" %in% fml)
stopifnot(!any(c("manc_cojo_bin","cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads",
                 "enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z") %in% fml))
stopifnot(exists(".cater_plink_ld",mode="function"))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))
'''
new = '''fml <- names(formals(cater_mr))
legacy_prefix <- c("grn","eqtl_dir","outcome","gene_annotation","ld_bfile",
                   "manc_cojo_bin","targets","cell_type","trait","cis_window","tf_window",
                   "cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads","qtl_n")
stopifnot(identical(fml[seq_along(legacy_prefix)],legacy_prefix))
plink_i <- match("plink_bin",fml)
stopifnot(identical(fml[(plink_i+1L):(plink_i+3L)],c("enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z")))
stopifnot("ld_threads" %in% fml)
stopifnot(exists(".cater_plink_ld",mode="function"))
stopifnot(!exists(".cater_read_manc_ldr",mode="function"))
stopifnot(!exists(".cater_susie_ld_diagnosis",mode="function"))
'''
if old not in t:
    raise SystemExit("documentation cleanup legacy block not found")
p.write_text(t.replace(old,new,1))

# Add an explicit regression that positional slots bind as before without executing the analysis.
p = Path("tests/positional_compat_smoke.R")
p.write_text('''source("CATER_MR.R")\n\nfml <- names(formals(cater_mr))\nexpected <- c("grn","eqtl_dir","outcome","gene_annotation","ld_bfile",\n              "manc_cojo_bin","targets","cell_type","trait","cis_window","tf_window",\n              "cojo_p","cojo_wind_kb","cojo_collinear","cojo_threads","qtl_n")\nstopifnot(identical(fml[seq_along(expected)],expected))\npi <- match("plink_bin",fml)\nstopifnot(identical(fml[(pi+1L):(pi+3L)],c("enable_ld_diagnosis","ld_diag_loglr","ld_diag_abs_z")))\n\n# Positional matching must assign the historical 6th/7th and COJO slots correctly.\nmc <- match.call(cater_mr, quote(cater_mr(G,E,O,A,L,"manc_cojo",c("X","Y"),"ct","trait",\n                                         1e6,1e6,5e-8,10000L,0.9,4L,500)), expand.dots=FALSE)\nstopifnot(identical(as.character(mc$manc_cojo_bin),"manc_cojo"))\nstopifnot(identical(mc$targets,quote(c("X","Y"))))\nstopifnot(identical(as.character(mc$cell_type),"ct"))\nstopifnot(identical(as.character(mc$trait),"trait"))\nstopifnot(identical(as.numeric(mc$cojo_p),5e-8))\nstopifnot(identical(as.integer(mc$cojo_wind_kb),10000L))\nstopifnot(identical(as.numeric(mc$cojo_collinear),0.9))\nstopifnot(identical(as.integer(mc$cojo_threads),4L))\nstopifnot(identical(as.numeric(mc$qtl_n),500))\n\ncat("CATER-MR positional compatibility tests passed\\n")\n''')

# Run the new regression in the main workflow.
p = Path(".github/workflows/r-smoke.yml")
t = p.read_text()
needle = '''      - name: Run documentation/API cleanup tests\n        run: Rscript tests/documentation_cleanup_smoke.R\n'''
repl = needle + '''      - name: Run positional compatibility tests\n        run: Rscript tests/positional_compat_smoke.R\n'''
if needle not in t:
    raise SystemExit("r-smoke insertion marker not found")
p.write_text(t.replace(needle,repl,1))
