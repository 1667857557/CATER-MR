from pathlib import Path

p = Path("R/CATER_MR_core.R")
s = p.read_text()
old = '''  legacy_args<-c(
    if(!missing(manc_cojo_bin)) "manc_cojo_bin",
    if(!missing(cojo_p)) "cojo_p",
    if(!missing(cojo_wind_kb)) "cojo_wind_kb",
    if(!missing(cojo_collinear)) "cojo_collinear",
    if(!missing(cojo_threads)) "cojo_threads",
    if(!missing(enable_ld_diagnosis)) "enable_ld_diagnosis",
    if(!missing(ld_diag_loglr)) "ld_diag_loglr",
    if(!missing(ld_diag_abs_z)) "ld_diag_abs_z")
  if(length(legacy_args)) warning(sprintf("Deprecated compatibility argument(s) ignored by CATER-MR v0.8: %s",paste(legacy_args,collapse=", ")),call.=FALSE)
'''
new = '''  mapped_cojo_threads<-!missing(cojo_threads)&&missing(ld_threads)
  if(mapped_cojo_threads){
    ld_threads<-cojo_threads
    warning("Deprecated cojo_threads is mapped to ld_threads for compatibility",call.=FALSE)
  }
  legacy_args<-c(
    if(!missing(manc_cojo_bin)) "manc_cojo_bin",
    if(!missing(cojo_p)) "cojo_p",
    if(!missing(cojo_wind_kb)) "cojo_wind_kb",
    if(!missing(cojo_collinear)) "cojo_collinear",
    if(!missing(cojo_threads)&&!mapped_cojo_threads) "cojo_threads",
    if(!missing(enable_ld_diagnosis)) "enable_ld_diagnosis",
    if(!missing(ld_diag_loglr)) "ld_diag_loglr",
    if(!missing(ld_diag_abs_z)) "ld_diag_abs_z")
  if(length(legacy_args)) warning(sprintf("Deprecated compatibility argument(s) ignored by CATER-MR v0.8: %s",paste(legacy_args,collapse=", ")),call.=FALSE)
'''
if old not in s:
    raise SystemExit("legacy compatibility block not found")
p.write_text(s.replace(old,new,1))

p = Path("tests/positional_compat_smoke.R")
t = p.read_text()
needle = 'stopifnot(identical(as.numeric(mc$qtl_n),500))\n'
add = '''stopifnot(identical(as.numeric(mc$qtl_n),500))
core_txt <- paste(readLines("R/CATER_MR_core.R",warn=FALSE),collapse="\\n")
stopifnot(grepl("ld_threads<-cojo_threads",core_txt,fixed=TRUE))
stopifnot(grepl("missing(ld_threads)",core_txt,fixed=TRUE))
'''
if needle not in t:
    raise SystemExit("positional regression insertion marker not found")
p.write_text(t.replace(needle,add,1))
