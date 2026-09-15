from pathlib import Path

p = Path("tests/documentation_cleanup_smoke.R")
t = p.read_text()
old = '''core_txt <- paste(readLines("R/CATER_MR_core.R",warn=FALSE),collapse="\\n")
stopifnot(!grepl("pre-COJO|Manc-COJO|CATER_MR_core_v05|LD_diagnosis_precojo|\\.cater_susie_ld_diagnosis|\\.cater_read_manc_ldr",
                 core_txt,ignore.case=TRUE))
'''
new = '''core_txt <- paste(readLines("R/CATER_MR_core.R",warn=FALSE),collapse="\\n")
stopifnot(!grepl("pre-COJO|Manc-COJO|CATER_MR_core_v05|LD_diagnosis_precojo",core_txt,ignore.case=TRUE))
stopifnot(!grepl(".cater_susie_ld_diagnosis",core_txt,fixed=TRUE))
stopifnot(!grepl(".cater_read_manc_ldr",core_txt,fixed=TRUE))
'''
if old not in t:
    raise SystemExit("generated cleanup scan block not found")
p.write_text(t.replace(old,new,1))
