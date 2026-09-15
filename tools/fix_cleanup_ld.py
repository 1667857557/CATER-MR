from pathlib import Path

core = Path("R/CATER_MR_core.R")
ldfile = Path("R/LD_diagnosis_precojo.R")
s = core.read_text()
ld = ldfile.read_text()

if '.cater_plink_ld <- function' not in s:
    a = ld.index('.cater_read_plink_bim <- function')
    b = ld.index('.cater_susie_ld_diagnosis <- function', a)
    block = ld[a:b]
    block = block.replace('.cater_validate_diag_ld(R,"pre-COJO PLINK signed LD")',
                          '.cater_validate_ld(R,"PLINK signed LD")')
    block = block.replace("required for pre-COJO LD diagnosis", "required for LD selection")
    block = block.replace("failed for LD diagnosis", "failed for LD selection")
    support = '''.cater_remove_prefix_outputs <- function(prefix){z<-Sys.glob(paste0(prefix,"*"));if(length(z))unlink(z,recursive=TRUE,force=TRUE);invisible(NULL)}\n.cater_quote_args <- function(args) vapply(as.character(args),function(z)if(grepl("[[:space:]]",z))shQuote(z)else z,character(1),USE.NAMES=FALSE)\n\n'''
    marker = '.cater_reorient_ld <- function'
    if marker not in s:
        raise SystemExit("LD helper insertion marker missing")
    s = s.replace(marker, support + block + marker, 1)

core.write_text(s)
