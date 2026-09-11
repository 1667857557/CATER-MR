# CATER-MR stable entrypoint
# Core functions are defined once; LD-diagnosis helpers are loaded separately
# and called directly by the core COJO implementation.

source(file.path("R", "CATER_MR_core_v05.R"), local = FALSE)
source(file.path("R", "LD_diagnosis_precojo.R"), local = FALSE)
