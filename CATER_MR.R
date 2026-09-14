# CATER-MR stable entrypoint
# Core functions are defined once. LD-reference helpers are loaded separately and
# used by the v0.8 significance + GRN + LD-selection path; legacy SuSiE-RSS/COJO
# helpers remain available for compatibility but are not part of default IV selection.

source(file.path("R", "CATER_MR_core_v05.R"), local = FALSE)
source(file.path("R", "LD_diagnosis_precojo.R"), local = FALSE)
