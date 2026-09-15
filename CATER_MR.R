# CATER-MR stable entrypoint
# Core functions are defined once. LD-reference and legacy SuSiE-RSS diagnostic
# helpers are loaded separately. The v0.8 default IV-selection path is
# significance + GRN/physical-locus gating + signed-LD pruning, not COJO.

source(file.path("R", "CATER_MR_core_v05.R"), local = FALSE)
source(file.path("R", "LD_diagnosis_precojo.R"), local = FALSE)
