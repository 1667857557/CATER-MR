# CATER-MR stable entrypoint
#
# v0.6 keeps the previously tested v0.5 estimator in a frozen module and applies
# a focused hardening layer. This makes the evidence-safety changes reviewable
# without duplicating or silently rewriting the complete estimator.

source(file.path("R", "CATER_MR_core_v05.R"), local = FALSE)
source(file.path("R", "CATER_MR_v06_overrides.R"), local = FALSE)
