# write_session_info.R — Capture exact package versions for reproducibility.
# Run AFTER sourcing 00_setup.R (so all analysis packages are loaded), e.g.:
#   source("00_setup.R"); source("write_session_info.R")
# Commit the resulting reproducibility/sessionInfo.txt to the repository.

dir.create(here::here("reproducibility"), showWarnings = FALSE)
writeLines(
  capture.output(sessionInfo()),
  here::here("reproducibility", "sessionInfo.txt")
)
message("Wrote reproducibility/sessionInfo.txt")
