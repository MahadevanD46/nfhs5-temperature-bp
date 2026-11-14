# =============================================================================
# run_all.R — run revision scripts 01-08 in order, each in a fresh R session
# (fresh sessions keep memory low on a laptop). From your project folder:
#     source("code/B_revision/run_all.R")
# Set CDS_KEY in your environment before 02a / 02c if those files still need
# downloading:  Sys.setenv(CDS_KEY = "your-token")   (never save it in a script)
# =============================================================================
steps <- c(
  "01_flow_validation.R",        # flow, exclusions, validation, date handling
  "02a_download_humidity_rain.R",# ERA5-Land dewpoint + precipitation (skip if done)
  "02b_extra_exposures.R",       # new exposures -> analysis_dataset_v2.rds
  "02c_timezone_check.R",        # hourly/daily/time-zone check
  "02d_imd_validation.R",        # optional: IMD validation (needs IMD files)
  "03_descriptives_diagnostics.R",
  "04_models.R",                 # ~1-2 h
  "05_standardisation.R",
  "06_bootstrap.R",              # ~2-3 h for 500 replicates; resumable
  "07_lisa.R",
  "08_figures_basic.R")
rs <- file.path(R.home("bin"), "Rscript")
dir.create("output_revision", showWarnings = FALSE)
for (s in steps) {
  f <- file.path("code", "B_revision", s)
  cat("\n=====", s, "=====\n")
  t0 <- Sys.time()
  st <- system2(rs, shQuote(f), stdout = file.path("output_revision", sub("\\.R$", ".log", s)),
                stderr = file.path("output_revision", sub("\\.R$", ".log", s)))
  cat(if (st == 0) "OK" else paste("FAILED (exit", st, ") - see output_revision/", sub("\\.R$", ".log", s)),
      "|", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
  if (st != 0 && !s %in% c("02a_download_humidity_rain.R", "02c_timezone_check.R", "02d_imd_validation.R")) break
}
suppressPackageStartupMessages(for (p in c("data.table", "fixest", "haven", "lubridate", "sf", "terra", "Matrix",
  "Hmisc", "metafor", "spdep", "ggplot2", "ecmwfr"))
  try(library(p, character.only = TRUE), silent = TRUE))
writeLines(capture.output(sessionInfo()), "output_revision/sessionInfo.txt")   # with all packages loaded
